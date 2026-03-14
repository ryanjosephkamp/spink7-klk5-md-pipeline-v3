"""CLI entry point for production molecular dynamics."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import openmm
from openmm import XmlSerializer, unit
from openmm.app import PDBFile, Simulation

from src.config import DATA_DIR, ProductionConfig
from src.simulate.production import run_production


logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for production MD."""

    defaults = ProductionConfig()
    parser = argparse.ArgumentParser(
        description="Resume from an equilibrated state and run streamed production molecular dynamics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--topology-pdb", type=Path, required=True, help="Solvated topology PDB path.")
    parser.add_argument("--system-xml", type=Path, required=True, help="Serialized OpenMM system XML path.")
    parser.add_argument("--state-xml", type=Path, required=True, help="Serialized OpenMM state XML path.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_DIR / "trajectories" / "production",
        help="Output directory for production outputs.",
    )
    parser.add_argument("--duration-ns", type=float, default=defaults.duration_ns)
    parser.add_argument("--temperature-k", type=float, default=defaults.temperature_k)
    parser.add_argument("--friction-per-ps", type=float, default=defaults.friction_per_ps)
    parser.add_argument("--timestep-ps", type=float, default=defaults.timestep_ps)
    parser.add_argument("--pressure-atm", type=float, default=defaults.pressure_atm)
    parser.add_argument("--save-interval-ps", type=float, default=defaults.save_interval_ps)
    parser.add_argument("--checkpoint-interval-ps", type=float, default=defaults.checkpoint_interval_ps)
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Python logging level.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the production CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(name)s: %(message)s")

    config = ProductionConfig(
        duration_ns=args.duration_ns,
        temperature_k=args.temperature_k,
        friction_per_ps=args.friction_per_ps,
        timestep_ps=args.timestep_ps,
        pressure_atm=args.pressure_atm,
        save_interval_ps=args.save_interval_ps,
        checkpoint_interval_ps=args.checkpoint_interval_ps,
    )

    topology_pdb = PDBFile(str(args.topology_pdb))
    system = XmlSerializer.deserialize(Path(args.system_xml).read_text(encoding="utf-8"))
    state = XmlSerializer.deserialize(Path(args.state_xml).read_text(encoding="utf-8"))
    integrator = openmm.LangevinMiddleIntegrator(
        config.temperature_k * unit.kelvin,
        config.friction_per_ps / unit.picosecond,
        config.timestep_ps * unit.picoseconds,
    )
    simulation = Simulation(topology_pdb.topology, system, integrator, openmm.Platform.getPlatformByName("CPU"))
    simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    simulation.context.setPositions(state.getPositions())
    simulation.context.setVelocities(state.getVelocities())

    result = run_production(simulation, config, Path(args.output_dir))
    logger.info("Production trajectory: %s", result["trajectory_path"])
    logger.info("Energy timeseries: %s", result["energy_timeseries_path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())