"""CLI entry point for NVT and NPT equilibration."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import openmm
from openmm import XmlSerializer, unit
from openmm.app import PDBFile, Simulation

from src.config import DATA_DIR, EquilibrationConfig, MinimizationConfig
from src.config import load_config
from src.physics.restraints import create_positional_restraints
from src.simulate.equilibrate import run_npt, run_nvt
from src.simulate.minimizer import minimize_energy


logger = logging.getLogger(__name__)

_SOLVENT_RESIDUES = {"HOH", "WAT", "NA", "CL"}


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for equilibration."""

    equilibration_defaults = EquilibrationConfig()
    minimization_defaults = MinimizationConfig()
    parser = argparse.ArgumentParser(
        description="Load a solvated system, apply heavy-atom restraints, minimize, and run NVT then NPT equilibration.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--topology-pdb", type=Path, required=True, help="Solvated topology PDB path.")
    parser.add_argument("--system-xml", type=Path, required=True, help="Serialized solvated OpenMM system XML path.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_DIR / "trajectories" / "equilibration",
        help="Output directory for equilibration artifacts.",
    )
    parser.add_argument("--minimize", action="store_true", help="Run energy minimization before equilibration.")
    parser.add_argument(
        "--max-minimization-iterations",
        type=int,
        default=minimization_defaults.max_iterations,
        help="Maximum minimization iterations if --minimize is set.",
    )
    parser.add_argument(
        "--minimization-tolerance-kj-mol-nm",
        type=float,
        default=minimization_defaults.tolerance_kj_mol_nm,
        help="Minimization tolerance in kJ/mol/nm.",
    )
    parser.add_argument("--nvt-duration-ps", type=float, default=equilibration_defaults.nvt_duration_ps)
    parser.add_argument("--npt-duration-ps", type=float, default=equilibration_defaults.npt_duration_ps)
    parser.add_argument("--temperature-k", type=float, default=equilibration_defaults.temperature_k)
    parser.add_argument("--friction-per-ps", type=float, default=equilibration_defaults.friction_per_ps)
    parser.add_argument("--timestep-ps", type=float, default=equilibration_defaults.timestep_ps)
    parser.add_argument("--pressure-atm", type=float, default=equilibration_defaults.pressure_atm)
    parser.add_argument("--barostat-interval", type=int, default=equilibration_defaults.barostat_interval)
    parser.add_argument(
        "--restraint-k-kj-mol-nm2",
        type=float,
        default=equilibration_defaults.restraint_k_kj_mol_nm2,
        help="Heavy-atom restraint force constant.",
    )
    parser.add_argument("--save-interval-ps", type=float, default=equilibration_defaults.save_interval_ps)
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Python logging level.",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to a YAML configuration file. Missing fields fall back to defaults.",
    )
    parser.add_argument(
        "--platform",
        default=None,
        choices=["CUDA", "OpenCL", "CPU"],
        help="OpenMM platform. Auto-detected (CUDA->OpenCL->CPU) if omitted.",
    )
    return parser


def _restrained_heavy_atoms(topology, positions) -> tuple[list[int], np.ndarray]:
    """Return heavy-atom indices and reference positions for positional restraints."""

    restrained_indices: list[int] = []
    reference_positions: list[list[float]] = []
    for atom_index, atom in enumerate(topology.atoms()):
        if atom.residue.name.upper() in _SOLVENT_RESIDUES:
            continue
        if atom.element is None or atom.element.symbol == "H":
            continue
        position_nm = positions[atom_index].value_in_unit(unit.nanometer)
        restrained_indices.append(atom_index)
        reference_positions.append([float(position_nm[0]), float(position_nm[1]), float(position_nm[2])])
    return restrained_indices, np.asarray(reference_positions, dtype=float)


def main(argv: list[str] | None = None) -> int:
    """Run the equilibration CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(name)s: %(message)s")

    topology_pdb = PDBFile(str(args.topology_pdb))
    system = XmlSerializer.deserialize(Path(args.system_xml).read_text(encoding="utf-8"))
    restrained_indices, reference_positions = _restrained_heavy_atoms(topology_pdb.topology, topology_pdb.positions)
    create_positional_restraints(system, restrained_indices, reference_positions, args.restraint_k_kj_mol_nm2)

    if args.config is not None:
        configs = load_config(args.config)
        equilibration_config = configs["equilibration"]
    else:
        equilibration_config = EquilibrationConfig(
            nvt_duration_ps=args.nvt_duration_ps,
            npt_duration_ps=args.npt_duration_ps,
            temperature_k=args.temperature_k,
            friction_per_ps=args.friction_per_ps,
            timestep_ps=args.timestep_ps,
            pressure_atm=args.pressure_atm,
            barostat_interval=args.barostat_interval,
            restraint_k_kj_mol_nm2=args.restraint_k_kj_mol_nm2,
            save_interval_ps=args.save_interval_ps,
        )
    integrator = openmm.LangevinMiddleIntegrator(
        equilibration_config.temperature_k * unit.kelvin,
        equilibration_config.friction_per_ps / unit.picosecond,
        equilibration_config.timestep_ps * unit.picoseconds,
    )
    from src.simulate.platform import select_platform
    simulation = Simulation(topology_pdb.topology, system, integrator, select_platform(args.platform))
    logger.info("Simulation platform: %s", simulation.context.getPlatform().getName())
    simulation.context.setPositions(topology_pdb.positions)

    if args.minimize:
        minimization_config = MinimizationConfig(
            max_iterations=args.max_minimization_iterations,
            tolerance_kj_mol_nm=args.minimization_tolerance_kj_mol_nm,
        )
        minimize_energy(simulation, minimization_config)

    output_dir = Path(args.output_dir)
    nvt_result = run_nvt(simulation, equilibration_config, output_dir / "nvt")
    npt_result = run_npt(simulation, equilibration_config, output_dir / "npt")

    logger.info("NVT trajectory: %s", nvt_result["trajectory_path"])
    logger.info("NVT final state: %s", nvt_result["final_state_path"])
    logger.info("NPT trajectory: %s", npt_result["trajectory_path"])
    logger.info("NPT final state: %s", npt_result["final_state_path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())