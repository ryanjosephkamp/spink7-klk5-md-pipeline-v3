"""CLI entry point for steered molecular dynamics campaigns."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config import DATA_DIR, SMDConfig
from src.simulate.smd import run_smd_campaign


logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for SMD campaigns."""

    defaults = SMDConfig()
    parser = argparse.ArgumentParser(
        description="Run a constant-velocity SMD campaign from an equilibrated state and serialized system.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--state-xml", type=Path, required=True, help="Serialized equilibrated OpenMM state XML path.")
    parser.add_argument("--system-xml", type=Path, required=True, help="Serialized OpenMM system XML path.")
    parser.add_argument("--pull-group-1", required=True, help="Comma-separated atom indices for the first COM group.")
    parser.add_argument("--pull-group-2", required=True, help="Comma-separated atom indices for the second COM group.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_DIR / "trajectories" / "smd",
        help="Output directory for SMD trajectories and work files.",
    )
    parser.add_argument("--spring-constant-kj-mol-nm2", type=float, default=defaults.spring_constant_kj_mol_nm2)
    parser.add_argument("--pulling-velocity-nm-per-ps", type=float, default=defaults.pulling_velocity_nm_per_ps)
    parser.add_argument("--pull-distance-nm", type=float, default=defaults.pull_distance_nm)
    parser.add_argument("--n-replicates", type=int, default=defaults.n_replicates)
    parser.add_argument("--save-interval-ps", type=float, default=defaults.save_interval_ps)
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Python logging level.",
    )
    return parser


def _parse_group(argument: str, name: str) -> list[int]:
    """Parse a comma-separated atom-index list."""

    try:
        indices = [int(token.strip()) for token in argument.split(",") if token.strip()]
    except ValueError as exc:
        raise ValueError(f"{name} must be a comma-separated list of integers") from exc
    if not indices:
        raise ValueError(f"{name} must not be empty")
    return indices


def main(argv: list[str] | None = None) -> int:
    """Run the SMD CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(name)s: %(message)s")

    config = SMDConfig(
        spring_constant_kj_mol_nm2=args.spring_constant_kj_mol_nm2,
        pulling_velocity_nm_per_ps=args.pulling_velocity_nm_per_ps,
        pull_distance_nm=args.pull_distance_nm,
        n_replicates=args.n_replicates,
        save_interval_ps=args.save_interval_ps,
    )
    pull_group_1 = _parse_group(args.pull_group_1, "pull_group_1")
    pull_group_2 = _parse_group(args.pull_group_2, "pull_group_2")

    results = run_smd_campaign(
        Path(args.state_xml),
        Path(args.system_xml),
        config,
        pull_group_1,
        pull_group_2,
        Path(args.output_dir),
    )
    logger.info("Completed %d SMD replicates", len(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())