"""CLI entry point for umbrella-sampling campaigns."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.config import DATA_DIR, UmbrellaConfig
from src.simulate.umbrella import generate_window_centers, run_umbrella_campaign


logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for umbrella sampling."""

    defaults = UmbrellaConfig()
    parser = argparse.ArgumentParser(
        description="Run an umbrella-sampling campaign from an equilibrated state and serialized system.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--state-xml", type=Path, required=True, help="Serialized equilibrated OpenMM state XML path.")
    parser.add_argument("--system-xml", type=Path, required=True, help="Serialized OpenMM system XML path.")
    parser.add_argument("--pull-group-1", required=True, help="Comma-separated atom indices for the first COM group.")
    parser.add_argument("--pull-group-2", required=True, help="Comma-separated atom indices for the second COM group.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DATA_DIR / "trajectories" / "umbrella",
        help="Output directory for umbrella trajectories and xi files.",
    )
    parser.add_argument("--xi-min-nm", type=float, default=defaults.xi_min_nm)
    parser.add_argument("--xi-max-nm", type=float, default=defaults.xi_max_nm)
    parser.add_argument("--window-spacing-nm", type=float, default=defaults.window_spacing_nm)
    parser.add_argument("--spring-constant-kj-mol-nm2", type=float, default=defaults.spring_constant_kj_mol_nm2)
    parser.add_argument("--per-window-duration-ns", type=float, default=defaults.per_window_duration_ns)
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
    """Run the umbrella-sampling CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(name)s: %(message)s")

    config = UmbrellaConfig(
        xi_min_nm=args.xi_min_nm,
        xi_max_nm=args.xi_max_nm,
        window_spacing_nm=args.window_spacing_nm,
        spring_constant_kj_mol_nm2=args.spring_constant_kj_mol_nm2,
        per_window_duration_ns=args.per_window_duration_ns,
        save_interval_ps=args.save_interval_ps,
    )
    pull_group_1 = _parse_group(args.pull_group_1, "pull_group_1")
    pull_group_2 = _parse_group(args.pull_group_2, "pull_group_2")

    results = run_umbrella_campaign(
        Path(args.state_xml),
        Path(args.system_xml),
        config,
        pull_group_1,
        pull_group_2,
        Path(args.output_dir),
    )
    logger.info("Completed %d umbrella windows", len(results))
    logger.info("Window centers: %s", generate_window_centers(config))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())