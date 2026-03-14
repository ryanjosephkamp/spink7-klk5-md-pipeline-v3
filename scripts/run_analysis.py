"""CLI entry point for post-processing and figure generation."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from src.config import DATA_DIR, ProductionConfig, WHAMConfig
from src.analyze.jarzynski import evaluate_convergence, jarzynski_free_energy
from src.analyze.structural import compute_radius_of_gyration, compute_rmsd, compute_rmsf, compute_sasa
from src.analyze.trajectory import load_trajectory
from src.analyze.wham import bootstrap_pmf_uncertainty, solve_wham
from src.visualization.plot_pmf import plot_pmf
from src.visualization.plot_timeseries import plot_energy_timeseries, plot_rmsd_timeseries, plot_temperature_timeseries


logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for analysis workflows."""

    defaults = ProductionConfig()
    wham_defaults = WHAMConfig()
    parser = argparse.ArgumentParser(
        description="Run structural analysis, free-energy estimation, and figure generation workflows.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Python logging level.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    structural = subparsers.add_parser(
        "structural",
        help="Compute RMSD, RMSF, radius of gyration, and SASA from a trajectory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    structural.add_argument("--trajectory", type=Path, required=True, help="Trajectory path (.dcd/.xtc/etc.).")
    structural.add_argument("--topology", type=Path, required=True, help="Topology path for MDTraj loading.")
    structural.add_argument("--stride", type=int, default=1, help="Frame stride for MDTraj loading.")
    structural.add_argument("--rmsd-selection", default="backbone", help="Atom selection used for RMSD.")
    structural.add_argument("--rmsf-selection", default="name CA", help="Atom selection used for RMSF.")
    structural.add_argument("--probe-radius-nm", type=float, default=0.14, help="Probe radius for SASA.")
    structural.add_argument(
        "--output",
        type=Path,
        default=DATA_DIR / "analysis" / "structural_metrics.npz",
        help="Output NPZ path for structural metrics.",
    )

    jarzynski = subparsers.add_parser(
        "jarzynski",
        help="Compute Jarzynski free-energy estimates from total-work samples.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    jarzynski.add_argument("--work-values", type=Path, required=True, help="CSV or NPY file of total work values.")
    jarzynski.add_argument("--temperature-k", type=float, default=defaults.temperature_k)
    jarzynski.add_argument("--n-subsets", type=int, default=10)
    jarzynski.add_argument(
        "--output",
        type=Path,
        default=DATA_DIR / "analysis" / "jarzynski_summary.npz",
        help="Output NPZ path for Jarzynski results.",
    )

    wham = subparsers.add_parser(
        "wham",
        help="Solve WHAM from umbrella xi timeseries and optionally bootstrap uncertainty.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    wham.add_argument("--xi-files", type=Path, nargs="+", required=True, help="Per-window xi timeseries .npy files.")
    wham.add_argument(
        "--window-centers",
        type=float,
        nargs="+",
        required=True,
        help="Window centers in nm, one per xi file.",
    )
    wham.add_argument(
        "--spring-constants",
        type=float,
        nargs="+",
        required=True,
        help="Spring constants in kJ/mol/nm^2; supply one value or one per window.",
    )
    wham.add_argument("--temperature-k", type=float, default=defaults.temperature_k)
    wham.add_argument("--tolerance", type=float, default=wham_defaults.tolerance)
    wham.add_argument("--max-iterations", type=int, default=wham_defaults.max_iterations)
    wham.add_argument("--n-bootstrap", type=int, default=wham_defaults.n_bootstrap)
    wham.add_argument("--histogram-bins", type=int, default=wham_defaults.histogram_bins)
    wham.add_argument("--bootstrap", action="store_true", help="Also compute bootstrap PMF uncertainty.")
    wham.add_argument(
        "--output",
        type=Path,
        default=DATA_DIR / "analysis" / "pmf" / "wham_pmf.npz",
        help="Output NPZ path for WHAM results.",
    )

    plot_pmf_parser = subparsers.add_parser(
        "plot-pmf",
        help="Render a PMF figure from a saved NPZ file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plot_pmf_parser.add_argument("--pmf-npz", type=Path, required=True, help="NPZ file containing xi_bins and pmf arrays.")
    plot_pmf_parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "figures" / "pmf.png",
        help="Output figure path.",
    )

    plot_timeseries = subparsers.add_parser(
        "plot-timeseries",
        help="Render energy, temperature, and RMSD time-series figures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plot_timeseries.add_argument("--energy-csv", type=Path, help="CSV with columns time_ps,potential_kj_mol,kinetic_kj_mol.")
    plot_timeseries.add_argument("--temperature-csv", type=Path, help="CSV with columns time_ps,temperature_k.")
    plot_timeseries.add_argument("--rmsd-npy", type=Path, help="NPY file of RMSD values in nm.")
    plot_timeseries.add_argument(
        "--rmsd-time-step-ns",
        type=float,
        default=0.01,
        help="Frame spacing in ns when plotting RMSD from an NPY array.",
    )
    plot_timeseries.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "figures",
        help="Output directory for generated figures.",
    )

    return parser


def _load_work_values(work_path: Path) -> np.ndarray:
    """Load total-work values from NPY or CSV."""

    if work_path.suffix.lower() == ".npy":
        values = np.load(work_path)
    else:
        values = np.loadtxt(work_path, delimiter=",", ndmin=1)
    return np.asarray(values, dtype=float).reshape(-1)


def _broadcast_spring_constants(values: list[float], n_windows: int) -> np.ndarray:
    """Expand a single spring constant to all windows when needed."""

    if len(values) == 1:
        return np.full(n_windows, float(values[0]), dtype=float)
    return np.asarray(values, dtype=float)


def _run_structural(args: argparse.Namespace) -> None:
    """Execute structural trajectory analysis."""

    trajectory = load_trajectory(args.trajectory, args.topology, stride=args.stride)
    rmsd_nm = compute_rmsd(trajectory, trajectory[0], atom_selection=args.rmsd_selection)
    rmsf_nm = compute_rmsf(trajectory, atom_selection=args.rmsf_selection)
    rg_nm = compute_radius_of_gyration(trajectory)
    sasa_nm2 = compute_sasa(trajectory, probe_radius_nm=args.probe_radius_nm)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output_path, rmsd_nm=rmsd_nm, rmsf_nm=rmsf_nm, radius_of_gyration_nm=rg_nm, sasa_nm2=sasa_nm2)
    logger.info("Structural metrics written to %s", output_path)


def _run_jarzynski(args: argparse.Namespace) -> None:
    """Execute Jarzynski post-processing."""

    work_values = _load_work_values(Path(args.work_values))
    estimates = jarzynski_free_energy(work_values, args.temperature_k)
    convergence = evaluate_convergence(work_values, args.temperature_k, n_subsets=args.n_subsets)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output_path, **estimates, **convergence)
    logger.info("Jarzynski results written to %s", output_path)


def _run_wham(args: argparse.Namespace) -> None:
    """Execute WHAM analysis and optional bootstrapping."""

    xi_timeseries_list = [np.asarray(np.load(path), dtype=float) for path in args.xi_files]
    window_centers = np.asarray(args.window_centers, dtype=float)
    spring_constants = _broadcast_spring_constants(args.spring_constants, len(xi_timeseries_list))
    config = WHAMConfig(
        tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        n_bootstrap=args.n_bootstrap,
        histogram_bins=args.histogram_bins,
    )
    result = solve_wham(xi_timeseries_list, window_centers, spring_constants, args.temperature_k, config)
    arrays: dict[str, np.ndarray | int | bool] = dict(result)
    if args.bootstrap:
        arrays.update(bootstrap_pmf_uncertainty(xi_timeseries_list, window_centers, spring_constants, args.temperature_k, config))
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output_path, **arrays)
    logger.info("WHAM results written to %s", output_path)


def _run_plot_pmf(args: argparse.Namespace) -> None:
    """Render a PMF figure from stored arrays."""

    arrays = np.load(args.pmf_npz)
    pmf_std = arrays["pmf_std"] if "pmf_std" in arrays.files else None
    figure = plot_pmf(arrays["xi_bins"], arrays["pmf_kcal_mol"], pmf_std_kcal_mol=pmf_std, output_path=args.output)
    logger.info("PMF figure written to %s", args.output)
    figure.clf()


def _run_plot_timeseries(args: argparse.Namespace) -> None:
    """Render time-series figures for available observables."""

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.energy_csv is not None:
        energy = np.loadtxt(args.energy_csv, delimiter=",", skiprows=1)
        figure = plot_energy_timeseries(energy[:, 0], energy[:, 1], energy[:, 2], output_dir / "energy_timeseries.png")
        figure.clf()
        logger.info("Energy figure written to %s", output_dir / "energy_timeseries.png")
    if args.temperature_csv is not None:
        temperature = np.loadtxt(args.temperature_csv, delimiter=",", skiprows=1)
        figure = plot_temperature_timeseries(temperature[:, 0], temperature[:, 1], output_dir / "temperature_timeseries.png")
        figure.clf()
        logger.info("Temperature figure written to %s", output_dir / "temperature_timeseries.png")
    if args.rmsd_npy is not None:
        rmsd_nm = np.asarray(np.load(args.rmsd_npy), dtype=float)
        time_ns = np.arange(rmsd_nm.size, dtype=float) * float(args.rmsd_time_step_ns)
        figure = plot_rmsd_timeseries(time_ns, rmsd_nm, output_dir / "rmsd_timeseries.png")
        figure.clf()
        logger.info("RMSD figure written to %s", output_dir / "rmsd_timeseries.png")


def main(argv: list[str] | None = None) -> int:
    """Run the analysis CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(name)s: %(message)s")

    if args.command == "structural":
        _run_structural(args)
    elif args.command == "jarzynski":
        _run_jarzynski(args)
    elif args.command == "wham":
        _run_wham(args)
    elif args.command == "plot-pmf":
        _run_plot_pmf(args)
    elif args.command == "plot-timeseries":
        _run_plot_timeseries(args)
    else:
        raise ValueError(f"Unsupported analysis command: {args.command}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())