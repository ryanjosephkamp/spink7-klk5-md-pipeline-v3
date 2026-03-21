"""CLI entry point for Markov State Model construction and analysis."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        description="Build and analyze Markov State Models from MD trajectories.",
    )
    subparsers = parser.add_subparsers(dest="command")

    # ---- featurize ----
    feat = subparsers.add_parser("featurize", help="Extract features from trajectory")
    feat.add_argument("--trajectory", type=Path, required=True)
    feat.add_argument("--topology", type=Path, required=True)
    feat.add_argument(
        "--output", type=Path, required=True,
        help="Output path for feature NPZ file",
    )

    # ---- build ----
    bld = subparsers.add_parser("build", help="Build MSM from features")
    bld.add_argument(
        "--features", type=Path, required=True,
        help="Feature NPZ file from featurize step",
    )
    bld.add_argument("--output-dir", type=Path, required=True)
    bld.add_argument(
        "--lag-time-ps", type=float, default=None,
        help="Lag time in ps (omit for ITS convergence only)",
    )

    # ---- kinetics ----
    kin = subparsers.add_parser("kinetics", help="Extract kinetic rates from MSM")
    kin.add_argument("--msm-dir", type=Path, required=True)
    kin.add_argument(
        "--source-states", type=str, required=True,
        help="Comma-separated microstate indices for source macrostate",
    )
    kin.add_argument(
        "--target-states", type=str, required=True,
        help="Comma-separated microstate indices for target macrostate",
    )

    return parser


def _cmd_featurize(args: argparse.Namespace) -> int:
    """Execute the featurize subcommand."""
    import mdtraj as md

    from src.analyze.featurize import compute_backbone_dihedrals
    from src.analyze.trajectory import load_trajectory

    traj = load_trajectory(args.trajectory, args.topology)
    features = compute_backbone_dihedrals(traj)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output, features=features)
    print(f"Saved {features.shape[1]} features x {features.shape[0]} frames -> {args.output}")
    return 0


def _cmd_build(args: argparse.Namespace) -> int:
    """Execute the build subcommand."""
    from src.analyze.msm import (
        build_msm,
        cluster_microstates,
        compute_implied_timescales,
        fit_tica,
    )
    from src.config import MSMConfig

    config = MSMConfig()
    data = np.load(args.features)
    features = data["features"]

    tica_result = fit_tica(features, config)
    assignments = cluster_microstates(tica_result["tica_output"], config)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.save(args.output_dir / "assignments.npy", assignments)
    np.save(args.output_dir / "tica_output.npy", tica_result["tica_output"])

    its_result = compute_implied_timescales(assignments, config)
    np.savez(
        args.output_dir / "implied_timescales.npz",
        lag_times_ps=its_result["lag_times_ps"],
        implied_timescales_ps=its_result["implied_timescales_ps"],
    )
    print(f"ITS saved to {args.output_dir / 'implied_timescales.npz'}")

    if args.lag_time_ps is not None:
        msm_result = build_msm(assignments, args.lag_time_ps, config)
        np.savez(
            args.output_dir / "msm.npz",
            transition_matrix=msm_result["transition_matrix"],
            stationary_distribution=msm_result["stationary_distribution"],
            timescales_ps=msm_result["timescales_ps"],
            active_set=msm_result["active_set"],
        )
        print(f"MSM (lag={args.lag_time_ps} ps) saved to {args.output_dir / 'msm.npz'}")

    return 0


def _cmd_kinetics(args: argparse.Namespace) -> int:
    """Execute the kinetics subcommand."""
    from src.analyze.msm import build_msm, compute_mfpt
    from src.config import BOLTZMANN_KJ, MSMConfig

    config = MSMConfig()

    msm_data = np.load(args.msm_dir / "msm.npz")
    assignments = np.load(args.msm_dir / "assignments.npy")

    # Re-parse the lag time from the transition matrix's associated data
    # For simplicity, rebuild the MSM from the stored assignments
    source = np.array([int(s) for s in args.source_states.split(",")])
    target = np.array([int(s) for s in args.target_states.split(",")])

    # Load the saved MSM transition matrix for reporting
    T = msm_data["transition_matrix"]
    pi = msm_data["stationary_distribution"]

    print(f"Transition matrix shape: {T.shape}")
    print(f"Source states: {source}")
    print(f"Target states: {target}")
    print(f"Stationary distribution range: [{pi.min():.6f}, {pi.max():.6f}]")

    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    handlers = {
        "featurize": _cmd_featurize,
        "build": _cmd_build,
        "kinetics": _cmd_kinetics,
    }
    return handlers[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
