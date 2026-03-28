#!/usr/bin/env python3
"""EXP-25: BPTI Conformational Variability — Local Analysis Script

Computes Cα RMSD fluctuation from apo BPTI production trajectory.
Requires: EXP-24 production trajectory and topology.

Usage:
    python EXP-25_local.py --data-dir /path/to/v3_gpu_results/EXP-24
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.analyze.structural import compute_rmsd

RMSD_PASS_LO = 0.025  # nm (0.25 Å)
RMSD_PASS_HI = 0.055  # nm (0.55 Å)


def main():
    parser = argparse.ArgumentParser(description='EXP-25: BPTI conformational variability')
    parser.add_argument('--data-dir', type=Path, required=True,
                        help='Path to EXP-24 GPU results directory')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Output directory (default: ./outputs)')
    args = parser.parse_args()

    output_dir = args.output_dir or Path(__file__).parent / 'outputs'
    figures_dir = Path(__file__).parent / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Load trajectory
    topo_path = args.data_dir / 'topology_reference.pdb'
    traj_candidates = [args.data_dir / 'production_trajectory.dcd',
                       args.data_dir / 'production_trajectory.xtc']
    traj_path = next((p for p in traj_candidates if p.exists()), None)
    if traj_path is None:
        print(f'ERROR: No trajectory found in {args.data_dir}')
        sys.exit(1)

    print(f'Loading trajectory: {traj_path}')
    traj = md.load(str(traj_path), top=str(topo_path))
    print(f'Loaded {traj.n_frames} frames, {traj.n_atoms} atoms')

    # Compute Cα RMSD to initial frame
    print('Computing Cα RMSD...')
    reference = traj[0]
    ca_indices = traj.topology.select('name CA')
    print(f'{len(ca_indices)} Cα atoms')

    rmsd_nm = compute_rmsd(traj, reference, atom_selection='name CA')
    rmsd_ang = rmsd_nm * 10.0  # convert nm to Å

    mean_rmsd_nm = float(np.mean(rmsd_nm))
    std_rmsd_nm = float(np.std(rmsd_nm))
    mean_rmsd_ang = mean_rmsd_nm * 10.0
    std_rmsd_ang = std_rmsd_nm * 10.0

    # The "RMSD fluctuation" is the standard deviation of the RMSD timeseries
    rmsd_fluctuation_nm = std_rmsd_nm
    rmsd_fluctuation_ang = std_rmsd_ang

    print(f'Mean Cα RMSD = {mean_rmsd_ang:.3f} Å')
    print(f'RMSD fluctuation (std) = {rmsd_fluctuation_ang:.3f} Å')
    print(f'Target: [{RMSD_PASS_LO * 10:.2f}, {RMSD_PASS_HI * 10:.2f}] Å')

    # Classification
    in_range = RMSD_PASS_LO <= rmsd_fluctuation_nm <= RMSD_PASS_HI
    verdict = 'PASS' if in_range else 'FAIL'

    results = {
        'experiment_id': 'EXP-25',
        'n_frames': int(traj.n_frames),
        'n_ca_atoms': int(len(ca_indices)),
        'mean_rmsd_nm': mean_rmsd_nm,
        'mean_rmsd_ang': mean_rmsd_ang,
        'rmsd_fluctuation_nm': rmsd_fluctuation_nm,
        'rmsd_fluctuation_ang': rmsd_fluctuation_ang,
        'pass_range_ang': [RMSD_PASS_LO * 10, RMSD_PASS_HI * 10],
        'verdict': verdict,
    }

    with open(output_dir / 'rmsd_analysis.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f'Verdict: {verdict}')

    # Figure: RMSD timeseries
    time_ns = np.arange(traj.n_frames) * (traj.timestep / 1000.0) if traj.timestep else np.arange(traj.n_frames)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(time_ns, rmsd_ang, 'b-', alpha=0.7, lw=0.5)
    ax1.axhline(y=mean_rmsd_ang, color='red', linestyle='--', label=f'Mean = {mean_rmsd_ang:.3f} Å')
    ax1.set_xlabel('Time (ns)' if traj.timestep else 'Frame')
    ax1.set_ylabel('Cα RMSD (Å)')
    ax1.set_title('EXP-25: Cα RMSD Timeseries')
    ax1.legend()

    ax2.hist(rmsd_ang, bins=50, edgecolor='black', alpha=0.7, density=True)
    ax2.axvline(mean_rmsd_ang, color='red', linestyle='--', label=f'Mean = {mean_rmsd_ang:.3f} Å')
    ax2.set_xlabel('Cα RMSD (Å)')
    ax2.set_ylabel('Density')
    ax2.set_title(f'EXP-25: RMSD Distribution (σ = {rmsd_fluctuation_ang:.3f} Å)')
    ax2.legend()

    fig.tight_layout()
    fig.savefig(figures_dir / 'rmsd_timeseries.png', dpi=150)
    plt.close(fig)
    print(f'Figure saved: {figures_dir / "rmsd_timeseries.png"}')


if __name__ == '__main__':
    main()
