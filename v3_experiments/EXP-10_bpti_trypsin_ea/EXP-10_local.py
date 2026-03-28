#!/usr/bin/env python3
"""EXP-10: BPTI-Trypsin Activation Energy — Local Analysis Script

Extracts barrier height from PMF and computes Ea with Eyring correction.
Requires: EXP-04 PMF data (pmf_data.npz).

Usage:
    python EXP-10_local.py --data-dir /path/to/v3_gpu_results/EXP-04
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import BOLTZMANN_KJ, KCAL_TO_KJ

TEMPERATURE_K = 310.0
EA_EXPERIMENTAL = 10.5  # kcal/mol
EA_PASS_LO = 4.6
EA_PASS_HI = 16.4


def main():
    parser = argparse.ArgumentParser(description='EXP-10: BPTI-trypsin activation energy from PMF')
    parser.add_argument('--data-dir', type=Path, required=True,
                        help='Path to EXP-04 GPU results directory')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Output directory (default: ./outputs)')
    args = parser.parse_args()

    output_dir = args.output_dir or Path(__file__).parent / 'outputs'
    figures_dir = Path(__file__).parent / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Load PMF data
    pmf_path = args.data_dir / 'analysis' / 'pmf_data.npz'
    if not pmf_path.exists():
        pmf_path = args.data_dir / 'pmf_data.npz'
    if not pmf_path.exists():
        print(f'ERROR: PMF data not found in {args.data_dir}')
        sys.exit(1)

    print(f'Loading PMF: {pmf_path}')
    data = np.load(str(pmf_path))
    xi = data['xi_bins'] if 'xi_bins' in data else data['xi']
    pmf_kj = data['pmf_kj_mol'] if 'pmf_kj_mol' in data else data['pmf']
    pmf_kcal = pmf_kj / KCAL_TO_KJ

    # Identify bound minimum and barrier
    bound_idx = np.argmin(pmf_kcal)
    bound_energy = float(pmf_kcal[bound_idx])

    barrier_region = pmf_kcal[bound_idx:]
    barrier_idx_rel = np.argmax(barrier_region)
    barrier_idx = bound_idx + barrier_idx_rel
    barrier_energy = float(pmf_kcal[barrier_idx])

    dg_barrier_kcal = barrier_energy - bound_energy  # ΔG‡

    # Eyring correction: Ea ≈ ΔG‡ + k_B T
    kbt_kcal = (BOLTZMANN_KJ * TEMPERATURE_K) / KCAL_TO_KJ
    ea_kcal = dg_barrier_kcal + kbt_kcal

    print(f'ΔG‡ = {dg_barrier_kcal:.2f} kcal/mol')
    print(f'k_B T = {kbt_kcal:.3f} kcal/mol')
    print(f'Ea = ΔG‡ + k_B T = {ea_kcal:.2f} kcal/mol')
    print(f'Experimental Ea = {EA_EXPERIMENTAL} kcal/mol')
    print(f'Pass range: [{EA_PASS_LO}, {EA_PASS_HI}] kcal/mol')

    # Classification
    in_range = EA_PASS_LO <= ea_kcal <= EA_PASS_HI
    verdict = 'PASS' if in_range else 'FAIL'

    results = {
        'experiment_id': 'EXP-10',
        'dg_barrier_kcal_mol': dg_barrier_kcal,
        'kbt_kcal_mol': kbt_kcal,
        'ea_kcal_mol': ea_kcal,
        'ea_experimental_kcal_mol': EA_EXPERIMENTAL,
        'pass_range_kcal_mol': [EA_PASS_LO, EA_PASS_HI],
        'verdict': verdict,
    }

    with open(output_dir / 'activation_energy.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f'Verdict: {verdict}')

    # Figure: PMF with barrier annotation
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xi, pmf_kcal, 'b-', lw=2)
    ax.fill_between([float(xi[bound_idx]), float(xi[barrier_idx])],
                    bound_energy, barrier_energy, alpha=0.15, color='red')
    ax.annotate(f'Ea ≈ {ea_kcal:.1f} kcal/mol\n(ΔG‡ + k_BT)',
                xy=(float(xi[barrier_idx]), barrier_energy),
                xytext=(float(xi[barrier_idx]) + 0.3, barrier_energy + 1),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10, color='red', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.axhline(y=bound_energy, color='green', linestyle=':', alpha=0.4)
    ax.set_xlabel('Reaction Coordinate ξ (nm)')
    ax.set_ylabel('PMF (kcal/mol)')
    ax.set_title(f'EXP-10: Activation Energy (Ea = {ea_kcal:.1f} kcal/mol, exp = {EA_EXPERIMENTAL})')
    fig.tight_layout()
    fig.savefig(figures_dir / 'pmf_barrier_annotation.png', dpi=150)
    plt.close(fig)
    print(f'Figure saved: {figures_dir / "pmf_barrier_annotation.png"}')


if __name__ == '__main__':
    main()
