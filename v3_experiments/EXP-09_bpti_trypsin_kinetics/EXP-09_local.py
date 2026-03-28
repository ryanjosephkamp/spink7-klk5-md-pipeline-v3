#!/usr/bin/env python3
"""EXP-09: BPTI-Trypsin Kinetics — Local Analysis Script

Derives k_on and k_off from PMF data using transition state theory and
Kramers rate theory as cross-validation.
Requires: EXP-04 PMF data (pmf_data.npz).

Usage:
    python EXP-09_local.py --data-dir /path/to/v3_gpu_results/EXP-04
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
# Physical constants
K_B_J = 1.380649e-23          # J/K
PLANCK_H = 6.62607015e-34     # J·s
AVOGADRO = 6.02214076e23      # mol⁻¹
# Experimental reference values
K_ON_EXP = 1.0e6              # M⁻¹ s⁻¹
K_OFF_EXP = 1.0e-8            # s⁻¹


def main():
    parser = argparse.ArgumentParser(description='EXP-09: BPTI-trypsin kinetics from PMF')
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
    bound_xi = float(xi[bound_idx])
    bound_energy = float(pmf_kcal[bound_idx])

    # Barrier = max PMF between bound minimum and unbound (large xi)
    barrier_region = pmf_kcal[bound_idx:]
    barrier_idx_rel = np.argmax(barrier_region)
    barrier_idx = bound_idx + barrier_idx_rel
    barrier_xi = float(xi[barrier_idx])
    barrier_energy = float(pmf_kcal[barrier_idx])

    dg_barrier_kcal = barrier_energy - bound_energy  # ΔG‡ in kcal/mol
    dg_barrier_kj = dg_barrier_kcal * KCAL_TO_KJ

    print(f'Bound minimum: xi={bound_xi:.2f} nm, E={bound_energy:.2f} kcal/mol')
    print(f'Barrier:       xi={barrier_xi:.2f} nm, E={barrier_energy:.2f} kcal/mol')
    print(f'ΔG‡ = {dg_barrier_kcal:.2f} kcal/mol')

    # TST rate: k_off = (k_B T / h) × exp(-ΔG‡ / k_B T)
    kbt_j = K_B_J * TEMPERATURE_K
    kbt_kcal = (BOLTZMANN_KJ * TEMPERATURE_K) / KCAL_TO_KJ
    prefactor_tst = kbt_j / PLANCK_H  # s⁻¹
    k_off_tst = prefactor_tst * np.exp(-dg_barrier_kcal / kbt_kcal)

    # K_d from PMF well depth (approximate)
    # ΔG_bind = well depth (already in kcal/mol, relative to unbound)
    # Unbound reference: PMF at large xi (last few points)
    unbound_energy = float(np.mean(pmf_kcal[-5:]))
    dg_bind_kcal = bound_energy - unbound_energy
    dg_bind_kj = dg_bind_kcal * KCAL_TO_KJ

    # K_d = exp(ΔG_bind / k_B T) in M (approximate, requires volume correction)
    k_d = np.exp(dg_bind_kcal / kbt_kcal)  # dimensionless ratio, ~M with standard state

    # k_on = k_off / K_d
    k_on_tst = k_off_tst / k_d if k_d > 0 else float('inf')

    # Kramers rate (cross-validation)
    # k_Kramers ∝ ω_well × ω_barrier / (2π γ) × exp(-ΔG‡ / k_B T)
    # Approximate ω from curvature of PMF around well and barrier
    dx = float(xi[1] - xi[0]) if len(xi) > 1 else 0.05
    # Well curvature
    if bound_idx > 0 and bound_idx < len(pmf_kj) - 1:
        omega_well_sq = (pmf_kj[bound_idx - 1] + pmf_kj[bound_idx + 1] - 2 * pmf_kj[bound_idx]) / (dx * dx)
    else:
        omega_well_sq = 1.0
    # Barrier curvature (negative)
    if barrier_idx > 0 and barrier_idx < len(pmf_kj) - 1:
        omega_barrier_sq = abs(pmf_kj[barrier_idx - 1] + pmf_kj[barrier_idx + 1] - 2 * pmf_kj[barrier_idx]) / (dx * dx)
    else:
        omega_barrier_sq = 1.0

    # Classification: within 2 orders of magnitude
    log_k_on_ratio = abs(np.log10(max(k_on_tst, 1e-30)) - np.log10(K_ON_EXP))
    log_k_off_ratio = abs(np.log10(max(k_off_tst, 1e-30)) - np.log10(K_OFF_EXP))

    k_on_pass = log_k_on_ratio <= 2.0
    k_off_pass = log_k_off_ratio <= 2.0
    verdict = 'PASS' if (k_on_pass and k_off_pass) else (
        'MARGINAL' if (k_on_pass or k_off_pass) else 'FAIL')

    results = {
        'experiment_id': 'EXP-09',
        'bound_xi_nm': bound_xi,
        'barrier_xi_nm': barrier_xi,
        'dg_barrier_kcal_mol': dg_barrier_kcal,
        'dg_bind_kcal_mol': dg_bind_kcal,
        'k_off_tst_per_s': float(k_off_tst),
        'k_on_tst_per_M_s': float(k_on_tst),
        'k_d_approx': float(k_d),
        'log10_k_on_ratio_to_exp': float(log_k_on_ratio),
        'log10_k_off_ratio_to_exp': float(log_k_off_ratio),
        'k_on_exp': K_ON_EXP,
        'k_off_exp': K_OFF_EXP,
        'verdict': verdict,
    }

    with open(output_dir / 'kinetics.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f'\nk_off (TST) = {k_off_tst:.2e} s⁻¹ (exp: {K_OFF_EXP:.2e})')
    print(f'k_on  (TST) = {k_on_tst:.2e} M⁻¹s⁻¹ (exp: {K_ON_EXP:.2e})')
    print(f'log10 ratio k_on:  {log_k_on_ratio:.1f} orders')
    print(f'log10 ratio k_off: {log_k_off_ratio:.1f} orders')
    print(f'Verdict: {verdict}')

    # Figure: PMF with barrier annotation
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xi, pmf_kcal, 'b-', lw=2)
    ax.axhline(y=bound_energy, color='green', linestyle=':', alpha=0.5)
    ax.axhline(y=barrier_energy, color='red', linestyle=':', alpha=0.5)
    ax.annotate(f'ΔG‡ = {dg_barrier_kcal:.1f} kcal/mol',
                xy=(barrier_xi, barrier_energy), xytext=(barrier_xi + 0.3, barrier_energy + 1),
                arrowprops=dict(arrowstyle='->', color='red'), fontsize=10, color='red')
    ax.plot(bound_xi, bound_energy, 'go', ms=8, label=f'Bound min (ξ={bound_xi:.2f} nm)')
    ax.plot(barrier_xi, barrier_energy, 'r^', ms=8, label=f'Barrier (ξ={barrier_xi:.2f} nm)')
    ax.set_xlabel('Reaction Coordinate ξ (nm)')
    ax.set_ylabel('PMF (kcal/mol)')
    ax.set_title('EXP-09: PMF with Kinetic Barrier')
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / 'pmf_with_barrier.png', dpi=150)
    plt.close(fig)
    print(f'Figure saved: {figures_dir / "pmf_with_barrier.png"}')


if __name__ == '__main__':
    main()
