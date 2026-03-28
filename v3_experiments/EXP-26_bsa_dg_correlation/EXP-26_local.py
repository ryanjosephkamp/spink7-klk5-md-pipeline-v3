#!/usr/bin/env python3
"""EXP-26: BSA–ΔG Correlation — Local Analysis Script

Performs linear regression of buried surface area vs. binding free energy
across 6 protein–protein systems.
Requires: results.json from EXP-04, EXP-05, EXP-06, EXP-13, EXP-28, EXP-29.

Usage:
    python EXP-26_local.py --results-dir /path/to/v3_gpu_results
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

# 6 systems with experimental ΔG values (kcal/mol)
SYSTEMS = {
    'EXP-04': {'name': 'BPTI-Trypsin', 'pdb': '2PTC', 'dg_exp': -18.0},
    'EXP-05': {'name': 'PSTI-Chymotrypsinogen', 'pdb': '1TGS', 'dg_exp': -14.7},
    'EXP-06': {'name': 'SPINK1-Trypsin', 'pdb': '1TGS', 'dg_exp': -11.1},
    'EXP-13': {'name': 'Barnase-Barstar', 'pdb': '1BRS', 'dg_exp': -19.0},
    'EXP-28': {'name': 'BPTI-loop–Trypsin', 'pdb': '2PTC', 'dg_exp': None},  # computed
    'EXP-29': {'name': 'SH3-p41', 'pdb': '4EIK', 'dg_exp': -7.99},
}

R2_PASS = 0.5
SLOPE_RANGE = (-0.020, -0.005)  # kcal/(mol·Å²)


def main():
    parser = argparse.ArgumentParser(description='EXP-26: BSA–ΔG correlation')
    parser.add_argument('--results-dir', type=Path, required=True,
                        help='Path to v3_gpu_results directory containing EXP-XX subdirs')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Output directory (default: ./outputs)')
    args = parser.parse_args()

    output_dir = args.output_dir or Path(__file__).parent / 'outputs'
    figures_dir = Path(__file__).parent / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Collect BSA and ΔG from each experiment
    bsa_values = []
    dg_values = []
    labels = []
    missing = []

    for exp_id, info in SYSTEMS.items():
        results_path = args.results_dir / exp_id / 'results.json'
        if not results_path.exists():
            print(f'WARNING: {results_path} not found — skipping {exp_id}')
            missing.append(exp_id)
            continue

        with open(results_path) as f:
            r = json.load(f)

        # Extract BSA (try multiple possible keys)
        bsa = r.get('mean_bsa_ang2', r.get('bsa_ang2', r.get('buried_surface_area_ang2')))
        if bsa is None:
            # Try nm² and convert
            bsa_nm2 = r.get('mean_bsa_nm2', r.get('bsa_nm2'))
            if bsa_nm2 is not None:
                bsa = bsa_nm2 * 100.0
        if bsa is None:
            print(f'WARNING: No BSA found in {exp_id} results — skipping')
            missing.append(exp_id)
            continue

        # Extract computed ΔG (use computed value, not experimental)
        dg = r.get('dg_bind_kcal_mol', r.get('delta_g_kcal_mol', r.get('dg_kcal_mol')))
        if dg is None:
            dg_kj = r.get('dg_bind_kj_mol', r.get('delta_g_kj_mol'))
            if dg_kj is not None:
                dg = dg_kj / 4.184
        if dg is None:
            print(f'WARNING: No ΔG found in {exp_id} results — skipping')
            missing.append(exp_id)
            continue

        bsa_values.append(float(bsa))
        dg_values.append(float(dg))
        labels.append(f"{info['name']}\n({exp_id})")

    if len(bsa_values) < 3:
        print(f'ERROR: Need at least 3 data points, got {len(bsa_values)}')
        sys.exit(1)

    bsa_arr = np.array(bsa_values)
    dg_arr = np.array(dg_values)

    # Linear regression: ΔG = slope × BSA + intercept
    coeffs = np.polyfit(bsa_arr, dg_arr, 1)
    slope, intercept = coeffs
    dg_fit = np.polyval(coeffs, bsa_arr)

    # R²
    ss_res = np.sum((dg_arr - dg_fit) ** 2)
    ss_tot = np.sum((dg_arr - np.mean(dg_arr)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # 95% CI for slope (approximate)
    n = len(bsa_arr)
    if n > 2:
        se_slope = np.sqrt(ss_res / (n - 2) / np.sum((bsa_arr - np.mean(bsa_arr)) ** 2))
        from scipy import stats as sp_stats
        t_val = sp_stats.t.ppf(0.975, n - 2)
        slope_ci = (slope - t_val * se_slope, slope + t_val * se_slope)
    else:
        slope_ci = (slope, slope)

    # Classification
    r2_pass = r_squared > R2_PASS
    slope_pass = SLOPE_RANGE[0] <= slope <= SLOPE_RANGE[1]
    verdict = 'PASS' if (r2_pass and slope_pass) else (
        'MARGINAL' if r2_pass else 'FAIL')

    print(f'\nLinear fit: ΔG = {slope:.5f} × BSA + {intercept:.2f}')
    print(f'R² = {r_squared:.3f} (target: >{R2_PASS})')
    print(f'Slope = {slope:.5f} kcal/(mol·Å²) (target: [{SLOPE_RANGE[0]}, {SLOPE_RANGE[1]}])')
    print(f'Slope 95% CI: [{slope_ci[0]:.5f}, {slope_ci[1]:.5f}]')
    print(f'Verdict: {verdict}')

    results = {
        'experiment_id': 'EXP-26',
        'n_systems': n,
        'missing_systems': missing,
        'slope_kcal_per_mol_ang2': float(slope),
        'intercept_kcal_mol': float(intercept),
        'r_squared': float(r_squared),
        'slope_95ci': [float(slope_ci[0]), float(slope_ci[1])],
        'data_points': [{'label': labels[i].replace('\n', ' '),
                         'bsa_ang2': float(bsa_arr[i]),
                         'dg_kcal_mol': float(dg_arr[i])}
                        for i in range(n)],
        'verdict': verdict,
    }

    with open(output_dir / 'bsa_dg_correlation.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Figure: BSA–ΔG scatter with regression line
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(bsa_arr, dg_arr, c='#1f77b4', s=80, zorder=3)
    for i, label in enumerate(labels):
        ax.annotate(label, (bsa_arr[i], dg_arr[i]),
                    textcoords='offset points', xytext=(8, 5), fontsize=7)

    bsa_range = np.linspace(bsa_arr.min() - 100, bsa_arr.max() + 100, 100)
    ax.plot(bsa_range, np.polyval(coeffs, bsa_range), 'r--', lw=1.5,
            label=f'y = {slope:.5f}x + {intercept:.1f}\nR² = {r_squared:.3f}')
    ax.set_xlabel('Buried Surface Area (Å²)')
    ax.set_ylabel('ΔG_bind (kcal/mol)')
    ax.set_title('EXP-26: BSA–ΔG Correlation')
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / 'bsa_dg_scatter.png', dpi=150)
    plt.close(fig)
    print(f'Figure saved: {figures_dir / "bsa_dg_scatter.png"}')


if __name__ == '__main__':
    main()
