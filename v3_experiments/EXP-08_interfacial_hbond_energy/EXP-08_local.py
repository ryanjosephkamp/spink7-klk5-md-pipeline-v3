#!/usr/bin/env python3
"""EXP-08: Interfacial H-Bond Energy — Local Analysis Script

Runs locally on CPU after downloading GPU results from Google Drive.
Requires: EXP-04 production trajectory and topology.

Usage:
    python EXP-08_local.py --data-dir /path/to/v3_gpu_results/EXP-04
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

from src.analyze.contacts import compute_hbonds
from src.config import BOLTZMANN_KJ, KCAL_TO_KJ

TEMPERATURE_K = 310.0


def main():
    parser = argparse.ArgumentParser(description='EXP-08: Interfacial H-bond energy analysis')
    parser.add_argument('--data-dir', type=Path, required=True,
                        help='Path to EXP-04 GPU results directory')
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

    # Identify chains
    chains = list(traj.topology.chains)
    chain_a_indices = np.array([a.index for a in chains[0].atoms])  # trypsin
    chain_b_indices = np.array([a.index for a in chains[1].atoms])  # BPTI

    # Compute interfacial H-bonds
    print('Computing interfacial H-bonds...')
    hbond_result = compute_hbonds(traj, chain_a_indices, chain_b_indices)

    triplets = hbond_result['hbond_triplets']  # [N_hbonds, 3]
    frequencies = hbond_result['hbond_frequency']  # [N_hbonds]

    if len(frequencies) == 0:
        print('WARNING: No interfacial H-bonds detected')
        results = {'experiment_id': 'EXP-08', 'n_hbonds': 0, 'verdict': 'FAIL'}
        with open(output_dir / 'hbond_analysis.json', 'w') as f:
            json.dump(results, f, indent=2)
        return

    # Compute per-H-bond energy from occupancy: E = -k_B T ln(occupancy)
    kbt_kj = BOLTZMANN_KJ * TEMPERATURE_K  # kJ/mol
    kbt_kcal = kbt_kj / KCAL_TO_KJ  # kcal/mol

    # Filter out zero-occupancy bonds
    nonzero_mask = frequencies > 0.0
    valid_freq = frequencies[nonzero_mask]
    valid_triplets = triplets[nonzero_mask]

    energies_kcal = -kbt_kcal * np.log(valid_freq)

    # H-bond labels
    hbond_labels = []
    for i, (d, h, a) in enumerate(valid_triplets):
        d_atom = traj.topology.atom(d)
        a_atom = traj.topology.atom(a)
        label = f'{d_atom.residue.name}{d_atom.residue.resSeq}:{d_atom.name}...{a_atom.residue.name}{a_atom.resSeq}:{a_atom.name}'
        hbond_labels.append(label)

    # Classification: per-H-bond energy in [0.7, 2.3] kcal/mol
    median_energy = float(np.median(energies_kcal))
    mean_energy = float(np.mean(energies_kcal))
    in_range = np.sum((energies_kcal >= 0.7) & (energies_kcal <= 2.3))
    fraction_in_range = float(in_range) / len(energies_kcal)

    verdict = 'PASS' if (0.7 <= median_energy <= 2.3 and fraction_in_range >= 0.5) else (
        'MARGINAL' if (0.5 <= median_energy <= 3.0) else 'FAIL')

    results = {
        'experiment_id': 'EXP-08',
        'n_hbonds_total': int(len(frequencies)),
        'n_hbonds_nonzero': int(len(valid_freq)),
        'median_energy_kcal_mol': median_energy,
        'mean_energy_kcal_mol': mean_energy,
        'fraction_in_target_range': fraction_in_range,
        'top_hbonds': [{'label': hbond_labels[i], 'occupancy': float(valid_freq[i]),
                        'energy_kcal_mol': float(energies_kcal[i])}
                       for i in np.argsort(valid_freq)[::-1][:10]],
        'verdict': verdict,
    }

    with open(output_dir / 'hbond_analysis.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\n{len(valid_freq)} interfacial H-bonds with nonzero occupancy')
    print(f'Median per-H-bond energy: {median_energy:.2f} kcal/mol (target: [0.7, 2.3])')
    print(f'Fraction in target range: {fraction_in_range:.1%}')
    print(f'Verdict: {verdict}')

    # Figure: H-bond occupancy
    sorted_idx = np.argsort(valid_freq)[::-1]
    top_n = min(20, len(valid_freq))
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(range(top_n), valid_freq[sorted_idx[:top_n]])
    ax.set_xlabel('H-bond index (ranked)')
    ax.set_ylabel('Occupancy')
    ax.set_title('EXP-08: Interfacial H-Bond Occupancy (Top 20)')
    fig.tight_layout()
    fig.savefig(figures_dir / 'hbond_occupancy.png', dpi=150)
    plt.close(fig)

    # Figure: energy distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(energies_kcal, bins=20, edgecolor='black', alpha=0.7)
    ax.axvspan(0.7, 2.3, alpha=0.15, color='green', label='Target range [0.7, 2.3]')
    ax.axvline(median_energy, color='red', linestyle='--', label=f'Median = {median_energy:.2f}')
    ax.set_xlabel('Per-H-bond Energy (kcal/mol)')
    ax.set_ylabel('Count')
    ax.set_title('EXP-08: H-Bond Energy Distribution')
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / 'hbond_energy_distribution.png', dpi=150)
    plt.close(fig)
    print(f'Figures saved to {figures_dir}')


if __name__ == '__main__':
    main()
