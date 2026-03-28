#!/usr/bin/env python3
"""EXP-07: P1 Energetic Contribution — Local Analysis Script

Runs locally on CPU after downloading GPU results from Google Drive.
Requires: EXP-04 production trajectory and topology.

Usage:
    python EXP-07_local.py --data-dir /path/to/v3_gpu_results/EXP-04
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

from src.analyze.contacts import compute_interface_contacts
from src.analyze.structural import compute_sasa


def main():
    parser = argparse.ArgumentParser(description='EXP-07: P1 energetic contribution analysis')
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
    if not topo_path.exists():
        print(f'ERROR: No topology found at {topo_path}')
        sys.exit(1)

    print(f'Loading trajectory: {traj_path}')
    traj = md.load(str(traj_path), top=str(topo_path))
    print(f'Loaded {traj.n_frames} frames, {traj.n_atoms} atoms')

    # Identify chains
    chains = list(traj.topology.chains)
    if len(chains) < 2:
        print('ERROR: Expected at least 2 chains (trypsin + BPTI)')
        sys.exit(1)

    chain_a_indices = np.array([a.index for a in chains[0].atoms])  # trypsin
    chain_b_indices = np.array([a.index for a in chains[1].atoms])  # BPTI

    # Compute interface contacts
    print('Computing interface contacts...')
    contact_result = compute_interface_contacts(
        traj, chain_a_indices, chain_b_indices, cutoff_nm=0.45)

    contact_freq = contact_result['contact_frequency']  # [N_res_a, N_res_b]

    # Sum contact frequency per BPTI residue (axis 0 = trypsin, axis 1 = BPTI)
    per_bpti_residue = np.sum(contact_freq, axis=0)

    # Map residue indices to residue names/numbers
    bpti_residues = []
    bpti_chain = chains[1]
    for res in bpti_chain.residues:
        bpti_residues.append(f'{res.name}{res.resSeq}')

    # Rank by contact fraction
    ranked_indices = np.argsort(per_bpti_residue)[::-1]
    ranked_residues = [(bpti_residues[i], float(per_bpti_residue[i])) for i in ranked_indices]

    # Find K15 (P1 residue)
    k15_rank = None
    k15_fraction = None
    for rank, (res_label, frac) in enumerate(ranked_residues):
        if '15' in res_label and ('LYS' in res_label or 'K' in res_label.upper()):
            k15_rank = rank + 1
            k15_fraction = frac
            break

    # Compute BSA
    print('Computing BSA...')
    sasa_complex = compute_sasa(traj)
    sasa_a = compute_sasa(traj.atom_slice(chain_a_indices))
    sasa_b = compute_sasa(traj.atom_slice(chain_b_indices))
    bsa = sasa_a + sasa_b - sasa_complex  # nm²
    mean_bsa_nm2 = float(np.mean(bsa))
    mean_bsa_ang2 = mean_bsa_nm2 * 100.0  # convert nm² to Å²

    # Classification
    k15_is_top = k15_rank == 1 if k15_rank is not None else False
    k15_above_40 = k15_fraction > 0.40 if k15_fraction is not None else False
    verdict = 'PASS' if (k15_is_top and k15_above_40) else (
        'MARGINAL' if (k15_rank is not None and k15_rank <= 3) else 'FAIL')

    results = {
        'experiment_id': 'EXP-07',
        'ranked_residues_top10': ranked_residues[:10],
        'k15_rank': k15_rank,
        'k15_contact_fraction': k15_fraction,
        'mean_bsa_nm2': mean_bsa_nm2,
        'mean_bsa_ang2': mean_bsa_ang2,
        'verdict': verdict,
    }

    with open(output_dir / 'residue_contacts.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f'\nResults: K15 rank={k15_rank}, fraction={k15_fraction:.3f}')
    print(f'Mean BSA = {mean_bsa_ang2:.1f} Å²')
    print(f'Verdict: {verdict}')

    # Figure: per-residue contact fractions (top 15)
    top_n = min(15, len(ranked_residues))
    labels = [r[0] for r in ranked_residues[:top_n]]
    fracs = [r[1] for r in ranked_residues[:top_n]]

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ['#d62728' if '15' in l else '#1f77b4' for l in labels]
    ax.bar(range(top_n), fracs, color=colors)
    ax.set_xticks(range(top_n))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Contact Fraction')
    ax.set_title('EXP-07: Per-Residue Interface Contact Fractions (BPTI)')
    ax.axhline(y=0.40, color='gray', linestyle='--', alpha=0.5, label='40% threshold')
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / 'per_residue_contacts.png', dpi=150)
    plt.close(fig)
    print(f'Figure saved: {figures_dir / "per_residue_contacts.png"}')


if __name__ == '__main__':
    main()
