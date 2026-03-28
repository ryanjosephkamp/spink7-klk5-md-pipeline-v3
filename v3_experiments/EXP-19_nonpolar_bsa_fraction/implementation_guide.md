# EXP-19: Nonpolar BSA Fraction — Implementation Guide

**Experiment ID:** EXP-19  
**Feature ID:** F-19 (benchmarks.md)  
**Category:** Structural (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies the fraction of BSA that is nonpolar (carbon atoms). BPTI–trypsin benchmark: 61% nonpolar (95% CI [55%, 67%]). Shared data with EXP-16. Reflects hydrophobic burial driving binding. Cross-ref EXP-26 (BSA–ΔG correlation).

### 2. Classification (§25.1)

- PASS: Nonpolar fraction in [55%, 67%]
- MARGINAL: In [48%, 74%]
- FAIL: Outside marginal range

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-19_nonpolar_bsa_fraction"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load Trajectory

```python
EXP_ROOT = PROJECT_ROOT / "v3_experiments"
traj_path = EXP_ROOT / "EXP-04_bpti_trypsin_dg_bind/outputs/production/production_trajectory.dcd"
top_path = EXP_ROOT / "EXP-04_bpti_trypsin_dg_bind/outputs/solvated_complex.pdb"

traj = load_trajectory(str(traj_path), str(top_path), stride=10)
topo = traj.topology

n_skip = traj.n_frames // 2
traj_analysis = traj[n_skip:]

protease_atoms = topo.select("chainid 0 and not water")
inhibitor_atoms = topo.select("chainid 1 and not water")
complex_atoms = np.concatenate([protease_atoms, inhibitor_atoms])
print(f"Analyzing {traj_analysis.n_frames} frames")
```

### Step 3: Atom Classification

```python
# Classify atoms as polar (N, O, S) or nonpolar (C)
def classify_atoms(topology_obj, atom_indices):
    """Return masks for polar and nonpolar atoms."""
    polar_mask = []
    nonpolar_mask = []
    for idx in atom_indices:
        elem = topology_obj.atom(idx).element.symbol
        if elem in ('N', 'O', 'S'):
            polar_mask.append(True)
            nonpolar_mask.append(False)
        elif elem == 'C':
            polar_mask.append(False)
            nonpolar_mask.append(True)
        else:
            polar_mask.append(False)
            nonpolar_mask.append(False)
    return np.array(polar_mask), np.array(nonpolar_mask)

polar_prot, nonpolar_prot = classify_atoms(topo, protease_atoms)
polar_inh, nonpolar_inh = classify_atoms(topo, inhibitor_atoms)
print(f"Protease: {np.sum(polar_prot)} polar, {np.sum(nonpolar_prot)} nonpolar atoms")
print(f"Inhibitor: {np.sum(polar_inh)} polar, {np.sum(nonpolar_inh)} nonpolar atoms")
```

### Step 4: Per-Frame Nonpolar/Polar BSA Decomposition

```python
nonpolar_fractions = []
bsa_total_ts = []
bsa_nonpolar_ts = []
bsa_polar_ts = []

for i in range(traj_analysis.n_frames):
    frame = traj_analysis[i]

    # SASA per atom for complex, protease, inhibitor
    sasa_complex = md.shrake_rupley(frame.atom_slice(complex_atoms), probe_radius=0.14, mode='atom')[0] * 100
    sasa_prot = md.shrake_rupley(frame.atom_slice(protease_atoms), probe_radius=0.14, mode='atom')[0] * 100
    sasa_inh = md.shrake_rupley(frame.atom_slice(inhibitor_atoms), probe_radius=0.14, mode='atom')[0] * 100

    # Per-atom BSA = SASA(free) - SASA(bound)
    # For protease atoms: BSA_prot = sasa_prot - sasa_complex[:len(protease_atoms)]
    bsa_prot_atoms = sasa_prot - sasa_complex[:len(protease_atoms)]
    bsa_inh_atoms = sasa_inh - sasa_complex[len(protease_atoms):]

    # Total BSA
    bsa_total = np.sum(np.maximum(bsa_prot_atoms, 0)) + np.sum(np.maximum(bsa_inh_atoms, 0))

    # Nonpolar BSA (carbon atoms only)
    bsa_nonpolar = (np.sum(np.maximum(bsa_prot_atoms[nonpolar_prot], 0)) +
                    np.sum(np.maximum(bsa_inh_atoms[nonpolar_inh], 0)))

    # Polar BSA
    bsa_polar = (np.sum(np.maximum(bsa_prot_atoms[polar_prot], 0)) +
                 np.sum(np.maximum(bsa_inh_atoms[polar_inh], 0)))

    bsa_total_ts.append(bsa_total)
    bsa_nonpolar_ts.append(bsa_nonpolar)
    bsa_polar_ts.append(bsa_polar)
    if bsa_total > 0:
        nonpolar_fractions.append(bsa_nonpolar / bsa_total * 100)
    else:
        nonpolar_fractions.append(0)

    if (i + 1) % 100 == 0:
        print(f"Frame {i+1}: BSA={bsa_total:.0f} Å², nonpolar={nonpolar_fractions[-1]:.1f}%")

nonpolar_fractions = np.array(nonpolar_fractions)
np_mean = np.mean(nonpolar_fractions)
np_std = np.std(nonpolar_fractions)
np_cv = np_std / np_mean * 100
print(f"Nonpolar fraction: {np_mean:.1f} ± {np_std:.1f}% (CV = {np_cv:.1f}%)")
```

### Step 5: Classification

```python
if 55 <= np_mean <= 67:
    classification = "PASS"
elif 48 <= np_mean <= 74:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-19", "feature_id": "F-19",
    "nonpolar_fraction_pct": float(np_mean),
    "nonpolar_fraction_std": float(np_std),
    "nonpolar_fraction_cv": float(np_cv),
    "mean_bsa_total": float(np.mean(bsa_total_ts)),
    "mean_bsa_nonpolar": float(np.mean(bsa_nonpolar_ts)),
    "mean_bsa_polar": float(np.mean(bsa_polar_ts)),
    "benchmark": "61 ± 6%",
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-19: {np_mean:.1f}% nonpolar → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Nonpolar Fraction Time Series

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 6))
time_ns = np.arange(len(nonpolar_fractions)) * 0.1 + 50
ax.plot(time_ns, nonpolar_fractions, 'b-', linewidth=0.5, alpha=0.7)
ax.axhline(y=61, color='green', linestyle='--', linewidth=2, label='Benchmark (61%)')
ax.axhspan(55, 67, alpha=0.1, color='green', label='PASS range')
ax.axhline(y=np_mean, color='red', linewidth=1.5, label=f'Mean = {np_mean:.1f}%')
ax.set_xlabel('Time (ns)', fontsize=14)
ax.set_ylabel('Nonpolar BSA Fraction (%)', fontsize=14)
ax.set_title(f'EXP-19: Nonpolar BSA Fraction — {classification}', fontsize=15)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-19_nonpolar_fraction_timeseries.png", dpi=300)
plt.close(fig)
```

### Figure 2: BSA Composition Pie Chart

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
sizes = [np.mean(bsa_nonpolar_ts), np.mean(bsa_polar_ts)]
labels = [f'Nonpolar\n{sizes[0]:.0f} Å² ({np_mean:.0f}%)',
          f'Polar\n{sizes[1]:.0f} Å² ({100-np_mean:.0f}%)']
colors = ['#FFD700', '#4169E1']
ax.pie(sizes, labels=labels, colors=colors, autopct='', startangle=90,
       textprops={'fontsize': 13}, wedgeprops={'edgecolor': 'black'})
ax.set_title(f'EXP-19: BSA Composition — {classification}', fontsize=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-19_bsa_pie.png", dpi=300)
plt.close(fig)
```

### Figure 3: Distribution Histogram

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(nonpolar_fractions, bins=40, color='gold', edgecolor='black', density=True)
ax.axvline(x=61, color='green', linestyle='--', linewidth=2, label='Benchmark (61%)')
ax.axvspan(55, 67, alpha=0.1, color='green', label='PASS range')
ax.set_xlabel('Nonpolar Fraction (%)', fontsize=14)
ax.set_ylabel('Density', fontsize=14)
ax.set_title('EXP-19: Nonpolar Fraction Distribution', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-19_distribution.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-19: Nonpolar BSA Fraction — Results Report

**Experiment ID:** EXP-19  **Feature ID:** F-19  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value |
|--------|-------|
| Nonpolar fraction | [val] ± [std]% |
| Nonpolar BSA | [val] Å² |
| Polar BSA | [val] Å² |
| Total BSA | [val] Å² |
| Benchmark | 61 ± 6% |

## Figures
1. Nonpolar fraction time series
2. BSA composition pie chart
3. Distribution histogram

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
