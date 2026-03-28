# EXP-17: Interfacial Hydrogen Bond Count — Implementation Guide

**Experiment ID:** EXP-17  
**Feature ID:** F-17 (benchmarks.md)  
**Category:** Structural (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies persistent intermolecular H-bonds at the protease–inhibitor interface. BPTI–trypsin benchmark: ~10 persistent H-bonds (95% CI [7, 13]). Uses Baker-Hubbard criteria: d(D···A) ≤ 3.5 Å, angle ≥ 135°. H-bond occupancy table, cross-ref EXP-08 (per-H-bond energy) and EXP-18 (water-mediated bridges).

### 2. Classification (§25.1)

- PASS: Persistent H-bond count in [7, 13]
- MARGINAL: Count in [5, 16]
- FAIL: < 5 or > 16

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-17_interfacial_hbond_count"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory
from src.analyze.structural import compute_hbonds
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load Trajectory

```python
EXP_ROOT = PROJECT_ROOT / "v3_experiments"
traj_path = EXP_ROOT / "EXP-04_bpti_trypsin_dg_bind/outputs/production/production_trajectory.dcd"
top_path = EXP_ROOT / "EXP-04_bpti_trypsin_dg_bind/outputs/solvated_complex.pdb"

traj = load_trajectory(str(traj_path), str(top_path), stride=1)
topo = traj.topology

protease_residues = set(r.index for r in topo.residues if r.chain.index == 0)
inhibitor_residues = set(r.index for r in topo.residues if r.chain.index == 1)
print(f"Protease: {len(protease_residues)} residues, Inhibitor: {len(inhibitor_residues)} residues")
```

### Step 3: Crystal Structure H-Bond Count

```python
from src.prep.pdb_fetch import fetch_pdb
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)
pdb_2ptc = fetch_pdb("2PTC", data_dir)
crystal = md.load(str(pdb_2ptc))

# Baker-Hubbard H-bonds in crystal
hbonds_crystal = md.baker_hubbard(crystal, freq=0.0)
interchain_crystal = []
for hb in hbonds_crystal:
    d_res = crystal.topology.atom(hb[0]).residue.chain.index
    a_res = crystal.topology.atom(hb[2]).residue.chain.index
    if d_res != a_res:
        interchain_crystal.append(hb)
print(f"Crystal interfacial H-bonds: {len(interchain_crystal)}")
```

### Step 4: Per-Frame H-Bond Analysis

```python
# Compute H-bonds for each frame
hbond_counts = []
hbond_pairs = {}  # Track unique D-A pairs and their occupancy

for i in range(traj.n_frames):
    frame = traj[i]
    hbonds = md.baker_hubbard(frame, freq=0.0)

    interchain = []
    for hb in hbonds:
        d_chain = topo.atom(hb[0]).residue.chain.index
        a_chain = topo.atom(hb[2]).residue.chain.index
        if d_chain != a_chain:
            interchain.append(hb)
            # Track pair
            d_name = f"{topo.atom(hb[0]).residue.name}{topo.atom(hb[0]).residue.resSeq}-{topo.atom(hb[0]).name}"
            a_name = f"{topo.atom(hb[2]).residue.name}{topo.atom(hb[2]).residue.resSeq}-{topo.atom(hb[2]).name}"
            pair_key = f"{d_name} → {a_name}"
            hbond_pairs[pair_key] = hbond_pairs.get(pair_key, 0) + 1

    hbond_counts.append(len(interchain))
    if (i + 1) % 1000 == 0:
        print(f"Frame {i+1}/{traj.n_frames}: {len(interchain)} interfacial H-bonds")

hbond_counts = np.array(hbond_counts)
```

### Step 5: Persistent H-Bond Identification

```python
# Persistent = present in >50% of frames
n_frames = traj.n_frames
persistent_threshold = 0.5

occupancy_table = []
for pair, count in sorted(hbond_pairs.items(), key=lambda x: -x[1]):
    occ = count / n_frames * 100
    occupancy_table.append({"pair": pair, "occupancy_pct": occ, "n_frames": count})

persistent_hbonds = [h for h in occupancy_table if h["occupancy_pct"] > 50]
n_persistent = len(persistent_hbonds)

print(f"\nTotal H-bond count: mean={np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f}")
print(f"Persistent H-bonds (>50% occupancy): {n_persistent}")
print("\nTop 15 H-bonds by occupancy:")
for h in occupancy_table[:15]:
    print(f"  {h['pair']}: {h['occupancy_pct']:.1f}%")
```

### Step 6: Classification

```python
if 7 <= n_persistent <= 13:
    classification = "PASS"
elif 5 <= n_persistent <= 16:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-17", "feature_id": "F-17",
    "n_persistent_hbonds": n_persistent,
    "mean_total_hbonds": float(np.mean(hbond_counts)),
    "std_total_hbonds": float(np.std(hbond_counts)),
    "crystal_hbonds": len(interchain_crystal),
    "benchmark": "10 ± 3",
    "persistent_hbond_table": persistent_hbonds[:20],
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-17: {n_persistent} persistent H-bonds → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: H-Bond Count Time Series

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 6))
time_ns = np.arange(len(hbond_counts)) * 0.01
ax.plot(time_ns, hbond_counts, 'b-', linewidth=0.3, alpha=0.5)
window = 100
running_avg = np.convolve(hbond_counts, np.ones(window)/window, mode='valid')
ax.plot(time_ns[:len(running_avg)], running_avg, 'r-', linewidth=1.5, label='Running avg (1 ns)')
ax.axhspan(7, 13, alpha=0.1, color='green', label='PASS range')
ax.axhline(y=10, color='green', linestyle='--', linewidth=1.5, label='Benchmark (10)')
ax.set_xlabel('Time (ns)', fontsize=14)
ax.set_ylabel('Interfacial H-bond Count', fontsize=14)
ax.set_title(f'EXP-17: Interfacial H-Bonds — {classification}', fontsize=15)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-17_hbond_timeseries.png", dpi=300)
plt.close(fig)
```

### Figure 2: H-Bond Occupancy Bar Chart

```python
fig, ax = plt.subplots(1, 1, figsize=(14, 8))
top_n = min(20, len(occupancy_table))
pairs = [h["pair"] for h in occupancy_table[:top_n]]
occs = [h["occupancy_pct"] for h in occupancy_table[:top_n]]
colors = ['green' if o > 50 else 'orange' if o > 25 else 'lightgray' for o in occs]
ax.barh(pairs[::-1], occs[::-1], color=colors[::-1], edgecolor='black')
ax.axvline(x=50, color='red', linestyle='--', label='Persistent threshold (50%)')
ax.set_xlabel('Occupancy (%)', fontsize=12)
ax.set_title(f'EXP-17: H-Bond Occupancy (Top {top_n})', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-17_hbond_occupancy.png", dpi=300)
plt.close(fig)
```

### Figure 3: H-Bond Count Distribution

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(hbond_counts, bins=range(0, max(hbond_counts)+2), color='steelblue',
        edgecolor='black', density=True, align='left')
ax.axvline(x=10, color='green', linestyle='--', linewidth=2, label='Benchmark (10)')
ax.axvspan(7, 13, alpha=0.1, color='green', label='PASS range')
ax.set_xlabel('Interfacial H-bond Count', fontsize=14)
ax.set_ylabel('Probability', fontsize=14)
ax.set_title('EXP-17: H-Bond Count Distribution', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-17_hbond_distribution.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-17: Interfacial H-Bond Count — Results Report

**Experiment ID:** EXP-17  **Feature ID:** F-17  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value |
|--------|-------|
| Persistent H-bonds (>50%) | [n] |
| Mean total H-bonds/frame | [val] ± [std] |
| Crystal H-bonds | [n] |
| Benchmark | 10 ± 3 |

## Key H-bonds
| Donor → Acceptor | Occupancy (%) |
|-------------------|---------------|
| [pair 1] | [val] |
| ... | ... |

## Figures
1. H-bond count time series
2. Occupancy bar chart
3. Count distribution

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
