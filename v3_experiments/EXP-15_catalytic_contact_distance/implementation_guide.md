# EXP-15: Catalytic Contact Distance — Implementation Guide

**Experiment ID:** EXP-15  
**Feature ID:** F-15 (benchmarks.md)  
**Category:** Structural (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates that MD simulations maintain the catalytic contact distance between Ser195 Oγ and the P1 carbonyl carbon at ~2.7 Å (95% CI [2.4, 3.0] Å). Tested on BPTI–trypsin (PDB 2PTC) primary, with PSTI (1TGS) and SPINK7–KLK5 (EXP-01) secondary. Also monitors complete catalytic triad geometry (Ser195–His57–Asp102).

### 2. Benchmark Values

| Metric | Value | Source |
|--------|-------|--------|
| d(Ser195 Oγ, P1 C) | 2.7 ± 0.3 Å | Radisky & Bhatt, Vincent & Bhatt 2007 |
| d(S195 Oγ, H57 Nε2) | ~2.8 Å | Catalytic triad |
| d(H57 Nδ1, D102 Oδ) | ~2.8 Å | Catalytic triad |
| Frames >3.5 Å | <10% | Structural stability |

### 3. Classification (§25.1)

- PASS: Mean in [2.4, 3.0] Å
- MARGINAL: Mean in [2.0, 3.5] Å
- FAIL: Outside [2.0, 3.5] Å or >20% frames >4.0 Å

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-15_catalytic_contact_distance"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory
from src.analyze.structural import compute_rmsd
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load Trajectories from Upstream Experiments

```python
EXP_ROOT = PROJECT_ROOT / "v3_experiments"

SYSTEMS = {
    "BPTI-trypsin": {
        "exp": "EXP-04_bpti_trypsin_dg_bind",
        "protease_chain": 0,
        "inhibitor_chain": 1,
        "p1_resname": "LYS",  # K15
        "ser_resname": "SER",  # Ser195
        "his_resname": "HIS",  # His57
        "asp_resname": "ASP",  # Asp102
    },
    "SPINK7-KLK5": {
        "exp": "EXP-01_spink7_klk5_binding",
        "protease_chain": 0,
        "inhibitor_chain": 1,
        "p1_resname": "ARG",
        "ser_resname": "SER",
        "his_resname": "HIS",
        "asp_resname": "ASP",
    },
}

trajectories = {}
for name, info in SYSTEMS.items():
    traj_path = EXP_ROOT / info["exp"] / "outputs/production/production_trajectory.dcd"
    top_path = EXP_ROOT / info["exp"] / "outputs/solvated_complex.pdb"
    if traj_path.exists() and top_path.exists():
        traj = load_trajectory(str(traj_path), str(top_path), stride=1)
        trajectories[name] = traj
        print(f"{name}: {traj.n_frames} frames loaded")
    else:
        print(f"WARNING: {name} trajectory not found")
```

### Step 3: Crystal Structure Reference

```python
from src.prep.pdb_fetch import fetch_pdb
from src.prep.pdb_clean import clean_structure

data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

pdb_2ptc = fetch_pdb("2PTC", data_dir)
crystal = md.load(str(pdb_2ptc))
topo_c = crystal.topology

# Measure crystal structure distance
ser_og_c = topo_c.select("resname SER and name OG and chainid 0")
lys_c_c = topo_c.select("resname LYS and name C and chainid 1")

if len(ser_og_c) > 0 and len(lys_c_c) > 0:
    d_crystal = md.compute_distances(crystal, [[ser_og_c[0], lys_c_c[0]]])[0, 0] * 10
    print(f"Crystal structure d(Ser Oγ, P1 C) = {d_crystal:.2f} Å")
```

### Step 4: Catalytic Distance Analysis on MD Trajectories

```python
results_all = {}

for name, traj in trajectories.items():
    info = SYSTEMS[name]
    topo = traj.topology

    # Ser195 Oγ (on protease chain)
    ser_og = topo.select(f"resname SER and name OG and chainid {info['protease_chain']}")
    # His57 Nε2
    his_ne2 = topo.select(f"resname HIS and name NE2 and chainid {info['protease_chain']}")
    # Asp102 Oδ
    asp_od = topo.select(f"resname ASP and name OD1 and chainid {info['protease_chain']}")
    # P1 carbonyl C (on inhibitor chain)
    p1_c = topo.select(f"resname {info['p1_resname']} and name C and chainid {info['inhibitor_chain']}")

    if len(ser_og) == 0 or len(p1_c) == 0:
        print(f"WARNING: {name} — could not locate atoms")
        continue

    # d1: Ser Oγ → P1 C (catalytic contact)
    d1 = md.compute_distances(traj, [[ser_og[0], p1_c[0]]])[:, 0] * 10  # Å

    # d2: Ser Oγ → His Nε2
    d2 = None
    if len(his_ne2) > 0:
        d2 = md.compute_distances(traj, [[ser_og[0], his_ne2[0]]])[:, 0] * 10

    # d3: His Nδ1 → Asp Oδ
    his_nd1 = topo.select(f"resname HIS and name ND1 and chainid {info['protease_chain']}")
    d3 = None
    if len(his_nd1) > 0 and len(asp_od) > 0:
        d3 = md.compute_distances(traj, [[his_nd1[0], asp_od[0]]])[:, 0] * 10

    # Block averaging (10 blocks)
    n_blocks = 10
    block_size = len(d1) // n_blocks
    block_means = [np.mean(d1[i*block_size:(i+1)*block_size]) for i in range(n_blocks)]
    block_se = np.std(block_means) / np.sqrt(n_blocks)

    frac_above_35 = np.mean(d1 > 3.5) * 100
    frac_above_40 = np.mean(d1 > 4.0) * 100

    results_all[name] = {
        "d1": d1, "d2": d2, "d3": d3,
        "mean_d1": float(np.mean(d1)),
        "std_d1": float(np.std(d1)),
        "se_d1": float(block_se),
        "frac_above_35": float(frac_above_35),
        "frac_above_40": float(frac_above_40),
        "block_means": block_means,
    }
    print(f"{name}: d(Ser Oγ, P1 C) = {np.mean(d1):.2f} ± {block_se:.2f} Å, "
          f">3.5 Å: {frac_above_35:.1f}%, >4.0 Å: {frac_above_40:.1f}%")
```

### Step 5: Classification

```python
# Primary system classification (BPTI-trypsin)
primary = results_all.get("BPTI-trypsin", {})
mean_d = primary.get("mean_d1", 0)
frac_40 = primary.get("frac_above_40", 100)

if 2.4 <= mean_d <= 3.0:
    classification = "PASS"
elif 2.0 <= mean_d <= 3.5 and frac_40 < 20:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-15", "feature_id": "F-15",
    "crystal_reference_ang": float(d_crystal) if 'd_crystal' in dir() else None,
    "systems": {name: {k: v for k, v in data.items() if k != 'd1' and k != 'd2' and k != 'd3'}
                for name, data in results_all.items()},
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-15 Classification: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Catalytic Distance Time Series

```python
fig, axes = plt.subplots(len(results_all), 1, figsize=(14, 4*len(results_all)), sharex=True)
if len(results_all) == 1: axes = [axes]

for ax, (name, data) in zip(axes, results_all.items()):
    d1 = data["d1"]
    time_ns = np.arange(len(d1)) * 0.01
    ax.plot(time_ns, d1, 'b-', linewidth=0.3, alpha=0.5)
    # Running average (1 ns window)
    window = min(100, len(d1) // 5)
    running_avg = np.convolve(d1, np.ones(window)/window, mode='valid')
    ax.plot(time_ns[:len(running_avg)], running_avg, 'r-', linewidth=1.5, label='1 ns running avg')
    ax.axhline(y=2.7, color='green', linestyle='--', linewidth=1.5, label='Crystal (2.7 Å)')
    ax.axhspan(2.4, 3.0, alpha=0.1, color='green', label='PASS [2.4, 3.0]')
    ax.set_ylabel('d(Ser Oγ, P1 C) (Å)', fontsize=11)
    ax.set_title(f'{name}: mean = {data["mean_d1"]:.2f} ± {data["se_d1"]:.2f} Å', fontsize=12)
    ax.legend(fontsize=9)

axes[-1].set_xlabel('Time (ns)', fontsize=12)
plt.suptitle(f'EXP-15: Catalytic Contact Distance — {classification}', fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-15_catalytic_distance_timeseries.png", dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Figure 2: Distance Histogram

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
for name, data in results_all.items():
    ax.hist(data["d1"], bins=80, alpha=0.5, density=True, label=name)
ax.axvline(x=2.7, color='green', linestyle='--', linewidth=2, label='Crystal (2.7 Å)')
ax.axvspan(2.4, 3.0, alpha=0.1, color='green')
ax.set_xlabel('d(Ser Oγ, P1 C) (Å)', fontsize=14)
ax.set_ylabel('Density', fontsize=14)
ax.set_title(f'EXP-15: Catalytic Contact Distance Distribution — {classification}', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-15_distance_histogram.png", dpi=300)
plt.close(fig)
```

### Figure 3: Complete Catalytic Triad Distances

```python
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
primary_data = results_all.get("BPTI-trypsin", {})
time_ns = np.arange(len(primary_data.get("d1", []))) * 0.01

labels_thresholds = [
    ("d(Ser Oγ → P1 C)", primary_data.get("d1"), 3.5, 'blue'),
    ("d(Ser Oγ → His Nε2)", primary_data.get("d2"), 3.5, 'red'),
    ("d(His Nδ1 → Asp Oδ)", primary_data.get("d3"), 3.0, 'green'),
]
for ax, (label, d_vals, thresh, color) in zip(axes, labels_thresholds):
    if d_vals is not None:
        ax.plot(time_ns[:len(d_vals)], d_vals, color=color, linewidth=0.4, alpha=0.5)
        ax.axhline(y=thresh, color='red', linestyle='--', label=f'{thresh} Å threshold')
        ax.set_ylabel(f'{label} (Å)', fontsize=10)
        mean_v = np.mean(d_vals)
        ax.set_title(f'{label}: mean = {mean_v:.2f} Å', fontsize=11)
        ax.legend(fontsize=9)

axes[-1].set_xlabel('Time (ns)', fontsize=12)
plt.suptitle('EXP-15: Catalytic Triad Geometry (BPTI–Trypsin)', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-15_triad_distances.png", dpi=300)
plt.close(fig)
```

### Figure 4: Block Average Convergence

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
for name, data in results_all.items():
    blocks = data["block_means"]
    ax.plot(range(1, len(blocks)+1), blocks, 'o-', label=name, markersize=8)
ax.axhspan(2.4, 3.0, alpha=0.1, color='green', label='PASS range')
ax.axhline(y=2.7, color='green', linestyle='--')
ax.set_xlabel('Block (10 ns each)', fontsize=12)
ax.set_ylabel('Mean d(Ser Oγ, P1 C) (Å)', fontsize=12)
ax.set_title('EXP-15: Block Convergence', fontsize=13)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-15_block_convergence.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-15: Catalytic Contact Distance — Results Report

**Experiment ID:** EXP-15  **Feature ID:** F-15  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results

| System | d(Ser Oγ, P1 C) (Å) | SE | >3.5 Å (%) | >4.0 Å (%) |
|--------|---------------------|----|------------|------------|
| BPTI-trypsin | [val] | [SE] | [pct] | [pct] |
| SPINK7-KLK5 | [val] | [SE] | [pct] | [pct] |

Crystal reference (2PTC): [val] Å
Benchmark: 2.7 ± 0.3 Å (95% CI [2.4, 3.0])

## Figures
1. Catalytic distance time series
2. Distance histogram
3. Catalytic triad geometry
4. Block convergence

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
