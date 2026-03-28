# EXP-18: Interfacial Water Molecule Count — Implementation Guide

**Experiment ID:** EXP-18  
**Feature ID:** F-18 (benchmarks.md)  
**Category:** Structural (Semi-Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies ordered water molecules at the protease–inhibitor interface. BPTI–trypsin: ~15 crystallographic waters within 3.5 Å of both partners (95% CI [10, 20]). Uses EXP-04 trajectory. Includes 3D water density mapping, residence time analysis, and comparison with crystal water positions.

### 2. Classification (§25.1)

- PASS: Mean count in [10, 20]
- MARGINAL: Count in [7, 25]
- FAIL: < 7 or > 25

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-18_interfacial_water_count"
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

# Use last 50 ns
n_skip = traj.n_frames // 2
traj_analysis = traj[n_skip:]

protease_heavy = topo.select("chainid 0 and not water and element != H")
inhibitor_heavy = topo.select("chainid 1 and not water and element != H")
water_oxygen = topo.select("water and name O")
print(f"Protease heavy: {len(protease_heavy)}, Inhibitor heavy: {len(inhibitor_heavy)}, "
      f"Water O: {len(water_oxygen)}, Frames: {traj_analysis.n_frames}")
```

### Step 3: Crystal Structure Water Count

```python
from src.prep.pdb_fetch import fetch_pdb
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)
pdb_2ptc = fetch_pdb("2PTC", data_dir)
crystal = md.load(str(pdb_2ptc))
ctopo = crystal.topology

crystal_water_o = ctopo.select("water and name O")
crystal_prot = ctopo.select("chainid 0 and not water and element != H")
crystal_inh = ctopo.select("chainid 1 and not water and element != H")

n_crystal_interfacial = 0
crystal_water_positions = []
for wo in crystal_water_o:
    # Distance to closest protease heavy atom
    d_prot = md.compute_distances(crystal, [[wo, pa] for pa in crystal_prot])[0]
    # Distance to closest inhibitor heavy atom
    d_inh = md.compute_distances(crystal, [[wo, ia] for ia in crystal_inh])[0]
    if np.min(d_prot) < 0.35 and np.min(d_inh) < 0.35:  # 3.5 Å in nm
        n_crystal_interfacial += 1
        crystal_water_positions.append(crystal.xyz[0, wo])

print(f"Crystal interfacial waters: {n_crystal_interfacial}")
crystal_water_positions = np.array(crystal_water_positions) if crystal_water_positions else np.array([])
```

### Step 4: Per-Frame Interfacial Water Count

```python
water_counts = []
cutoff = 0.35  # 3.5 Å in nm

for i in range(traj_analysis.n_frames):
    frame = traj_analysis[i]
    xyz = frame.xyz[0]

    n_interfacial = 0
    for wo in water_oxygen:
        pos_w = xyz[wo]
        # Min distance to protease
        d_prot = np.sqrt(np.min(np.sum((xyz[protease_heavy] - pos_w)**2, axis=1)))
        if d_prot > cutoff:
            continue
        # Min distance to inhibitor
        d_inh = np.sqrt(np.min(np.sum((xyz[inhibitor_heavy] - pos_w)**2, axis=1)))
        if d_inh <= cutoff:
            n_interfacial += 1

    water_counts.append(n_interfacial)
    if (i + 1) % 50 == 0:
        print(f"Frame {i+1}/{traj_analysis.n_frames}: {n_interfacial} interfacial waters")

water_counts = np.array(water_counts)
mean_count = np.mean(water_counts)
std_count = np.std(water_counts)
print(f"Interfacial waters: {mean_count:.1f} ± {std_count:.1f}")
```

### Step 5: Water Density Map (Optional Advanced)

```python
# 3D histogram of water oxygen positions near interface
# Define grid centered on interface COM
interface_com = 0.5 * (np.mean(traj_analysis.xyz[0][protease_heavy], axis=0) +
                        np.mean(traj_analysis.xyz[0][inhibitor_heavy], axis=0))

grid_range = 1.5  # nm, ±15 Å
n_bins = 30
bins = [np.linspace(interface_com[d] - grid_range, interface_com[d] + grid_range, n_bins+1)
        for d in range(3)]

# Accumulate water oxygen positions near interface
water_positions_all = []
for i in range(traj_analysis.n_frames):
    xyz = traj_analysis.xyz[i]
    for wo in water_oxygen:
        pos_w = xyz[wo]
        d_prot = np.sqrt(np.min(np.sum((xyz[protease_heavy] - pos_w)**2, axis=1)))
        d_inh = np.sqrt(np.min(np.sum((xyz[inhibitor_heavy] - pos_w)**2, axis=1)))
        if d_prot < cutoff * 2 and d_inh < cutoff * 2:  # within 7 Å of both
            water_positions_all.append(pos_w)

water_positions_all = np.array(water_positions_all)
density, edges = np.histogramdd(water_positions_all, bins=bins)
density /= traj_analysis.n_frames  # Average per frame

# Bulk density reference
box_vol = np.prod([b[-1] - b[0] for b in bins]) * 1e-27  # nm³
n_water_bulk = 33.3 * box_vol * 1e27  # ~33.3 waters/nm³ for TIP3P
density_ratio = density / (n_water_bulk / n_bins**3) if n_water_bulk > 0 else density

print(f"Max density ratio: {np.max(density_ratio):.1f}× bulk")
print(f"Sites >2× bulk: {np.sum(density_ratio > 2)}")
```

### Step 6: Classification

```python
if 10 <= mean_count <= 20:
    classification = "PASS"
elif 7 <= mean_count <= 25:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-18", "feature_id": "F-18",
    "mean_interfacial_waters": float(mean_count),
    "std_interfacial_waters": float(std_count),
    "crystal_interfacial_waters": n_crystal_interfacial,
    "benchmark": "15 ± 5",
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-18: {mean_count:.1f} interfacial waters → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Interfacial Water Count Time Series

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 6))
time_ns = np.arange(len(water_counts)) * 0.1 + 50
ax.plot(time_ns, water_counts, 'b-', linewidth=0.5, alpha=0.7)
ax.axhline(y=15, color='green', linestyle='--', linewidth=2, label='Crystal (15)')
ax.axhspan(10, 20, alpha=0.1, color='green', label='PASS range')
ax.axhline(y=mean_count, color='red', linestyle='-', linewidth=1.5,
           label=f'Mean = {mean_count:.1f}')
ax.set_xlabel('Time (ns)', fontsize=14)
ax.set_ylabel('Interfacial Water Count', fontsize=14)
ax.set_title(f'EXP-18: Interfacial Waters — {classification}', fontsize=15)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-18_water_count_timeseries.png", dpi=300)
plt.close(fig)
```

### Figure 2: Water Count Distribution

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(water_counts, bins=range(0, max(water_counts)+2), color='cyan',
        edgecolor='black', density=True, align='left')
ax.axvline(x=15, color='green', linestyle='--', linewidth=2, label='Crystal (15)')
ax.axvspan(10, 20, alpha=0.1, color='green', label='PASS range')
ax.set_xlabel('Interfacial Water Count', fontsize=14)
ax.set_ylabel('Probability', fontsize=14)
ax.set_title('EXP-18: Interfacial Water Distribution', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-18_water_distribution.png", dpi=300)
plt.close(fig)
```

### Figure 3: Summary Dashboard

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
# Left: comparison bar
ax1.bar(['Crystal\n(PDB 2PTC)', 'MD Mean'], [n_crystal_interfacial, mean_count],
        color=['steelblue', 'coral'], edgecolor='black', yerr=[0, std_count], capsize=5)
ax1.axhspan(10, 20, alpha=0.1, color='green', label='PASS range')
ax1.set_ylabel('Interfacial Water Count', fontsize=12)
ax1.set_title('Crystal vs MD Comparison', fontsize=13)
ax1.legend()

# Right: classification box
summary = (f"Crystal: {n_crystal_interfacial} waters\n"
           f"MD: {mean_count:.1f} ± {std_count:.1f}\n"
           f"Benchmark: 15 ± 5\n"
           f"Classification: {classification}")
color = 'lightgreen' if classification == 'PASS' else 'lightyellow'
ax2.text(0.5, 0.5, summary, transform=ax2.transAxes, fontsize=14,
         va='center', ha='center', bbox=dict(boxstyle='round', facecolor=color))
ax2.axis('off')
ax2.set_title('Summary', fontsize=13)

plt.suptitle('EXP-18: Interfacial Water Count', fontsize=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-18_summary.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-18: Interfacial Water Count — Results Report

**Experiment ID:** EXP-18  **Feature ID:** F-18  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value |
|--------|-------|
| MD mean interfacial waters | [val] ± [std] |
| Crystal interfacial waters | [val] |
| Benchmark | 15 ± 5 |
| High-density sites (>2× bulk) | [n] |

## Figures
1. Water count time series
2. Count distribution
3. Summary dashboard

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
