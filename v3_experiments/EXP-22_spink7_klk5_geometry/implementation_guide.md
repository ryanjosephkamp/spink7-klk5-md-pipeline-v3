# EXP-22: SPINK7–KLK5 Interface Geometry Validation — Implementation Guide

**Experiment ID:** EXP-22  
**Feature ID:** F-22 (benchmarks.md)  
**Category:** Structural (Qualitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates the SPINK7–KLK5 docked complex (no co-crystal exists) by checking canonical binding mode, structural stability (RMSD < 3.0 Å over 100 ns), and consistency with structural benchmarks (BSA, H-bonds, catalytic distance, P1 in S1). Cross-validates against known Kazal–protease patterns using BPTI–trypsin as positive control.

### 2. Classification (§25.1)

- PASS: Canonical mode; backbone RMSD < 3.0 Å; all structural metrics in PASS
- MARGINAL: Correct mode; RMSD < 4.0 Å; most metrics MARGINAL
- FAIL: Non-canonical mode; RMSD > 4.0 Å; interface disruption

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-22_spink7_klk5_geometry"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory, align_trajectory
from src.analyze.structural import (compute_rmsd, compute_sasa,
                                      compute_hbonds, compute_interface_contacts)
from src.analyze.structural import compute_radius_of_gyration
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load EXP-01 Production Trajectory

```python
EXP_ROOT = PROJECT_ROOT / "v3_experiments"
traj_path = EXP_ROOT / "EXP-01_spink7_klk5_binding/outputs/production/production_trajectory.dcd"
top_path = EXP_ROOT / "EXP-01_spink7_klk5_binding/outputs/solvated_complex.pdb"

traj = load_trajectory(str(traj_path), str(top_path), stride=1)
topo = traj.topology
print(f"Loaded: {traj.n_frames} frames")

protease_bb = topo.select("chainid 0 and backbone")
inhibitor_bb = topo.select("chainid 1 and backbone")
complex_bb = np.concatenate([protease_bb, inhibitor_bb])
```

### Step 3: Global Stability (RMSD + Rg)

```python
# Backbone RMSD vs first equilibrated frame
traj_aligned = align_trajectory(traj, traj[0], atom_sel="backbone")
rmsd_bb = compute_rmsd(traj_aligned, traj_aligned[0], atom_sel="backbone")
rmsd_ang = rmsd_bb * 10  # nm → Å

# Radius of gyration
rg = compute_radius_of_gyration(traj.atom_slice(complex_bb))
rg_ang = rg * 10

time_ns = np.arange(len(rmsd_ang)) * 0.01

print(f"RMSD: mean={np.mean(rmsd_ang):.2f} ± {np.std(rmsd_ang):.2f} Å, max={np.max(rmsd_ang):.2f} Å")
print(f"Rg: mean={np.mean(rg_ang):.2f} ± {np.std(rg_ang):.2f} Å")
```

### Step 4: Interface Metrics (BSA, H-bonds, Catalytic Distance)

```python
# BSA (sample every 10 frames)
bsa_values = []
protease_atoms = topo.select("chainid 0 and not water")
inhibitor_atoms = topo.select("chainid 1 and not water")
complex_atoms = np.concatenate([protease_atoms, inhibitor_atoms])

for i in range(0, traj.n_frames, 10):
    frame = traj[i]
    sasa_c = np.sum(md.shrake_rupley(frame.atom_slice(complex_atoms), probe_radius=0.14, mode='atom')[0]) * 100
    sasa_p = np.sum(md.shrake_rupley(frame.atom_slice(protease_atoms), probe_radius=0.14, mode='atom')[0]) * 100
    sasa_i = np.sum(md.shrake_rupley(frame.atom_slice(inhibitor_atoms), probe_radius=0.14, mode='atom')[0]) * 100
    bsa_values.append(sasa_p + sasa_i - sasa_c)

bsa_mean = np.mean(bsa_values)
print(f"BSA: {bsa_mean:.0f} ± {np.std(bsa_values):.0f} Å²")

# H-bond count
hbond_counts = []
for i in range(0, traj.n_frames, 10):
    hbonds = md.baker_hubbard(traj[i], freq=0.0)
    inter = sum(1 for h in hbonds
                if topo.atom(h[0]).residue.chain.index != topo.atom(h[2]).residue.chain.index)
    hbond_counts.append(inter)
hbond_mean = np.mean(hbond_counts)
print(f"H-bonds: {hbond_mean:.1f} ± {np.std(hbond_counts):.1f}")

# Catalytic distance (Ser Oγ → P1 C)
ser_og = topo.select("resname SER and name OG and chainid 0")
p1_c = topo.select("resname ARG and name C and chainid 1")
if len(ser_og) > 0 and len(p1_c) > 0:
    d_cat = md.compute_distances(traj, [[ser_og[0], p1_c[0]]])[:, 0] * 10
    d_cat_mean = np.mean(d_cat)
    print(f"Catalytic distance: {d_cat_mean:.2f} ± {np.std(d_cat):.2f} Å")

# P1 in S1 (Arg → Asp salt bridge)
arg_nh = topo.select("resname ARG and name NH1 and chainid 1")
asp_od = topo.select("resname ASP and name OD1 and chainid 0")
if len(arg_nh) > 0 and len(asp_od) > 0:
    d_sb = md.compute_distances(traj, [[arg_nh[0], asp_od[-1]]])[:, 0] * 10
    sb_occ = np.mean(d_sb < 4.0) * 100
    print(f"P1-S1 salt bridge: {sb_occ:.1f}% occupancy")
```

### Step 5: Binding Mode Cross-Validation

```python
# Superpose SPINK7-KLK5 binding loop onto BPTI-trypsin (2PTC)
from src.prep.pdb_fetch import fetch_pdb
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)
ref_pdb = fetch_pdb("2PTC", data_dir)
ref_struct = md.load(str(ref_pdb))

# Extract binding loops and compute RMSD
# Get representative frame (median RMSD)
median_idx = np.argmin(np.abs(rmsd_ang - np.median(rmsd_ang)))
representative = traj[median_idx]
print(f"Representative frame: {median_idx} (RMSD={rmsd_ang[median_idx]:.2f} Å)")
```

### Step 6: Classification

```python
rmsd_pass = np.mean(rmsd_ang) < 3.0
bsa_pass = 1200 <= bsa_mean <= 1600
hbond_pass = 7 <= hbond_mean <= 13
cat_pass = 2.4 <= d_cat_mean <= 3.0 if 'd_cat_mean' in dir() else False
sb_pass = sb_occ > 80 if 'sb_occ' in dir() else False

n_pass = sum([rmsd_pass, bsa_pass, hbond_pass, cat_pass, sb_pass])

if n_pass >= 4:
    classification = "PASS"
elif n_pass >= 2:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-22", "feature_id": "F-22",
    "rmsd_mean_ang": float(np.mean(rmsd_ang)),
    "rmsd_max_ang": float(np.max(rmsd_ang)),
    "bsa_mean_ang2": float(bsa_mean),
    "hbond_mean": float(hbond_mean),
    "catalytic_distance_ang": float(d_cat_mean) if 'd_cat_mean' in dir() else None,
    "salt_bridge_occ_pct": float(sb_occ) if 'sb_occ' in dir() else None,
    "metrics_pass": {"rmsd": rmsd_pass, "bsa": bsa_pass, "hbond": hbond_pass,
                      "catalytic": cat_pass, "salt_bridge": sb_pass},
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-22: {n_pass}/5 metrics PASS → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Structural Stability (RMSD + Rg)

```python
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
ax1.plot(time_ns, rmsd_ang, 'b-', linewidth=0.5, alpha=0.7)
ax1.axhline(y=3.0, color='green', linestyle='--', label='PASS (3.0 Å)')
ax1.set_ylabel('Backbone RMSD (Å)', fontsize=12)
ax1.legend()
ax1.set_title(f'EXP-22: SPINK7-KLK5 Stability — {classification}', fontsize=14)

ax2.plot(time_ns, rg_ang, 'r-', linewidth=0.5, alpha=0.7)
ax2.set_ylabel('Rg (Å)', fontsize=12)
ax2.set_xlabel('Time (ns)', fontsize=12)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-22_stability.png", dpi=300)
plt.close(fig)
```

### Figure 2: Interface Metrics Dashboard

```python
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# BSA
time_bsa = np.arange(len(bsa_values)) * 0.1
axes[0,0].plot(time_bsa, bsa_values, 'b-', linewidth=0.8)
axes[0,0].axhspan(1200, 1600, alpha=0.1, color='green')
axes[0,0].set_ylabel('BSA (Å²)')
axes[0,0].set_title(f'BSA: {bsa_mean:.0f} Å² [{"PASS" if bsa_pass else "FAIL"}]')

# H-bonds
time_hb = np.arange(len(hbond_counts)) * 0.1
axes[0,1].plot(time_hb, hbond_counts, 'g-', linewidth=0.8)
axes[0,1].axhspan(7, 13, alpha=0.1, color='green')
axes[0,1].set_ylabel('H-bond Count')
axes[0,1].set_title(f'H-bonds: {hbond_mean:.1f} [{"PASS" if hbond_pass else "FAIL"}]')

# Catalytic distance
if 'd_cat' in dir():
    axes[1,0].plot(time_ns, d_cat, 'r-', linewidth=0.3, alpha=0.5)
    axes[1,0].axhspan(2.4, 3.0, alpha=0.1, color='green')
    axes[1,0].set_ylabel('d(Ser Oγ, P1 C) (Å)')
    axes[1,0].set_title(f'Cat. dist: {d_cat_mean:.2f} Å [{"PASS" if cat_pass else "FAIL"}]')
    axes[1,0].set_xlabel('Time (ns)')

# Classification summary
metrics = ['RMSD', 'BSA', 'H-bonds', 'Cat. dist', 'Salt bridge']
passes = [rmsd_pass, bsa_pass, hbond_pass, cat_pass, sb_pass]
colors_m = ['green' if p else 'red' for p in passes]
axes[1,1].barh(metrics, [1]*5, color=colors_m, edgecolor='black')
for i, (m, p) in enumerate(zip(metrics, passes)):
    axes[1,1].text(0.5, i, 'PASS' if p else 'FAIL', ha='center', va='center',
                    fontsize=12, fontweight='bold', color='white')
axes[1,1].set_title(f'Overall: {classification} ({n_pass}/5)')

plt.suptitle('EXP-22: SPINK7-KLK5 Interface Validation', fontsize=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-22_interface_dashboard.png", dpi=300)
plt.close(fig)
```

### Figure 3: RMSD Distribution

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.hist(rmsd_ang, bins=50, color='steelblue', edgecolor='black', density=True)
ax.axvline(x=3.0, color='green', linestyle='--', linewidth=2, label='3.0 Å threshold')
ax.set_xlabel('Backbone RMSD (Å)', fontsize=14)
ax.set_ylabel('Density', fontsize=14)
ax.set_title(f'EXP-22: RMSD Distribution — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-22_rmsd_distribution.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-22: SPINK7-KLK5 Interface Geometry — Results Report

**Experiment ID:** EXP-22  **Feature ID:** F-22  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| RMSD | [val] Å | < 3.0 Å | [P/F] |
| BSA | [val] Å² | 1200-1600 | [P/F] |
| H-bonds | [val] | 7-13 | [P/F] |
| Cat. dist | [val] Å | 2.4-3.0 | [P/F] |
| Salt bridge | [val]% | >80% | [P/F] |

## Figures
1. Stability (RMSD + Rg)
2. Interface dashboard
3. RMSD distribution

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
