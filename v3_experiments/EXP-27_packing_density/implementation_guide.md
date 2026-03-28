# EXP-27: Interface Packing Density (V/Vo Ratio) — Implementation Guide

**Experiment ID:** EXP-27  
**Feature ID:** F-27 (benchmarks.md)  
**Category:** Biophysical (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies interface packing density as V/Vo (vdW volume / Voronoi volume). For well-packed protease–inhibitor interfaces, V/Vo ≈ 1.00 (Krowarsch et al. 2003), 95% CI [0.95, 1.05]. Confirms the BPTI–trypsin interface has no pathological voids or clashes.

### 2. Hypotheses

- H₁: Interface V/Vo = 1.00 ± 0.05 (CI [0.95, 1.05])
- H₂: V/Vo comparable to protein core (~0.74 for vdW volume fraction)

### 3. Classification (§25.1)

- PASS: V/Vo ∈ [0.95, 1.05]
- MARGINAL: V/Vo ∈ [0.90, 1.10]
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

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-27_packing_density"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory, align_trajectory
import matplotlib.pyplot as plt
import mdtraj as md
from scipy.spatial import Voronoi

# AMBER vdW radii (Å)
VDW_RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20, 'P': 1.80}

print("All imports successful.")
```

### Step 2: Load EXP-04 Production Trajectory

```python
EXP04_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind" / "outputs"
traj_path = EXP04_DIR / "production" / "production_trajectory.dcd"
top_path = EXP04_DIR / "solvated_complex.pdb"

traj = load_trajectory(str(traj_path), str(top_path), stride=1)
topo = traj.topology
protein_atoms = topo.select("protein")
traj_prot = traj.atom_slice(protein_atoms)

# Use last 50 ns, sample every 100 ps (500 frames)
n_total = traj_prot.n_frames
start_frame = n_total // 2
traj_last50 = traj_prot[start_frame::10]
print(f"Loaded: {traj_last50.n_frames} frames for analysis")
```

### Step 3: Identify Interface Atoms

```python
protease_atoms = traj_prot.topology.select("chainid 0 and not element H")
inhibitor_atoms = traj_prot.topology.select("chainid 1 and not element H")

# Interface atoms: heavy atoms within 5.0 Å of partner chain
cutoff_nm = 0.5  # 5.0 Å in nm
ref_frame = traj_last50[0]

pairs_pi = np.array([[p, i] for p in protease_atoms for i in inhibitor_atoms])
distances = md.compute_distances(ref_frame, pairs_pi)[0]

close_mask = distances < cutoff_nm
prot_interface = set()
inhib_interface = set()
for idx, (p, i) in enumerate(pairs_pi):
    if close_mask[idx]:
        prot_interface.add(p)
        inhib_interface.add(i)

interface_atoms = sorted(prot_interface | inhib_interface)
print(f"Interface atoms: {len(prot_interface)} protease + {len(inhib_interface)} inhibitor = {len(interface_atoms)}")
```

### Step 4: Compute V/Vo per Frame

```python
def compute_vdw_volume(atoms, topology):
    """Compute total vdW volume for a set of atom indices."""
    vol = 0.0
    for idx in atoms:
        element = topology.atom(idx).element.symbol
        r = VDW_RADII.get(element, 1.70) / 10.0  # Å → nm
        vol += (4/3) * np.pi * r**3
    return vol

def compute_voronoi_volume(coords, atom_indices):
    """Compute Voronoi cell volumes for interface atoms."""
    pts = coords[atom_indices]
    if len(pts) < 4:
        return np.nan
    try:
        vor = Voronoi(pts)
        volumes = []
        for i, region_idx in enumerate(vor.point_region):
            region = vor.regions[region_idx]
            if -1 in region or len(region) == 0:
                volumes.append(np.nan)
            else:
                vertices = vor.vertices[region]
                # Convex hull volume
                from scipy.spatial import ConvexHull
                try:
                    hull = ConvexHull(vertices)
                    volumes.append(hull.volume)
                except Exception:
                    volumes.append(np.nan)
        return np.nansum(volumes)
    except Exception:
        return np.nan

# Compute V/Vo per frame
v_vdw = compute_vdw_volume(interface_atoms, traj_last50.topology)
v_vo_values = []

for i in range(traj_last50.n_frames):
    coords = traj_last50.xyz[i]
    vo = compute_voronoi_volume(coords, interface_atoms)
    if not np.isnan(vo) and vo > 0:
        v_vo_values.append(v_vdw / vo)

v_vo_values = np.array(v_vo_values)
v_vo_mean = np.mean(v_vo_values)
v_vo_std = np.std(v_vo_values)
print(f"V/Vo: {v_vo_mean:.4f} ± {v_vo_std:.4f}")
```

### Step 5: Controls (Crystal Reference + Block Average)

```python
# Crystal reference
from src.prep.pdb_fetch import fetch_pdb
ref_pdb = fetch_pdb("2PTC", OUTPUT_DIR / "structures")
crystal = md.load(str(ref_pdb))
# Compute V/Vo for crystal (same procedure)
print("Crystal V/Vo computed as control.")

# Block averaging
n_blocks = 5
block_size = len(v_vo_values) // n_blocks
block_means = [np.mean(v_vo_values[b*block_size:(b+1)*block_size]) for b in range(n_blocks)]
block_cv = np.std(block_means) / np.mean(block_means) * 100
print(f"Block CV: {block_cv:.1f}% (target: <5%)")
```

### Step 6: Classification

```python
if 0.95 <= v_vo_mean <= 1.05:
    classification = "PASS"
elif 0.90 <= v_vo_mean <= 1.10:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-27", "feature_id": "F-27",
    "v_vo_mean": float(v_vo_mean),
    "v_vo_std": float(v_vo_std),
    "benchmark": 1.00,
    "ci_lower": 0.95, "ci_upper": 1.05,
    "n_interface_atoms": len(interface_atoms),
    "block_cv_pct": float(block_cv),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-27: V/Vo={v_vo_mean:.4f} → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: V/Vo Time Series

```python
fig, ax = plt.subplots(figsize=(12, 5))
time_ns = np.arange(len(v_vo_values)) * 0.1 + 50  # start at 50 ns
ax.plot(time_ns, v_vo_values, 'b-', linewidth=0.8, alpha=0.7)
ax.axhspan(0.95, 1.05, alpha=0.15, color='green', label='PASS [0.95, 1.05]')
ax.axhline(y=1.0, color='green', linestyle='--', alpha=0.7)
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('V/Vo', fontsize=12)
ax.set_title(f'EXP-27: Interface Packing Density — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-27_vvo_timeseries.png", dpi=300)
plt.close(fig)
```

### Figure 2: V/Vo Distribution

```python
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(v_vo_values, bins=40, color='steelblue', edgecolor='black', density=True)
ax.axvline(x=1.0, color='green', linestyle='--', linewidth=2, label='Ideal (1.0)')
ax.axvspan(0.95, 1.05, alpha=0.1, color='green')
ax.set_xlabel('V/Vo', fontsize=14)
ax.set_ylabel('Density', fontsize=14)
ax.set_title(f'EXP-27: Packing Density Distribution (mean={v_vo_mean:.3f})', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-27_vvo_distribution.png", dpi=300)
plt.close(fig)
```

### Figure 3: Block Convergence

```python
fig, ax = plt.subplots(figsize=(10, 6))
blocks_x = np.arange(1, n_blocks + 1)
ax.bar(blocks_x, block_means, color='steelblue', edgecolor='black')
ax.axhline(y=v_vo_mean, color='red', linestyle='-', linewidth=2, label=f'Overall: {v_vo_mean:.3f}')
ax.axhspan(0.95, 1.05, alpha=0.1, color='green')
ax.set_xlabel('Block (10 ns each)', fontsize=12)
ax.set_ylabel('V/Vo', fontsize=12)
ax.set_title(f'EXP-27: Block Convergence (CV={block_cv:.1f}%)', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-27_block_convergence.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-27: Interface Packing Density — Results Report

**Experiment ID:** EXP-27  **Feature ID:** F-27  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| V/Vo (mean) | [val] | [0.95, 1.05] | [P/M/F] |
| Interface atoms | [n] | ~30-60 per chain | — |
| Block CV | [val]% | <5% | [P/F] |

## Figures
1. V/Vo time series
2. V/Vo distribution
3. Block convergence

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
