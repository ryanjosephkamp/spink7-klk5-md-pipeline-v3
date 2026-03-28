# EXP-20: Binding Loop Conservation Across Kazal Inhibitors — Implementation Guide

**Experiment ID:** EXP-20  
**Feature ID:** F-20 (benchmarks.md)  
**Category:** Structural (Qualitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates that the canonical binding loop (P3–P3') is structurally conserved across Kazal-family inhibitors despite low sequence identity (20–40%). Pairwise backbone RMSD < 1.0 Å for crystal structures, < 1.5 Å for MD-averaged (Laskowski & Qasim 2000). Tests BPTI, SPINK1/PSTI, SPINK7, and LEKTI D6.

### 2. Classification (§25.1)

- PASS: All pairwise RMSD < 1.0 Å (crystal), < 1.5 Å (MD)
- MARGINAL: All < 1.5 Å (crystal), < 2.0 Å (MD)
- FAIL: Any pair > 2.0 Å

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path
from itertools import combinations

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-20_loop_conservation"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory, align_trajectory
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Extract Binding Loops from Crystal/Starting Structures

```python
from src.prep.pdb_fetch import fetch_pdb
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

# Define inhibitor binding loop residue ranges (P3-P3', 7 residues)
# These are approximate and need verification per structure
INHIBITORS = {
    "BPTI": {
        "source": "EXP-04_bpti_trypsin_dg_bind",
        "chain": 1,
        "p1_offset": 0,  # Will find by residue name
        "p1_resname": "LYS",
    },
    "SPINK1": {
        "source": "EXP-06_spink1_trypsin_dg_bind",
        "chain": 1,
        "p1_resname": "LYS",
    },
    "SPINK7": {
        "source": "EXP-01_spink7_klk5_binding",
        "chain": 1,
        "p1_resname": "ARG",
    },
}

EXP_ROOT = PROJECT_ROOT / "v3_experiments"

loop_structures = {}
for name, info in INHIBITORS.items():
    top_path = EXP_ROOT / info["source"] / "outputs/solvated_complex.pdb"
    if top_path.exists():
        struct = md.load(str(top_path))
        topo = struct.topology

        # Find P1 residue on inhibitor chain
        inh_residues = [r for r in topo.residues if r.chain.index == info["chain"]]
        p1_res = [r for r in inh_residues if r.name == info["p1_resname"]]

        if p1_res:
            p1_idx = p1_res[0].index
            loop_resids = list(range(p1_idx - 3, p1_idx + 4))
            loop_bb = topo.select(f"resid {' '.join(str(r) for r in loop_resids)} and name CA C N O")

            if len(loop_bb) > 0:
                loop_structures[name] = {
                    "struct": struct,
                    "loop_atoms": loop_bb,
                    "p1_idx": p1_idx,
                    "loop_resids": loop_resids,
                }
                print(f"{name}: loop = resids {loop_resids}, {len(loop_bb)} backbone atoms")
    else:
        print(f"WARNING: {name} structure not found")
```

### Step 3: Pairwise Crystal RMSD

```python
crystal_rmsd_matrix = {}

inhibitor_names = list(loop_structures.keys())
for i, j in combinations(range(len(inhibitor_names)), 2):
    name_i = inhibitor_names[i]
    name_j = inhibitor_names[j]

    struct_i = loop_structures[name_i]["struct"]
    struct_j = loop_structures[name_j]["struct"]
    atoms_i = loop_structures[name_i]["loop_atoms"]
    atoms_j = loop_structures[name_j]["loop_atoms"]

    # Both loops must have same number of backbone atoms
    n_atoms = min(len(atoms_i), len(atoms_j))
    loop_i = struct_i.atom_slice(atoms_i[:n_atoms])
    loop_j = struct_j.atom_slice(atoms_j[:n_atoms])

    # Superpose j onto i
    loop_j_sup = loop_j.superpose(loop_i)
    rmsd = md.rmsd(loop_j_sup, loop_i, 0)[0] * 10  # Å

    pair_key = f"{name_i}-{name_j}"
    crystal_rmsd_matrix[pair_key] = float(rmsd)
    print(f"Crystal {pair_key}: RMSD = {rmsd:.2f} Å")
```

### Step 4: MD-Averaged Loop RMSD

```python
md_rmsd_matrix = {}

# Use BPTI as reference
ref_name = "BPTI"
ref_struct = loop_structures[ref_name]["struct"]
ref_atoms = loop_structures[ref_name]["loop_atoms"]

md_rmsd_per_system = {}
for name, info_dict in loop_structures.items():
    traj_path = EXP_ROOT / INHIBITORS[name]["source"] / "outputs/production/production_trajectory.dcd"
    top_path = EXP_ROOT / INHIBITORS[name]["source"] / "outputs/solvated_complex.pdb"

    if traj_path.exists():
        traj = load_trajectory(str(traj_path), str(top_path), stride=10)
        n_skip = traj.n_frames // 2
        traj_analysis = traj[n_skip:]

        loop_atoms = info_dict["loop_atoms"]
        n_atoms = min(len(loop_atoms), len(ref_atoms))

        rmsd_vals = []
        for frame_idx in range(traj_analysis.n_frames):
            frame_loop = traj_analysis[frame_idx].atom_slice(loop_atoms[:n_atoms])
            ref_loop = ref_struct.atom_slice(ref_atoms[:n_atoms])
            frame_sup = frame_loop.superpose(ref_loop)
            r = md.rmsd(frame_sup, ref_loop, 0)[0] * 10
            rmsd_vals.append(r)

        md_rmsd_per_system[name] = np.array(rmsd_vals)
        print(f"MD {name} vs BPTI ref: {np.mean(rmsd_vals):.2f} ± {np.std(rmsd_vals):.2f} Å")
```

### Step 5: Classification

```python
max_crystal_rmsd = max(crystal_rmsd_matrix.values()) if crystal_rmsd_matrix else 999
max_md_rmsd = max(np.mean(v) for v in md_rmsd_per_system.values()) if md_rmsd_per_system else 999

if max_crystal_rmsd < 1.0 and max_md_rmsd < 1.5:
    classification = "PASS"
elif max_crystal_rmsd < 1.5 and max_md_rmsd < 2.0:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-20", "feature_id": "F-20",
    "crystal_pairwise_rmsd": crystal_rmsd_matrix,
    "md_mean_rmsd_vs_bpti": {n: float(np.mean(v)) for n, v in md_rmsd_per_system.items()},
    "max_crystal_rmsd": float(max_crystal_rmsd),
    "max_md_rmsd": float(max_md_rmsd),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-20: max crystal RMSD = {max_crystal_rmsd:.2f} Å, max MD = {max_md_rmsd:.2f} Å → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Pairwise RMSD Heatmap

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 7))
names = inhibitor_names
n = len(names)
matrix = np.zeros((n, n))
for (i, j) in combinations(range(n), 2):
    key = f"{names[i]}-{names[j]}"
    val = crystal_rmsd_matrix.get(key, 0)
    matrix[i, j] = val
    matrix[j, i] = val

im = ax.imshow(matrix, cmap='RdYlGn_r', vmin=0, vmax=2.0)
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(names, fontsize=11, rotation=45, ha='right')
ax.set_yticklabels(names, fontsize=11)
for i in range(n):
    for j in range(n):
        ax.text(j, i, f'{matrix[i,j]:.2f}', ha='center', va='center', fontsize=12,
                color='white' if matrix[i,j] > 1.5 else 'black')
plt.colorbar(im, ax=ax, label='RMSD (Å)')
ax.set_title(f'EXP-20: Binding Loop Pairwise RMSD — {classification}', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-20_rmsd_heatmap.png", dpi=300)
plt.close(fig)
```

### Figure 2: MD RMSD vs Reference

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 6))
for name, vals in md_rmsd_per_system.items():
    time_ns = np.arange(len(vals)) * 0.1 + 50
    ax.plot(time_ns, vals, linewidth=0.8, label=f'{name} (mean={np.mean(vals):.2f} Å)')
ax.axhline(y=1.5, color='green', linestyle='--', label='PASS threshold (1.5 Å)')
ax.set_xlabel('Time (ns)', fontsize=14)
ax.set_ylabel('Loop RMSD vs BPTI ref (Å)', fontsize=14)
ax.set_title('EXP-20: Loop Conservation During MD', fontsize=14)
ax.legend(fontsize=10)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-20_md_rmsd.png", dpi=300)
plt.close(fig)
```

### Figure 3: Summary Bar Chart

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
systems = list(md_rmsd_per_system.keys())
means = [np.mean(md_rmsd_per_system[s]) for s in systems]
stds = [np.std(md_rmsd_per_system[s]) for s in systems]
colors = ['green' if m < 1.5 else 'orange' if m < 2.0 else 'red' for m in means]
ax.bar(systems, means, yerr=stds, color=colors, edgecolor='black', capsize=5)
ax.axhline(y=1.0, color='green', linestyle='--', label='Crystal PASS (1.0 Å)')
ax.axhline(y=1.5, color='orange', linestyle=':', label='MD PASS (1.5 Å)')
ax.set_ylabel('Mean Loop RMSD vs BPTI (Å)', fontsize=12)
ax.set_title(f'EXP-20: Loop Conservation Summary — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-20_summary.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-20: Binding Loop Conservation — Results Report

**Experiment ID:** EXP-20  **Feature ID:** F-20  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Crystal Pairwise RMSD
| Pair | RMSD (Å) |
|------|----------|
| BPTI-SPINK1 | [val] |
| BPTI-SPINK7 | [val] |
| SPINK1-SPINK7 | [val] |

## MD Average RMSD vs BPTI Reference
| System | Mean RMSD (Å) | SD |
|--------|---------------|----|
| BPTI | [val] | [std] |
| SPINK1 | [val] | [std] |
| SPINK7 | [val] | [std] |

## Figures
1. Pairwise RMSD heatmap
2. MD RMSD time series
3. Summary bar chart

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
