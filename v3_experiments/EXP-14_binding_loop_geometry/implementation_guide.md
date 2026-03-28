# EXP-14: Canonical Binding Loop Geometry — Implementation Guide

**Experiment ID:** EXP-14  
**Feature ID:** F-14 (benchmarks.md)  
**Category:** Structural/Quantitative  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates that the canonical P3–P3' binding loop of Kazal-type serine protease inhibitors maintains proper extended β-sheet geometry in MD simulations. Tests across three systems: BPTI–trypsin (EXP-04), SPINK7–KLK5 (EXP-01), SPINK1–trypsin (EXP-06). Metrics: loop RMSD < 1.5 Å, P1 φ/ψ in extended β-sheet (φ ≈ -120°, ψ ≈ +130°, > 90% frames), oxyanion hole H-bond occupancy > 85%, disulfide S–S distance 1.8–2.2 Å.

### 2. Benchmark Values

| Metric | Criterion | Source |
|--------|-----------|--------|
| P3–P3' loop RMSD | < 1.5 Å | Laskowski & Qasim 2000 |
| P1 φ angle | -120° ± 30° | β-sheet Ramachandran |
| P1 ψ angle | +130° ± 30° | β-sheet Ramachandran |
| β-sheet fraction | > 90% of frames | Standard |
| Oxyanion H-bond | > 85% occupancy | Radisky & Bhatt 2002 |
| Disulfide S–S | 1.8–2.2 Å | Standard |

### 3. Classification Criteria (§25.1)

- PASS: RMSD < 1.5 Å AND β-sheet > 90% (for ≥ 2/3 systems)
- MARGINAL: RMSD 1.5–2.5 Å AND β-sheet 70–90%
- FAIL: RMSD > 2.5 Å OR β-sheet < 70%

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-14_binding_loop_geometry"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory, align_trajectory
from src.analyze.structural import compute_rmsd, compute_rmsf, compute_hbonds
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Define Systems and Load Trajectories

```python
EXP_ROOT = PROJECT_ROOT / "v3_experiments"

SYSTEMS = {
    "BPTI-trypsin": {
        "exp": "EXP-04_bpti_trypsin",
        "traj_pattern": "outputs/production/production_trajectory.dcd",
        "top": "outputs/solvated_complex.pdb",
        "inhibitor_chain": 1,  # chainid for BPTI
        "p1_resname": "LYS",
        "p1_resid_range": (13, 17),  # P3(C14)–P3'(A18) in BPTI (0-indexed approx)
        "disulfides": [(4, 54), (13, 37), (29, 50)],  # BPTI numbering (0-indexed)
    },
    "SPINK7-KLK5": {
        "exp": "EXP-01_spink7_klk5_binding",
        "traj_pattern": "outputs/production/production_trajectory.dcd",
        "top": "outputs/solvated_complex.pdb",
        "inhibitor_chain": 1,
        "p1_resname": "ARG",
        "p1_resid_range": None,  # Determine from structure
        "disulfides": [],
    },
    "SPINK1-trypsin": {
        "exp": "EXP-06_spink1_trypsin",
        "traj_pattern": "outputs/production/production_trajectory.dcd",
        "top": "outputs/solvated_complex.pdb",
        "inhibitor_chain": 1,
        "p1_resname": "LYS",
        "p1_resid_range": None,
        "disulfides": [],
    },
}

trajectories = {}
for name, info in SYSTEMS.items():
    traj_path = EXP_ROOT / info["exp"] / info["traj_pattern"]
    top_path = EXP_ROOT / info["exp"] / info["top"]
    if traj_path.exists() and top_path.exists():
        traj = load_trajectory(str(traj_path), str(top_path), stride=10)
        trajectories[name] = traj
        print(f"{name}: {traj.n_frames} frames loaded")
    else:
        print(f"WARNING: {name} trajectory not found at {traj_path}")
```

### Step 3: Binding Loop RMSD Analysis

```python
loop_rmsd_results = {}

for name, traj in trajectories.items():
    info = SYSTEMS[name]
    topo = traj.topology

    # Select inhibitor chain binding loop (P3-P3', ~7 residues around P1)
    inh_residues = [r for r in topo.residues if r.chain.index == info["inhibitor_chain"]]

    # Find P1 residue (the key contact residue)
    p1_candidates = [r for r in inh_residues if r.name == info["p1_resname"]]
    if p1_candidates:
        p1_res = p1_candidates[0]
        p1_idx = p1_res.index
        # P3-P3' = P1 ± 3 residues
        loop_resids = list(range(max(0, p1_idx - 3), p1_idx + 4))
        loop_atoms = topo.select(f"resid {' '.join(str(r) for r in loop_resids)} and name CA C N O")

        if len(loop_atoms) > 0:
            loop_traj = traj.atom_slice(loop_atoms)
            loop_traj_aligned = align_trajectory(loop_traj, loop_traj[0])
            rmsd_vals = md.rmsd(loop_traj_aligned, loop_traj_aligned, 0) * 10  # nm → Å
            loop_rmsd_results[name] = rmsd_vals
            print(f"{name}: loop RMSD = {np.mean(rmsd_vals):.2f} ± {np.std(rmsd_vals):.2f} Å")
```

### Step 4: Ramachandran Analysis for P1 Residue

```python
phi_psi_results = {}

for name, traj in trajectories.items():
    info = SYSTEMS[name]
    topo = traj.topology

    inh_residues = [r for r in topo.residues if r.chain.index == info["inhibitor_chain"]]
    p1_candidates = [r for r in inh_residues if r.name == info["p1_resname"]]

    if p1_candidates:
        p1_res = p1_candidates[0]
        p1_idx = p1_res.index

        # Compute phi/psi for P1
        phi_indices = md.compute_phi(traj)[0]
        psi_indices = md.compute_psi(traj)[0]

        # Find which phi/psi indices correspond to P1
        phi_angles, psi_angles = None, None
        phi_all = md.compute_phi(traj)
        psi_all = md.compute_psi(traj)

        for i, idx_set in enumerate(phi_all[0]):
            # phi angle defined by C(i-1)-N(i)-CA(i)-C(i)
            atom = topo.atom(idx_set[2])  # CA
            if atom.residue.index == p1_idx:
                phi_angles = np.degrees(phi_all[1][:, i])
                break

        for i, idx_set in enumerate(psi_all[0]):
            atom = topo.atom(idx_set[1])  # N
            if atom.residue.index == p1_idx:
                psi_angles = np.degrees(psi_all[1][:, i])
                break

        if phi_angles is not None and psi_angles is not None:
            # β-sheet criteria: φ ≈ -120° ± 30°, ψ ≈ +130° ± 30°
            in_beta = ((-150 < phi_angles) & (phi_angles < -90) &
                       (100 < psi_angles) & (psi_angles < 160))
            beta_fraction = np.mean(in_beta) * 100
            phi_psi_results[name] = {
                "phi": phi_angles, "psi": psi_angles,
                "beta_fraction": beta_fraction,
                "phi_mean": float(np.mean(phi_angles)),
                "psi_mean": float(np.mean(psi_angles)),
            }
            print(f"{name}: P1 β-sheet fraction = {beta_fraction:.1f}%, "
                  f"φ = {np.mean(phi_angles):.1f}°, ψ = {np.mean(psi_angles):.1f}°")
```

### Step 5: Oxyanion Hole H-Bond Occupancy

```python
hbond_results = {}

for name, traj in trajectories.items():
    # Compute all H-bonds
    hbonds = compute_hbonds(traj, freq=0.0)  # all H-bonds

    # Filter for oxyanion hole interactions
    # Oxyanion hole: backbone NH of protease residues flanking active site Ser
    # (typically Gly193 NH, Ser195 NH in chymotrypsin numbering)
    # → P1 carbonyl oxygen
    topo = traj.topology
    info = SYSTEMS[name]

    inh_residues = [r for r in topo.residues if r.chain.index == info["inhibitor_chain"]]
    p1_candidates = [r for r in inh_residues if r.name == info["p1_resname"]]

    if p1_candidates:
        p1_res = p1_candidates[0]
        p1_O = [a.index for a in p1_res.atoms if a.name == "O"]

        if p1_O:
            # Count frames where P1 O is H-bonded to protease backbone N
            protease_N = topo.select(f"chainid 0 and name N")
            oxyanion_occ = 0
            for frame_idx in range(traj.n_frames):
                for p_n in protease_N:
                    dist = md.compute_distances(traj[frame_idx], [[p1_O[0], p_n]])[0, 0]
                    if dist < 0.35:  # 3.5 Å in nm
                        oxyanion_occ += 1
                        break
            occupancy = oxyanion_occ / traj.n_frames * 100
            hbond_results[name] = occupancy
            print(f"{name}: oxyanion H-bond occupancy = {occupancy:.1f}%")
```

### Step 6: Disulfide Bond Integrity

```python
disulfide_results = {}

for name, traj in trajectories.items():
    info = SYSTEMS[name]
    if not info["disulfides"]:
        disulfide_results[name] = "N/A"
        continue

    topo = traj.topology
    inh_residues = [r for r in topo.residues if r.chain.index == info["inhibitor_chain"]]

    for res_i, res_j in info["disulfides"]:
        sg_i = topo.select(f"resid {res_i} and name SG and chainid {info['inhibitor_chain']}")
        sg_j = topo.select(f"resid {res_j} and name SG and chainid {info['inhibitor_chain']}")
        if len(sg_i) > 0 and len(sg_j) > 0:
            ss_dist = md.compute_distances(traj, [[sg_i[0], sg_j[0]]])[:, 0] * 10  # Å
            intact = np.mean((ss_dist > 1.8) & (ss_dist < 2.2)) * 100
            print(f"{name}: SS {res_i}-{res_j}: mean={np.mean(ss_dist):.2f} Å, intact={intact:.1f}%")
            disulfide_results.setdefault(name, {})[f"{res_i}-{res_j}"] = {
                "mean_dist": float(np.mean(ss_dist)), "intact_pct": float(intact)
            }
```

### Step 7: Classification

```python
system_classifications = {}

for name in SYSTEMS.keys():
    rmsd_pass = name in loop_rmsd_results and np.mean(loop_rmsd_results[name]) < 1.5
    beta_pass = name in phi_psi_results and phi_psi_results[name]["beta_fraction"] > 90
    rmsd_marginal = name in loop_rmsd_results and np.mean(loop_rmsd_results[name]) < 2.5
    beta_marginal = name in phi_psi_results and phi_psi_results[name]["beta_fraction"] > 70

    if rmsd_pass and beta_pass:
        system_classifications[name] = "PASS"
    elif rmsd_marginal and beta_marginal:
        system_classifications[name] = "MARGINAL"
    else:
        system_classifications[name] = "FAIL"
    print(f"{name}: {system_classifications[name]}")

n_pass = sum(1 for c in system_classifications.values() if c == "PASS")
overall = "PASS" if n_pass >= 2 else "MARGINAL" if n_pass >= 1 else "FAIL"

results = {
    "experiment_id": "EXP-14", "feature_id": "F-14",
    "loop_rmsd_means": {n: float(np.mean(v)) for n, v in loop_rmsd_results.items()},
    "beta_fractions": {n: v["beta_fraction"] for n, v in phi_psi_results.items()},
    "oxyanion_occupancy": hbond_results,
    "system_classifications": system_classifications,
    "overall_classification": overall,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-14 overall: {overall}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Loop RMSD Time Series (All Systems)

```python
fig, axes = plt.subplots(len(loop_rmsd_results), 1, figsize=(12, 4*len(loop_rmsd_results)), sharex=True)
if len(loop_rmsd_results) == 1:
    axes = [axes]

for ax, (name, rmsd_vals) in zip(axes, loop_rmsd_results.items()):
    time_ns = np.arange(len(rmsd_vals)) * 0.1  # stride=10
    ax.plot(time_ns, rmsd_vals, 'b-', linewidth=0.5)
    ax.axhline(y=1.5, color='green', linestyle='--', label='PASS (1.5 Å)')
    ax.axhline(y=2.5, color='orange', linestyle=':', label='MARGINAL (2.5 Å)')
    cls = system_classifications.get(name, "?")
    ax.set_ylabel(f'{name}\nRMSD (Å)', fontsize=11)
    ax.set_title(f'{name} — {cls}', fontsize=12)
    ax.legend(fontsize=9)

axes[-1].set_xlabel('Time (ns)', fontsize=12)
plt.suptitle('EXP-14: Binding Loop RMSD', fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-14_loop_rmsd.png", dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Figure 2: P1 Ramachandran Scatter Plot

```python
fig, axes = plt.subplots(1, len(phi_psi_results), figsize=(6*len(phi_psi_results), 5))
if len(phi_psi_results) == 1:
    axes = [axes]

for ax, (name, data) in zip(axes, phi_psi_results.items()):
    phi = data["phi"]
    psi = data["psi"]
    ax.scatter(phi, psi, s=1, alpha=0.3, c='blue')
    # Mark β-sheet region
    from matplotlib.patches import Rectangle
    rect = Rectangle((-150, 100), 60, 60, linewidth=2, edgecolor='green',
                     facecolor='green', alpha=0.15, label=f'β-sheet ({data["beta_fraction"]:.0f}%)')
    ax.add_patch(rect)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel('φ (°)', fontsize=12)
    ax.set_ylabel('ψ (°)', fontsize=12)
    ax.set_title(f'{name}: P1 Ramachandran', fontsize=12)
    ax.legend()
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5)

plt.suptitle('EXP-14: P1 Residue Ramachandran Analysis', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-14_ramachandran.png", dpi=300)
plt.close(fig)
```

### Figure 3: Summary Dashboard

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
systems_list = list(SYSTEMS.keys())

# Panel A: Mean loop RMSD
rmsd_means = [np.mean(loop_rmsd_results.get(s, [0])) for s in systems_list]
colors_rmsd = ['green' if v < 1.5 else 'orange' if v < 2.5 else 'red' for v in rmsd_means]
axes[0].bar(systems_list, rmsd_means, color=colors_rmsd, edgecolor='black')
axes[0].axhline(1.5, color='green', linestyle='--')
axes[0].set_ylabel('Mean RMSD (Å)', fontsize=12)
axes[0].set_title('Loop RMSD', fontsize=13)
axes[0].tick_params(axis='x', rotation=20)

# Panel B: Beta-sheet fraction
beta_vals = [phi_psi_results.get(s, {}).get("beta_fraction", 0) for s in systems_list]
colors_beta = ['green' if v > 90 else 'orange' if v > 70 else 'red' for v in beta_vals]
axes[1].bar(systems_list, beta_vals, color=colors_beta, edgecolor='black')
axes[1].axhline(90, color='green', linestyle='--')
axes[1].set_ylabel('β-sheet fraction (%)', fontsize=12)
axes[1].set_title('P1 Geometry', fontsize=13)
axes[1].tick_params(axis='x', rotation=20)

# Panel C: Classification
cls_colors = {'PASS': 'green', 'MARGINAL': 'orange', 'FAIL': 'red'}
cls_vals = [system_classifications.get(s, 'FAIL') for s in systems_list]
axes[2].bar(systems_list, [1]*len(systems_list),
            color=[cls_colors[c] for c in cls_vals], edgecolor='black')
for i, (s, c) in enumerate(zip(systems_list, cls_vals)):
    axes[2].text(i, 0.5, c, ha='center', va='center', fontsize=14, fontweight='bold')
axes[2].set_ylabel('')
axes[2].set_title(f'Classification (Overall: {overall})', fontsize=13)
axes[2].tick_params(axis='x', rotation=20)

plt.suptitle('EXP-14: Binding Loop Geometry Summary', fontsize=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-14_summary_dashboard.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-14: Canonical Binding Loop Geometry — Results Report

**Experiment ID:** EXP-14  **Feature ID:** F-14  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results

| System | Loop RMSD (Å) | β-sheet (%) | Oxyanion Occ. (%) | Classification |
|--------|---------------|-------------|-------------------|----------------|
| BPTI-trypsin | [value] | [value] | [value] | [PASS/MARGINAL/FAIL] |
| SPINK7-KLK5 | [value] | [value] | [value] | [PASS/MARGINAL/FAIL] |
| SPINK1-trypsin | [value] | [value] | [value] | [PASS/MARGINAL/FAIL] |

Overall: [X]/3 PASS → [classification]

## Figures
1. Loop RMSD time series
2. P1 Ramachandran scatter
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
