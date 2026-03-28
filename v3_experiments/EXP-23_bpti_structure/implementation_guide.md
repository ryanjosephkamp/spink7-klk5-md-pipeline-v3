# EXP-23: BPTI Crystal Structure Reproduction — Implementation Guide

**Experiment ID:** EXP-23  
**Feature ID:** F-23 (benchmarks.md)  
**Category:** Structural (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates that the pipeline's MD simulations of free BPTI (PDB 4PTI, 1.0 Å resolution) faithfully reproduce the crystal structure geometry. Backbone RMSD < 1.5 Å over 100 ns, disulfide bonds maintained (2.0 ± 0.15 Å), secondary structure >90% occupancy.

### 2. Hypotheses

- H₁: Time-averaged backbone RMSD < 1.5 Å vs crystal (PASS criterion)
- H₂: All three disulfide bonds (C5-C55, C14-C38, C30-C51) maintain S-S distance 2.0 ± 0.15 Å
- H₃: Secondary structure (α-helix, β-sheet) persists >90% via DSSP

### 3. Classification (§25.1)

- PASS: Backbone RMSD < 1.5 Å (mean)
- MARGINAL: 1.5–2.5 Å
- FAIL: > 2.5 Å

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup and Structure Fetch

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-23_bpti_structure"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.prep.pdb_fetch import fetch_pdb
from src.prep.structure_prep import clean_structure, assign_protonation
from src.prep.topology import build_topology, solvate_system
from src.simulate.equilibration import minimize_energy, run_nvt, run_npt
from src.simulate.production import run_production
from src.analyze.trajectory import load_trajectory, align_trajectory
from src.analyze.structural import compute_rmsd, compute_rmsf
from src.config import (SystemConfig, MinimizationConfig, EquilibrationConfig,
                         ProductionConfig)
from src.simulate.platform_util import select_platform
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Structure Preparation

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)
pdb_path = fetch_pdb("4PTI", data_dir)
clean_path = clean_structure(pdb_path, chains_to_keep=None,
                              remove_heteroatoms=True, remove_waters=True)
prot_path = assign_protonation(clean_path, ph=7.4, force_field="AMBER", use_propka=True)
print(f"Protonated: {prot_path}")
```

### Step 3: Build System and Solvate

```python
sys_config = SystemConfig(
    force_field="amber14-all.xml",
    water_model="amber14/tip3p.xml",
    box_padding_nm=1.2,
    ionic_strength_M=0.15,
    temperature_K=310.0,
)
topology, system, modeller = build_topology(prot_path, sys_config,
                                             nonbonded_method="PME",
                                             nonbonded_cutoff_nm=1.0)
modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)
print(f"Solvated: {n_waters} waters, {n_pos} Na+, {n_neg} Cl-")

# Save solvated PDB
from openmm.app import PDBFile
solvated_pdb = OUTPUT_DIR / "solvated_bpti.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
```

### Step 4: Equilibration

```python
import openmm as mm
from openmm import unit

platform = select_platform()
integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picosecond,
                                          0.002*unit.picoseconds)
simulation = mm.app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# Minimize
min_config = MinimizationConfig(max_iterations=10000, tolerance_kj_per_mol_nm=10.0)
min_result = minimize_energy(simulation, min_config)
print(f"Minimized: {min_result}")

# NVT
equil_dir = OUTPUT_DIR / "equilibration"
equil_dir.mkdir(exist_ok=True)
nvt_config = EquilibrationConfig(ensemble="NVT", duration_ps=500, temperature_K=310,
                                   restraint_strength_kj=1000.0)
nvt_result = run_nvt(simulation, nvt_config, equil_dir)
print(f"NVT done: {nvt_result}")

# NPT
npt_config = EquilibrationConfig(ensemble="NPT", duration_ps=1000, temperature_K=310,
                                   pressure_atm=1.0, barostat_interval=25)
npt_result = run_npt(simulation, npt_config, equil_dir)
print(f"NPT done: {npt_result}")
```

### Step 5: Production MD (100 ns)

```python
prod_dir = OUTPUT_DIR / "production"
prod_dir.mkdir(exist_ok=True)
prod_config = ProductionConfig(duration_ns=100, temperature_K=310,
                                save_interval_ps=10, checkpoint_interval_ps=100)
prod_result = run_production(simulation, prod_config, prod_dir)
print(f"Production done: {prod_result}")
```

### Step 6: Structural Analysis

```python
# Load trajectory
traj = load_trajectory(str(prod_dir / "production_trajectory.dcd"),
                        str(solvated_pdb), stride=1)
topo = traj.topology
crystal = md.load(str(prot_path))
print(f"Loaded: {traj.n_frames} frames, {crystal.n_atoms} crystal atoms")

# Backbone RMSD vs crystal
protein_atoms = topo.select("protein")
traj_prot = traj.atom_slice(protein_atoms)
traj_aligned = traj_prot.superpose(crystal, atom_indices=crystal.topology.select("backbone"))
bb_idx = crystal.topology.select("backbone")
rmsd_bb = md.rmsd(traj_aligned, crystal, atom_indices=bb_idx) * 10  # nm → Å

time_ns = np.arange(len(rmsd_bb)) * 0.01
print(f"RMSD: mean={np.mean(rmsd_bb):.2f} ± {np.std(rmsd_bb):.2f} Å")

# Per-residue RMSF
ca_idx = crystal.topology.select("name CA")
rmsf_ca = md.rmsf(traj_aligned, crystal, atom_indices=ca_idx) * 10
print(f"RMSF: mean={np.mean(rmsf_ca):.2f} Å, max={np.max(rmsf_ca):.2f} Å")

# B-factor comparison
bfactors_xtal = np.array([crystal.topology.atom(i).element for i in ca_idx])  # placeholder
# RMSF_xtal = sqrt(3*B/(8*pi^2))
# Will use crystal B-factors from PDB file

# Disulfide distances
disulfides = [("CYS5", "CYS55"), ("CYS14", "CYS38"), ("CYS30", "CYS51")]
ss_results = {}
for res_a, res_b in disulfides:
    idx_a = topo.select(f"resname CYS and resid {int(res_a[3:])} and name SG")
    idx_b = topo.select(f"resname CYS and resid {int(res_b[3:])} and name SG")
    if len(idx_a) > 0 and len(idx_b) > 0:
        d_ss = md.compute_distances(traj_prot, [[idx_a[0], idx_b[0]]])[:, 0] * 10
        ss_results[f"{res_a}-{res_b}"] = {"mean": float(np.mean(d_ss)),
                                            "std": float(np.std(d_ss))}
        print(f"  {res_a}-{res_b}: {np.mean(d_ss):.3f} ± {np.std(d_ss):.3f} Å")

# Secondary structure via DSSP
dssp = md.compute_dssp(traj_prot)
helix_residues = list(range(46, 56))  # residues 47-56
sheet_residues = list(range(17, 35))   # residues 18-24, 29-35
helix_occ = np.mean(dssp[:, helix_residues] == 'H') * 100
sheet_occ = np.mean(np.isin(dssp[:, sheet_residues], ['E', 'B'])) * 100
print(f"Helix occupancy: {helix_occ:.1f}%, Sheet occupancy: {sheet_occ:.1f}%")
```

### Step 7: Classification

```python
rmsd_mean = np.mean(rmsd_bb)
if rmsd_mean < 1.5:
    classification = "PASS"
elif rmsd_mean < 2.5:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-23", "feature_id": "F-23",
    "rmsd_mean_ang": float(rmsd_mean),
    "rmsd_std_ang": float(np.std(rmsd_bb)),
    "rmsd_max_ang": float(np.max(rmsd_bb)),
    "rmsf_mean_ang": float(np.mean(rmsf_ca)),
    "disulfide_distances": ss_results,
    "helix_occupancy_pct": float(helix_occ),
    "sheet_occupancy_pct": float(sheet_occ),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-23: RMSD={rmsd_mean:.2f} Å → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Backbone RMSD Time Series

```python
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(time_ns, rmsd_bb, 'b-', linewidth=0.5, alpha=0.7)
ax.axhline(y=1.5, color='green', linestyle='--', label='PASS (1.5 Å)')
ax.axhline(y=2.5, color='orange', linestyle='--', label='MARGINAL (2.5 Å)')
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Backbone RMSD (Å)', fontsize=12)
ax.set_title(f'EXP-23: BPTI RMSD vs Crystal (4PTI) — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-23_rmsd_timeseries.png", dpi=300)
plt.close(fig)
```

### Figure 2: Per-Residue RMSF Profile

```python
fig, ax = plt.subplots(figsize=(14, 5))
residues = np.arange(1, len(rmsf_ca) + 1)
ax.bar(residues, rmsf_ca, color='steelblue', edgecolor='none', width=0.8)
# Highlight disulfide-bonded residues
for r in [5, 14, 30, 38, 51, 55]:
    ax.axvline(x=r, color='gold', linestyle=':', alpha=0.7)
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Core/loop threshold')
ax.set_xlabel('Residue Number', fontsize=12)
ax.set_ylabel('Cα RMSF (Å)', fontsize=12)
ax.set_title('EXP-23: Per-Residue RMSF', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-23_rmsf_profile.png", dpi=300)
plt.close(fig)
```

### Figure 3: Disulfide Distances + Secondary Structure Summary

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Disulfides
labels = list(ss_results.keys())
means = [ss_results[k]["mean"] for k in labels]
stds = [ss_results[k]["std"] for k in labels]
ax1.bar(labels, means, yerr=stds, color=['#2ecc71', '#3498db', '#e74c3c'],
         edgecolor='black', capsize=5)
ax1.axhspan(1.85, 2.15, alpha=0.1, color='green')
ax1.set_ylabel('S-S Distance (Å)', fontsize=12)
ax1.set_title('Disulfide Bond Distances', fontsize=13)

# Secondary structure
ss_labels = ['α-helix', 'β-sheet']
ss_vals = [helix_occ, sheet_occ]
colors_ss = ['#e74c3c' if v < 90 else '#2ecc71' for v in ss_vals]
ax2.bar(ss_labels, ss_vals, color=colors_ss, edgecolor='black')
ax2.axhline(y=90, color='green', linestyle='--', label='90% threshold')
ax2.set_ylabel('Occupancy (%)', fontsize=12)
ax2.set_title('Secondary Structure Persistence', fontsize=13)
ax2.set_ylim(0, 105)
ax2.legend()

plt.suptitle(f'EXP-23: BPTI Structural Integrity — {classification}', fontsize=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-23_structural_summary.png", dpi=300)
plt.close(fig)
```

### Figure 4: RMSD Distribution

```python
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(rmsd_bb, bins=50, color='steelblue', edgecolor='black', density=True)
ax.axvline(x=1.5, color='green', linestyle='--', linewidth=2, label='PASS threshold')
ax.set_xlabel('Backbone RMSD (Å)', fontsize=14)
ax.set_ylabel('Density', fontsize=14)
ax.set_title(f'EXP-23: RMSD Distribution (mean={rmsd_mean:.2f} Å)', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-23_rmsd_histogram.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-23: BPTI Crystal Structure Reproduction — Results Report

**Experiment ID:** EXP-23  **Feature ID:** F-23  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| Backbone RMSD (mean) | [val] Å | < 1.5 Å | [P/M/F] |
| Backbone RMSD (max) | [val] Å | — | — |
| RMSF (mean) | [val] Å | — | — |
| Disulfide C5-C55 | [val] Å | 1.85-2.15 | [P/F] |
| Disulfide C14-C38 | [val] Å | 1.85-2.15 | [P/F] |
| Disulfide C30-C51 | [val] Å | 1.85-2.15 | [P/F] |
| α-Helix occupancy | [val]% | >90% | [P/F] |
| β-Sheet occupancy | [val]% | >90% | [P/F] |

## Figures
1. RMSD time series
2. RMSF per-residue profile
3. Structural summary (disulfides + SS)
4. RMSD distribution

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
