# EXP-11: KLK5 Substrate Kinetics — Implementation Guide

**Experiment ID:** EXP-11  
**Feature ID:** F-11 (benchmarks.md)  
**Category:** Kinetic (Qualitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates KLK5 model by confirming substrate (Boc-VPR-AMC) adopts productive binding geometry in S1-S3 subsites with correct catalytic triad distances (Ser195 Oγ–scissile bond < 3.5 Å). Km = 0.20 mM, kcat = 196.8 min⁻¹ (Brattsand 2005, Michael 2005) as reference. Qualitative/structural validation — not direct kinetics computation.

### 2. Hypothesis

**H₁:** Substrate P1 Arg–S1 salt bridge maintained. **H₂:** Catalytic distances < 3.5 Å in 10 ns MD.

### 3. Protocol

Use KLK5 from EXP-01 (PDB 2PSX). Dock VPR tripeptide. 10 ns short MD. Monitor catalytic distances, P1 position, substrate RMSD. PASS: productive geometry maintained; MARGINAL: transient departures < 20%; FAIL: substrate dissociates.

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-11_klk5_substrate_kinetics"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import SystemConfig, MinimizationConfig, EquilibrationConfig, ProductionConfig
from src.prep.pdb_fetch import fetch_pdb
from src.prep.pdb_clean import clean_structure
from src.prep.protonate import assign_protonation
from src.prep.topology import build_topology
from src.prep.solvate import solvate_system
from src.simulate.minimizer import minimize_energy
from src.simulate.equilibrate import run_nvt, run_npt
from src.simulate.production import run_production
from src.simulate.platform import select_platform
from src.analyze.trajectory import load_trajectory
from src.analyze.structural import compute_rmsd
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Prepare KLK5 and Substrate

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

# Use KLK5 structure from EXP-01
klk5_pdb = fetch_pdb("2PSX", data_dir)
klk5_clean = clean_structure(klk5_pdb, chains_to_keep=["A"],
                              remove_heteroatoms=True, remove_waters=True, model_index=1)

# Use simplified VPR tripeptide (standard amino acids)
# Build substrate-KLK5 complex by manual placement or from EXP-01 docking
# P1 Arg → S1 (Asp189), P2 Pro → S2, P3 Val → S3
substrate_complex = data_dir / "klk5_vpr_complex.pdb"
if not substrate_complex.exists():
    print(f"Place docked KLK5-VPR complex at: {substrate_complex}")
    print("Dock P1 Arg into S1 pocket using known serine protease substrate orientation")
```

### Step 3: System Setup and Short MD

```python
from openmm.app import PME, Simulation, PDBFile
from openmm import LangevinMiddleIntegrator, XmlSerializer
import openmm

complex_protonated = assign_protonation(substrate_complex, ph=7.4, force_field="AMBER", use_propka=True)
sys_config = SystemConfig()
topology, system, modeller = build_topology(complex_protonated, sys_config, nonbonded_method=PME, nonbonded_cutoff_nm=1.0)
modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)

solvated_pdb = OUTPUT_DIR / "solvated_complex.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

platform = select_platform()
integrator = LangevinMiddleIntegrator(310*openmm.unit.kelvin, 1.0/openmm.unit.picosecond, 0.002*openmm.unit.picoseconds)
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

minimize_energy(sim, MinimizationConfig())
eq_output = OUTPUT_DIR / "equilibration"
eq_output.mkdir(exist_ok=True)
run_nvt(sim, EquilibrationConfig(), eq_output)
run_npt(sim, EquilibrationConfig(), eq_output)

# Short 10 ns production
short_prod_config = ProductionConfig()
short_prod_config.duration_ns = 10.0
prod_output = OUTPUT_DIR / "production"
prod_output.mkdir(exist_ok=True)
run_production(sim, short_prod_config, prod_output)
```

### Step 4: Catalytic Distance Analysis

```python
traj = load_trajectory(prod_output / "production_trajectory.dcd", solvated_pdb, stride=1)
topology_md = traj.topology

# Identify catalytic triad atoms (KLK5 numbering)
# Ser195 Oγ, His57 Nε2, Asp102 Oδ (renumbered for KLK5)
ser_og = topology_md.select("resname SER and name OG and chainid 0")
his_ne2 = topology_md.select("resname HIS and name NE2 and chainid 0")
asp_od = topology_md.select("resname ASP and name OD1 and chainid 0")

# P1 Arg scissile carbon
arg_p1_c = topology_md.select("resname ARG and name C and chainid 1")  # substrate chain

# Compute catalytic distances
if len(ser_og) > 0 and len(arg_p1_c) > 0:
    d_ser_scissile = md.compute_distances(traj, [[ser_og[0], arg_p1_c[0]]])[:, 0] * 10  # nm→Å
    print(f"Ser-Oγ to scissile C: mean={np.mean(d_ser_scissile):.2f} Å, max={np.max(d_ser_scissile):.2f} Å")
    productive_fraction = np.mean(d_ser_scissile < 3.5) * 100
    print(f"Productive geometry (< 3.5 Å): {productive_fraction:.1f}% of frames")

# P1 Arg guanidinium - S1 Asp salt bridge
arg_cz = topology_md.select("resname ARG and name CZ and chainid 1")
asp189_cg = topology_md.select("resname ASP and name CG and chainid 0")
if len(arg_cz) > 0 and len(asp189_cg) > 0:
    d_salt_bridge = md.compute_distances(traj, [[arg_cz[0], asp189_cg[-1]]])[:, 0] * 10
    salt_bridge_fraction = np.mean(d_salt_bridge < 4.0) * 100
    print(f"P1 Arg-S1 Asp salt bridge: mean={np.mean(d_salt_bridge):.2f} Å, fraction={salt_bridge_fraction:.1f}%")
```

### Step 5: Classification

```python
if productive_fraction > 80 and salt_bridge_fraction > 80:
    classification = "PASS"
elif productive_fraction > 60:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-11", "feature_id": "F-11",
    "productive_fraction_pct": float(productive_fraction),
    "salt_bridge_fraction_pct": float(salt_bridge_fraction),
    "mean_ser_scissile_dist_ang": float(np.mean(d_ser_scissile)),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-11 Classification: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Catalytic Distance Time Series

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
time_ns = np.arange(len(d_ser_scissile)) * 0.01
axes[0].plot(time_ns, d_ser_scissile, 'b-', linewidth=0.5)
axes[0].axhline(y=3.5, color='red', linestyle='--', label='3.5 Å threshold')
axes[0].set_ylabel('Ser-Oγ → scissile C (Å)', fontsize=12)
axes[0].legend()
axes[1].plot(time_ns, d_salt_bridge, 'g-', linewidth=0.5)
axes[1].axhline(y=4.0, color='red', linestyle='--', label='4.0 Å threshold')
axes[1].set_ylabel('P1 Arg → S1 Asp (Å)', fontsize=12)
axes[1].set_xlabel('Time (ns)', fontsize=12)
axes[1].legend()
plt.suptitle(f'EXP-11: KLK5 Substrate Catalytic Distances — {classification}', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-11_catalytic_distances.png", dpi=300)
plt.close(fig)
```

### Figure 2: Substrate RMSD

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
substrate_atoms = topology_md.select("chainid 1")
rmsd_sub = compute_rmsd(traj.atom_slice(substrate_atoms), traj.atom_slice(substrate_atoms)[0], atom_sel="all")
ax.plot(time_ns, rmsd_sub * 10, 'b-', linewidth=0.8)
ax.axhline(y=2.0, color='red', linestyle='--', label='2.0 Å threshold')
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Substrate RMSD (Å)', fontsize=12)
ax.set_title('EXP-11: Substrate Heavy Atom RMSD', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-11_substrate_rmsd.png", dpi=300)
plt.close(fig)
```

### Figure 3: Active Site Geometry Summary

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
metrics = ['Ser-Oγ→scissile\n(< 3.5 Å)', 'P1 Arg→S1 Asp\n(< 4.0 Å)']
fractions = [productive_fraction, salt_bridge_fraction]
colors = ['steelblue' if f > 80 else 'orange' if f > 60 else 'red' for f in fractions]
ax.bar(metrics, fractions, color=colors, edgecolor='black')
ax.axhline(y=80, color='green', linestyle='--', label='PASS threshold')
ax.set_ylabel('% frames in criterion', fontsize=12)
ax.set_title(f'EXP-11: Active Site Geometry — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-11_geometry_summary.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-11: KLK5 Substrate Kinetics — Results Report

**Experiment ID:** EXP-11  **Feature ID:** F-11  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## 1–11. [Standard sections per §28]
---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
