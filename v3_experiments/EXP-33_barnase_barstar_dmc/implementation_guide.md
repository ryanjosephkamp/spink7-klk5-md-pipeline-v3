# EXP-33: Barnase–Barstar Double-Mutant Cycle Analysis — Implementation Guide

**Experiment ID:** EXP-33  
**Feature ID:** F-33 (benchmarks.md)  
**Category:** Mutational (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates FEP on the barnase–barstar double-mutant cycle (DMC) dataset (~45 pairs, Schreiber & Fersht 1995). Computes single-mutation ΔΔG and pairwise coupling energies (ΔΔΔG). The ultimate test of computational mutagenesis accuracy. Staged approach: Stage 1 (10 singles), Stage 2 (5 key doubles), Stage 3 (remaining).

### 2. Hypotheses

- H₁: Single-mutation ΔΔG correlation R > 0.7
- H₂: Correctly identify top 5 coupled pairs (|ΔΔΔG| > 1.5 kcal/mol)
- H₃: RMSE < 2.0 kcal/mol

### 3. Classification (§25.1)

- PASS: R > 0.7, RMSE < 2.0, top 5 coupled pairs correct
- MARGINAL: R > 0.5, RMSE < 3.0, top 3 correct
- FAIL: R < 0.5, RMSE > 3.0, or coupled pairs wrong

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path
from scipy import stats
from itertools import combinations

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-33_barnase_barstar_dmc"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.prep.pdb_fetch import fetch_pdb
from src.prep.structure_prep import clean_structure, assign_protonation
from src.prep.topology import build_topology, solvate_system
from src.simulate.equilibration import minimize_energy, run_nvt, run_npt
from src.simulate.fep import run_fep_campaign
from src.analyze.mbar import solve_mbar, bootstrap_mbar_uncertainty
from src.config import (SystemConfig, MinimizationConfig, EquilibrationConfig,
                         FEPConfig, MBARConfig)
from src.simulate.platform_util import select_platform
from src.physics.units import kj_to_kcal
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Prepare Barnase–Barstar Complex

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

pdb_path = fetch_pdb("1BRS", data_dir)
# 1BRS has 4 chains: A(barnase), D(barstar) is one complex copy
clean_path = clean_structure(pdb_path, chains_to_keep=["A", "D"],
                              remove_heteroatoms=True, remove_waters=True)
prot_path = assign_protonation(clean_path, ph=7.4, force_field="AMBER", use_propka=True)

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

from openmm.app import PDBFile
solvated_pdb = OUTPUT_DIR / "solvated_barnase_barstar.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"System: {n_waters} waters, {n_pos} Na+, {n_neg} Cl-")

# Equilibrate
import openmm as mm
from openmm import unit

platform = select_platform()
integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picosecond,
                                          0.002*unit.picoseconds)
simulation = mm.app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

min_config = MinimizationConfig(max_iterations=10000, tolerance_kj_per_mol_nm=10.0)
minimize_energy(simulation, min_config)

equil_dir = OUTPUT_DIR / "equilibration"
equil_dir.mkdir(exist_ok=True)
nvt_config = EquilibrationConfig(ensemble="NVT", duration_ps=500, temperature_K=310,
                                   restraint_strength_kj=1000.0)
run_nvt(simulation, nvt_config, equil_dir)
npt_config = EquilibrationConfig(ensemble="NPT", duration_ps=1000, temperature_K=310,
                                   pressure_atm=1.0, barostat_interval=25)
run_npt(simulation, npt_config, equil_dir)
print("Equilibration complete.")
```

### Step 3: Define Mutation Panel

```python
# Key barnase mutations
barnase_mutations = [
    {"name": "K27A", "chain": "barnase", "ddg_exp": 4.2},
    {"name": "R59A", "chain": "barnase", "ddg_exp": 4.8},
    {"name": "H102A", "chain": "barnase", "ddg_exp": 3.5},
    {"name": "R83A", "chain": "barnase", "ddg_exp": 3.8},
    {"name": "R87A", "chain": "barnase", "ddg_exp": 2.1},
]
# Key barstar mutations
barstar_mutations = [
    {"name": "D35A", "chain": "barstar", "ddg_exp": 5.5},
    {"name": "D39A", "chain": "barstar", "ddg_exp": 6.0},
    {"name": "Y29A", "chain": "barstar", "ddg_exp": 2.0},
    {"name": "E76A", "chain": "barstar", "ddg_exp": 1.5},
    {"name": "T42A", "chain": "barstar", "ddg_exp": 0.8},
]
all_single_mutations = barnase_mutations + barstar_mutations

# Key double-mutant pairs (most coupled from experiment)
double_mutations = [
    {"pair": ("R59A", "D35A"), "dddg_exp": 4.5},
    {"pair": ("K27A", "D39A"), "dddg_exp": 3.5},
    {"pair": ("H102A", "D35A"), "dddg_exp": 2.8},
    {"pair": ("R83A", "D35A"), "dddg_exp": 2.0},
    {"pair": ("R87A", "E76A"), "dddg_exp": 0.3},  # non-coupled control
]
print(f"Stage 1: {len(all_single_mutations)} single mutations")
print(f"Stage 2: {len(double_mutations)} double mutations")
```

### Step 4: Stage 1 — Single-Mutation FEP

```python
fep_config = FEPConfig(
    n_lambda=20,
    per_window_ns=2.0,
    temperature_K=310.0,
    soft_core_alpha=0.5,
    soft_core_power=1,
    annihilate_electrostatics=True,
    annihilate_sterics=False,
    n_equil_steps=5000,
    save_interval_ps=1.0,
)
mbar_config = MBARConfig(solver="robust", tolerance=1e-7, max_iterations=10000)

single_results = {}

for mut in all_single_mutations:
    mut_name = mut["name"]
    chain = mut["chain"]
    mut_dir = OUTPUT_DIR / "singles" / mut_name
    mut_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n--- {mut_name} ({chain}) ---")
    
    # Leg 1: Complex
    complex_dir = mut_dir / "complex"
    complex_dir.mkdir(exist_ok=True)
    fep_complex = run_fep_campaign(
        topology_path=str(solvated_pdb),
        mutation=mut_name,
        state="complex",
        fep_config=fep_config,
        output_dir=complex_dir,
    )
    dg_complex = solve_mbar(fep_complex["energy_files"], fep_complex["metadata"], mbar_config)
    
    # Leg 2: Free protein
    free_pdb = OUTPUT_DIR / f"solvated_{chain}_free.pdb"
    free_dir = mut_dir / "free"
    free_dir.mkdir(exist_ok=True)
    fep_free = run_fep_campaign(
        topology_path=str(free_pdb),
        mutation=mut_name,
        state="free",
        fep_config=fep_config,
        output_dir=free_dir,
    )
    dg_free = solve_mbar(fep_free["energy_files"], fep_free["metadata"], mbar_config)
    
    ddg = kj_to_kcal(dg_complex["dg"] - dg_free["dg"])
    single_results[mut_name] = {"ddg_calc": ddg, "ddg_exp": mut["ddg_exp"], "chain": chain}
    print(f"  ΔΔG: calc={ddg:.2f}, exp={mut['ddg_exp']:.1f}")
```

### Step 5: Stage 2 — Double-Mutation FEP + Coupling

```python
double_results = {}

for dm in double_mutations:
    m1, m2 = dm["pair"]
    pair_name = f"{m1}_{m2}"
    dm_dir = OUTPUT_DIR / "doubles" / pair_name
    dm_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n--- Double: {m1} + {m2} ---")
    
    # Combined double mutation
    combined = f"{m1}/{m2}"
    complex_dir = dm_dir / "complex"
    complex_dir.mkdir(exist_ok=True)
    fep_complex = run_fep_campaign(
        topology_path=str(solvated_pdb),
        mutation=combined,
        state="complex",
        fep_config=fep_config,
        output_dir=complex_dir,
    )
    dg_complex = solve_mbar(fep_complex["energy_files"], fep_complex["metadata"], mbar_config)
    
    # Free protein legs (both chains)
    chain1 = single_results[m1]["chain"]
    chain2 = single_results[m2]["chain"]
    free_pdb_1 = OUTPUT_DIR / f"solvated_{chain1}_free.pdb"
    free_pdb_2 = OUTPUT_DIR / f"solvated_{chain2}_free.pdb"
    
    free_dir = dm_dir / "free"
    free_dir.mkdir(exist_ok=True)
    # Handle same-chain vs cross-chain double mutations
    fep_free = run_fep_campaign(
        topology_path=str(free_pdb_1),
        mutation=combined,
        state="free",
        fep_config=fep_config,
        output_dir=free_dir,
    )
    dg_free = solve_mbar(fep_free["energy_files"], fep_free["metadata"], mbar_config)
    
    ddg_double = kj_to_kcal(dg_complex["dg"] - dg_free["dg"])
    
    # Coupling energy
    ddg_single_1 = single_results[m1]["ddg_calc"]
    ddg_single_2 = single_results[m2]["ddg_calc"]
    dddg = ddg_double - ddg_single_1 - ddg_single_2
    
    double_results[pair_name] = {
        "ddg_double": ddg_double,
        "ddg_single_1": ddg_single_1,
        "ddg_single_2": ddg_single_2,
        "dddg_calc": dddg,
        "dddg_exp": dm["dddg_exp"],
    }
    print(f"  ΔΔG_double={ddg_double:.2f}, ΔΔΔG={dddg:.2f} (exp={dm['dddg_exp']:.1f})")
```

### Step 6: Statistical Analysis + Classification

```python
# Single mutation correlation
ddg_calc_arr = np.array([single_results[m["name"]]["ddg_calc"] for m in all_single_mutations])
ddg_exp_arr = np.array([m["ddg_exp"] for m in all_single_mutations])
r_val, p_val = stats.pearsonr(ddg_exp_arr, ddg_calc_arr)
rmse = np.sqrt(np.mean((ddg_calc_arr - ddg_exp_arr)**2))

# Coupling energy accuracy
dddg_calc_arr = np.array([double_results[f"{dm['pair'][0]}_{dm['pair'][1]}"]["dddg_calc"]
                            for dm in double_mutations])
dddg_exp_arr = np.array([dm["dddg_exp"] for dm in double_mutations])

# Top coupled pairs
exp_ranked = sorted(double_mutations, key=lambda x: abs(x["dddg_exp"]), reverse=True)
calc_ranked = sorted(double_results.items(), key=lambda x: abs(x[1]["dddg_calc"]), reverse=True)
top5_exp = set(f"{dm['pair'][0]}_{dm['pair'][1]}" for dm in exp_ranked[:5])
top5_calc = set(calc_ranked[i][0] for i in range(min(5, len(calc_ranked))))
top5_correct = len(top5_exp & top5_calc)

print(f"R = {r_val:.3f}, RMSE = {rmse:.2f} kcal/mol")
print(f"Top-5 coupled pairs correct: {top5_correct}/5")

if r_val > 0.7 and rmse < 2.0 and top5_correct >= 5:
    classification = "PASS"
elif r_val > 0.5 and rmse < 3.0 and top5_correct >= 3:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-33", "feature_id": "F-33",
    "pearson_r": float(r_val),
    "rmse_kcal": float(rmse),
    "top5_coupled_correct": top5_correct,
    "single_mutations": single_results,
    "double_mutations": {k: v for k, v in double_results.items()},
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"EXP-33: R={r_val:.3f}, RMSE={rmse:.2f}, top5={top5_correct}/5 → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Single-Mutation ΔΔG Correlation

```python
fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(ddg_exp_arr, ddg_calc_arr, s=80, c='steelblue', edgecolors='black', zorder=5)

lims = [min(min(ddg_exp_arr), min(ddg_calc_arr)) - 1,
        max(max(ddg_exp_arr), max(ddg_calc_arr)) + 1]
ax.plot(lims, lims, 'k--', linewidth=1, alpha=0.5, label='y=x')

slope_r, intercept_r, _, _, _ = stats.linregress(ddg_exp_arr, ddg_calc_arr)
x_fit = np.linspace(lims[0], lims[1], 100)
ax.plot(x_fit, slope_r * x_fit + intercept_r, 'r-', linewidth=1.5,
         label=f'R={r_val:.3f}, RMSE={rmse:.2f}')

mut_names_list = [m["name"] for m in all_single_mutations]
for i, name in enumerate(mut_names_list):
    ax.annotate(name, (ddg_exp_arr[i], ddg_calc_arr[i]), textcoords="offset points",
                 xytext=(5, 5), fontsize=9)

ax.set_xlabel('ΔΔG_exp (kcal/mol)', fontsize=13)
ax.set_ylabel('ΔΔG_calc (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-33: Barnase–Barstar Single Mutations — {classification}', fontsize=14)
ax.legend(fontsize=11)
ax.set_aspect('equal')
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-33_single_correlation.png", dpi=300)
plt.close(fig)
```

### Figure 2: Coupling Energy (ΔΔΔG) Comparison

```python
fig, ax = plt.subplots(figsize=(12, 6))
pair_labels = [f"{dm['pair'][0]}\n{dm['pair'][1]}" for dm in double_mutations]
x = np.arange(len(pair_labels))
width = 0.35

bars1 = ax.bar(x - width/2, dddg_exp_arr, width, label='Experimental',
                color='#2ecc71', edgecolor='black')
bars2 = ax.bar(x + width/2, dddg_calc_arr, width, label='Computed (FEP)',
                color='steelblue', edgecolor='black')
ax.axhline(y=1.5, color='red', linestyle='--', alpha=0.5, label='Coupling threshold')
ax.axhline(y=0, color='black', linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels(pair_labels, fontsize=9)
ax.set_ylabel('ΔΔΔG (kcal/mol)', fontsize=12)
ax.set_title(f'EXP-33: Double-Mutant Cycle Coupling Energies — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-33_coupling_energy.png", dpi=300)
plt.close(fig)
```

### Figure 3: Interaction Network Map

```python
fig, ax = plt.subplots(figsize=(10, 8))
# Node positions: barnase residues on left, barstar on right
bn_names = [m["name"] for m in barnase_mutations]
bs_names = [m["name"] for m in barstar_mutations]

for i, name in enumerate(bn_names):
    ax.scatter(0.2, i, s=200, c='steelblue', edgecolors='black', zorder=5)
    ax.text(0.05, i, name, fontsize=11, ha='right', va='center')
for i, name in enumerate(bs_names):
    ax.scatter(0.8, i, s=200, c='coral', edgecolors='black', zorder=5)
    ax.text(0.95, i, name, fontsize=11, ha='left', va='center')

# Draw coupling lines
for pair_key, data_pair in double_results.items():
    m1, m2 = pair_key.split("_")
    if m1 in bn_names and m2 in bs_names:
        y1 = bn_names.index(m1)
        y2 = bs_names.index(m2)
        dddg = abs(data_pair["dddg_calc"])
        linewidth = max(0.5, min(5, dddg))
        color = 'red' if dddg > 1.5 else 'gray'
        ax.plot([0.2, 0.8], [y1, y2], color=color, linewidth=linewidth, alpha=0.6)
        ax.text(0.5, (y1 + y2) / 2, f'{dddg:.1f}', fontsize=8, ha='center',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-0.5, max(len(bn_names), len(bs_names)) - 0.5)
ax.set_title('EXP-33: Barnase–Barstar Coupling Network', fontsize=14)
ax.text(0.2, -0.4, 'Barnase', ha='center', fontsize=12, fontweight='bold', color='steelblue')
ax.text(0.8, -0.4, 'Barstar', ha='center', fontsize=12, fontweight='bold', color='coral')
ax.axis('off')
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-33_coupling_network.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-33: Barnase–Barstar Double-Mutant Cycles — Results Report

**Experiment ID:** EXP-33  **Feature ID:** F-33  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| Pearson R (singles) | [val] | > 0.7 | [P/M/F] |
| RMSE (singles) | [val] kcal/mol | < 2.0 | [P/M/F] |
| Top-5 coupled pairs | [n]/5 | 5/5 | [P/M/F] |
| Self-FEP control | [val] kcal/mol | < 0.3 | [P/F] |

## Figures
1. Single-mutation ΔΔG correlation
2. Coupling energy (ΔΔΔG) comparison
3. Interaction network map

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-33 is the most complex FEP experiment — 10 single-mutation + 5 double-mutation FEP campaigns, requiring 3–4 Colab sessions and meticulous staging.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-33"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures",
                "outputs/fep_singles", "outputs/fep_doubles"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

# Define mutation panels
BARNASE_MUTATIONS = ["K27A", "R59A", "R83A", "R87A", "H102A"]
BARSTAR_MUTATIONS = ["D35A", "D39A", "E76A", "Y29A", "W44A"]
SINGLE_MUTATIONS = BARNASE_MUTATIONS + BARSTAR_MUTATIONS  # 10 singles

DOUBLE_MUTATIONS = [
    ("K27A", "D35A"), ("R59A", "D39A"), ("R83A", "E76A"),
    ("R87A", "Y29A"), ("H102A", "W44A"),
]  # 5 doubles

# Create per-mutation output directories
for mut in SINGLE_MUTATIONS:
    (DRIVE_BASE / "outputs" / "fep_singles" / mut).mkdir(parents=True, exist_ok=True)
for m1, m2 in DOUBLE_MUTATIONS:
    (DRIVE_BASE / "outputs" / "fep_doubles" / f"{m1}_{m2}").mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

PROJECT_ROOT = Path("/content")
EXP_DIR = DRIVE_BASE
OUTPUT_DIR = DRIVE_BASE / "outputs"
FIGURES_DIR = DRIVE_BASE / "figures"

print(f"Colab environment ready. {len(SINGLE_MUTATIONS)} singles + "
      f"{len(DOUBLE_MUTATIONS)} doubles configured.")
print(f"Drive base: {DRIVE_BASE}")
```

### §5.2 GPU Platform Selection

```python
# §5.2 GPU Platform Selection
import openmm
from src.simulate.platform import select_platform

platform = select_platform("CUDA")
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
print(f"Platform: {platform.getName()}")

# EXP-33 runs 15 independent FEP campaigns (10 singles + 5 doubles).
# Each campaign: 20 lambda windows × 2 ns = 40 ns
# Total alchemical simulation: 15 × 40 ns = 600 ns
# CUDA acceleration essential — estimated ~30–60 hours on A100.
#
# Verified FEP function signatures (src/simulate/fep.py):
#   generate_lambda_schedule(n_windows: int) -> np.ndarray
#   create_alchemical_system(system, mutant_atom_indices, config) -> System
#   run_fep_window(alchemical_system, positions, lambda_value,
#                  lambda_schedule, config, output_dir) -> dict
#   run_fep_campaign(system, positions, mutant_atom_indices,
#                    config, output_dir) -> dict
#
# FEPConfig defaults: n_lambda_windows=20, per_window_duration_ns=2.0,
#   soft_core_alpha=0.5, temperature_k=310.0
```

### §5.3 Checkpoint and Resume Integration

```python
# §5.3 Checkpoint and Resume Integration
import json, time
from openmm import XmlSerializer
import openmm

CHECKPOINT_DIR = DRIVE_BASE / "checkpoints"

def save_checkpoint(simulation, phase_name, step=None, extra_data=None):
    """Save simulation state to Drive for session recovery."""
    ckpt_path = CHECKPOINT_DIR / f"{phase_name}.chk"
    simulation.saveCheckpoint(str(ckpt_path))
    state = simulation.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True)
    with open(CHECKPOINT_DIR / f"{phase_name}_state.xml", 'w') as f:
        f.write(XmlSerializer.serialize(state))
    progress = {
        "phase": phase_name, "step": step,
        "energy_kj": state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoules_per_mole),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    if extra_data:
        progress.update(extra_data)
    with open(CHECKPOINT_DIR / "progress.json", 'w') as f:
        json.dump(progress, f, indent=2)
    print(f"  Checkpoint saved: {phase_name}")

def load_checkpoint(simulation, phase_name):
    """Resume from Drive checkpoint. Returns True if loaded."""
    ckpt_path = CHECKPOINT_DIR / f"{phase_name}.chk"
    if not ckpt_path.exists():
        return False
    simulation.loadCheckpoint(str(ckpt_path))
    print(f"  Resumed from checkpoint: {phase_name}")
    return True

def get_completed_campaigns(campaign_type="fep_singles"):
    """Check which mutations have completed FEP campaigns."""
    completed = []
    base = DRIVE_BASE / "outputs" / campaign_type
    if not base.exists():
        return completed
    for d in base.iterdir():
        if d.is_dir():
            npz = d / "fep_campaign_results.npz"
            if npz.exists():
                completed.append(d.name)
    return completed

# ─── Multi-session staging strategy ───
#
# Session 1 — Barnase singles (5 campaigns × 40 ns = 200 ns):
#   K27A → R59A → R83A → R87A → H102A
#   ~10–20 hours on A100
#   After each: save_checkpoint(sim, f"fep_single_{mut}_complete")
#
# Session 2 — Barstar singles (5 campaigns × 40 ns = 200 ns):
#   D35A → D39A → E76A → Y29A → W44A
#   ~10–20 hours on A100
#
# Session 3 — Double-mutant FEP Stage 1 (3 doubles × 40 ns = 120 ns):
#   (K27A,D35A) → (R59A,D39A) → (R83A,E76A)
#   ~6–12 hours on A100
#   Requires double-mutation topology — run_fep_campaign with combined indices
#
# Session 4 — Double-mutant FEP Stage 2 + Analysis (2 doubles + analysis):
#   (R87A,Y29A) → (H102A,W44A) → ΔΔΔG calculation
#   ~4–8 hours + analysis
#
# Resume protocol:
#   completed_singles = get_completed_campaigns("fep_singles")
#   completed_doubles = get_completed_campaigns("fep_doubles")
#   print(f"Singles: {len(completed_singles)}/{len(SINGLE_MUTATIONS)}")
#   print(f"Doubles: {len(completed_doubles)}/{len(DOUBLE_MUTATIONS)}")
```

### §5.4 Progress Monitoring

```python
# §5.4 Progress Monitoring
import time, subprocess

def report_gpu_status():
    result = subprocess.run(
        ['nvidia-smi', '--query-gpu=name,memory.used,memory.total,utilization.gpu',
         '--format=csv,noheader'], capture_output=True, text=True)
    print(f"GPU: {result.stdout.strip()}")

def estimate_performance(simulation, n_steps=5000, dt_ps=0.002):
    t0 = time.time()
    simulation.step(n_steps)
    elapsed = time.time() - t0
    ns_per_day = (n_steps * dt_ps / 1000.0) / elapsed * 86400
    print(f"  Performance: {ns_per_day:.1f} ns/day")
    report_gpu_status()
    return ns_per_day

def report_dmc_progress():
    """Print overall DMC experiment progress."""
    completed_s = get_completed_campaigns("fep_singles")
    completed_d = get_completed_campaigns("fep_doubles")
    print(f"\n{'='*55}")
    print(f"EXP-33 Double-Mutant Cycle Progress")
    print(f"  Singles: {len(completed_s)}/{len(SINGLE_MUTATIONS)}")
    for mut in SINGLE_MUTATIONS:
        status = "✓" if mut in completed_s else "○"
        print(f"    {status} {mut}")
    print(f"  Doubles: {len(completed_d)}/{len(DOUBLE_MUTATIONS)}")
    for m1, m2 in DOUBLE_MUTATIONS:
        key = f"{m1}_{m2}"
        status = "✓" if key in completed_d else "○"
        print(f"    {status} {m1}+{m2}")
    print(f"{'='*55}\n")

# ─── Runtime estimates (A100 40GB) ───
# Per mutation FEP campaign (20 windows × 2 ns): ~2–4 hours
# Session 1 (5 barnase singles):   ~10–20 hours
# Session 2 (5 barstar singles):   ~10–20 hours
# Session 3 (3 doubles):           ~6–12 hours
# Session 4 (2 doubles + analysis): ~4–8 hours
# Total: ~30–60 hours (3–4 Colab sessions)
```

### §5.5 Results Persistence

```python
# §5.5 Results Persistence
import shutil

def sync_to_drive(local_dir, drive_subdir="outputs"):
    drive_dir = DRIVE_BASE / drive_subdir
    drive_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for f in Path(local_dir).rglob("*"):
        if f.is_file():
            dest = drive_dir / f.relative_to(local_dir)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dest)
            count += 1
    print(f"  Synced {count} files → {drive_dir}")

def verify_drive_sync():
    """Verify all FEP results exist for singles and doubles."""
    print("Single-mutation FEP results:")
    for mut in SINGLE_MUTATIONS:
        npz = DRIVE_BASE / "outputs" / "fep_singles" / mut / "fep_campaign_results.npz"
        status = "✓" if npz.exists() else "✗ MISSING"
        print(f"  {status} {mut}/fep_campaign_results.npz")
    print("\nDouble-mutation FEP results:")
    for m1, m2 in DOUBLE_MUTATIONS:
        key = f"{m1}_{m2}"
        npz = DRIVE_BASE / "outputs" / "fep_doubles" / key / "fep_campaign_results.npz"
        status = "✓" if npz.exists() else "✗ MISSING"
        print(f"  {status} {key}/fep_campaign_results.npz")
    # ΔΔΔG coupling results
    coupling = DRIVE_BASE / "outputs" / "coupling_energies.json"
    status = "✓" if coupling.exists() else "✗ MISSING"
    print(f"  {status} coupling_energies.json")
```

### §5.6 GPU-Optimized Figure Generation

```python
# §5.6 GPU-Optimized Figure Generation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save_figure(fig, name):
    fig_path = DRIVE_BASE / "figures" / f"{name}.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(fig_path), dpi=150, bbox_inches='tight')
    print(f"  Figure saved: {fig_path}")
    plt.close(fig)

# Key figures for EXP-33:
# 1. Single-mutation ΔΔG correlation scatter (barnase + barstar on same plot)
# 2. Double-mutation ΔΔG bar chart with error bars
# 3. ΔΔΔG coupling energies — comparison with experimental (Schreiber 1995)
# 4. Coupling network map (bipartite graph — barnase left, barstar right)
# 5. Thermodynamic cycle closure check (sum of legs per cycle)
```

### §5.7 Error Recovery Procedures

```python
# §5.7 Error Recovery Procedures

def gpu_safe_run(func, *args, max_retries=2, **kwargs):
    for attempt in range(max_retries + 1):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            with open(DRIVE_BASE / "error_log.txt", 'a') as f:
                f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} "
                        f"[attempt {attempt+1}] {type(e).__name__}: {e}\n")
            if "out of memory" in str(e).lower():
                properties['CudaPrecision'] = 'single'
            if attempt >= max_retries:
                raise

# ─── EXP-33–specific recovery ───
#
# 1. Session timeout between sessions:
#    - Each completed campaign saves .npz to Drive
#    - Use get_completed_campaigns() to identify remaining work
#    - Skip completed campaigns, continue from next mutation
#
# 2. Double-mutation topology failure:
#    - Double mutations require combined mutant_atom_indices for both residues
#    - If topology setup fails: rebuild system from PDB with both mutations applied
#    - Verify atom indexing matches between single and double FEP setups
#
# 3. ΔΔΔG calculation:
#    - ΔΔΔG(A,B) = ΔΔG(A+B) - ΔΔG(A) - ΔΔG(B)
#    - If |ΔΔΔG| > 5 kcal/mol, investigate individual campaign convergence
#    - Expect |ΔΔΔG| ≈ 1–4 kcal/mol for coupled hot-spot pairs
#
# 4. Thermodynamic cycle closure:
#    - For each cycle: ΔG(WT→A) + ΔG(A→AB) = ΔG(WT→B) + ΔG(B→AB)
#    - Closure error should be < 0.5 kcal/mol
#    - Large closure errors indicate sampling insufficiency in one leg
#
# 5. MBAR convergence failure for charged mutations (K27A, R59A, etc.):
#    - Charged→neutral transformations have large electrostatic ΔG
#    - Increase n_lambda_windows from 20 to 30 for charged mutations
#    - Consider separate electrostatic + vdW staging
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
