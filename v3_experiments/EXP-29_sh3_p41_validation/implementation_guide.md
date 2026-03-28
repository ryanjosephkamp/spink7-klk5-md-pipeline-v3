# EXP-29: SH3–p41 Method Validation (ΔG_bind) — Implementation Guide

**Experiment ID:** EXP-29  
**Feature ID:** F-29 (benchmarks.md)  
**Category:** Biophysical / Method Validation (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates the pipeline's free-energy methods on Fyn SH3–p41 (PDB 4EIK, ΔG_bind = −7.99 kcal/mol). An out-of-family positive control: protein-peptide interaction with different binding geometry than protease-inhibitor systems. Both US/WHAM and SMD/Jarzynski must independently reproduce ΔG within CI [−12.7, −3.3].

### 2. Hypotheses

- H₁: US/WHAM ΔG within [−12.7, −3.3] kcal/mol
- H₂: SMD/Jarzynski ΔG within same CI
- H₃: Both methods agree within respective uncertainties

### 3. Classification (§25.1)

- PASS: [−12.7, −3.3] kcal/mol
- MARGINAL: [−16.0, −1.0] kcal/mol
- FAIL: Outside marginal range or wrong sign

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-29_sh3_p41_validation"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.prep.pdb_fetch import fetch_pdb
from src.prep.structure_prep import clean_structure, assign_protonation
from src.prep.topology import build_topology, solvate_system
from src.simulate.equilibration import minimize_energy, run_nvt, run_npt
from src.simulate.production import run_production
from src.simulate.smd import run_smd_campaign
from src.simulate.umbrella import run_umbrella_campaign, generate_window_centers
from src.analyze.wham import solve_wham, bootstrap_pmf_uncertainty
from src.analyze.mbar import solve_mbar, bootstrap_mbar_uncertainty
from src.analyze.jarzynski import jarzynski_free_energy, diagnose_dissipation
from src.config import (SystemConfig, MinimizationConfig, EquilibrationConfig,
                         ProductionConfig, SMDConfig, UmbrellaConfig, WHAMConfig, MBARConfig)
from src.simulate.platform_util import select_platform
from src.physics.units import kj_to_kcal
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Structure Preparation

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

pdb_path = fetch_pdb("4EIK", data_dir)
clean_path = clean_structure(pdb_path, chains_to_keep=None,
                              remove_heteroatoms=True, remove_waters=True)
prot_path = assign_protonation(clean_path, ph=7.4, force_field="AMBER", use_propka=True)
print(f"Prepared: {prot_path}")
```

### Step 3: Build System + Equilibrate

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
print(f"Solvated: {n_waters} waters")

from openmm.app import PDBFile
solvated_pdb = OUTPUT_DIR / "solvated_sh3_p41.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

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

### Step 4: SMD Campaign (50 replicates)

```python
smd_dir = OUTPUT_DIR / "smd"
smd_dir.mkdir(exist_ok=True)

smd_config = SMDConfig(
    spring_constant_kj=1000.0,
    velocity_nm_ps=0.001,
    pull_distance_nm=3.0,
    n_replicates=50,
    save_interval_ps=1.0,
)
smd_result = run_smd_campaign(simulation, smd_config, smd_dir,
                               group1_selection="chainid 0",
                               group2_selection="chainid 1")

# Jarzynski analysis
work_values = np.array(smd_result["work_values"])
jarz = jarzynski_free_energy(work_values, T=310.0)
dg_smd_kj = jarz["dg_jarzynski"]
dg_smd_kcal = kj_to_kcal(dg_smd_kj)
dissipation = diagnose_dissipation(work_values, T=310.0)
print(f"SMD ΔG (Jarzynski): {dg_smd_kcal:.2f} kcal/mol")
print(f"Dissipation: {dissipation}")
```

### Step 5: US/WHAM Campaign

```python
us_dir = OUTPUT_DIR / "umbrella"
us_dir.mkdir(exist_ok=True)

us_config = UmbrellaConfig(
    xi_start_nm=1.5, xi_end_nm=4.0, xi_spacing_nm=0.05,
    spring_constant_kj=1000.0,
    per_window_ns=10.0, equilibration_ps=200,
    pre_position_velocity_nm_ps=0.01,
    save_interval_ps=1.0,
)
us_result = run_umbrella_campaign(simulation, us_config, us_dir,
                                   group1_selection="chainid 0",
                                   group2_selection="chainid 1")

# WHAM
wham_config = WHAMConfig(tolerance=1e-6, max_iterations=100000, n_bins=200)
pmf, bin_centers = solve_wham(us_result["xi_files"], us_result["metadata"], wham_config)
pmf_boot, ci_lower, ci_upper = bootstrap_pmf_uncertainty(
    us_result["xi_files"], us_result["metadata"], wham_config, n_bootstrap=200)

# ΔG extraction
dg_wham_kj = float(np.min(pmf) - pmf[-1])
dg_wham_kcal = kj_to_kcal(dg_wham_kj)

# MBAR cross-check
mbar_config = MBARConfig(solver="robust", tolerance=1e-7, max_iterations=10000)
mbar_result = solve_mbar(us_result["xi_files"], us_result["metadata"], mbar_config)
dg_mbar_kcal = kj_to_kcal(mbar_result["dg"])

print(f"WHAM ΔG: {dg_wham_kcal:.2f} kcal/mol")
print(f"MBAR ΔG: {dg_mbar_kcal:.2f} kcal/mol")
```

### Step 6: Classification

```python
# Method agreement
method_agree = abs(dg_wham_kcal - dg_smd_kcal) < 3.0
wham_mbar_agree = abs(dg_wham_kcal - dg_mbar_kcal) < 1.5

if -12.7 <= dg_wham_kcal <= -3.3:
    classification = "PASS"
elif -16.0 <= dg_wham_kcal <= -1.0:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-29", "feature_id": "F-29",
    "dg_wham_kcal": float(dg_wham_kcal),
    "dg_mbar_kcal": float(dg_mbar_kcal),
    "dg_smd_kcal": float(dg_smd_kcal),
    "dg_exp_kcal": -7.99,
    "ci_lower": -12.7, "ci_upper": -3.3,
    "method_agreement": method_agree,
    "wham_mbar_agreement": wham_mbar_agree,
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-29: WHAM={dg_wham_kcal:.2f}, SMD={dg_smd_kcal:.2f} → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: PMF Profile

```python
fig, ax = plt.subplots(figsize=(12, 6))
bc_ang = bin_centers * 10
ax.plot(bc_ang, pmf_boot, 'b-', linewidth=2, label='WHAM PMF')
ax.fill_between(bc_ang, ci_lower, ci_upper, alpha=0.2, color='blue')
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.set_xlabel('ξ (Å)', fontsize=13)
ax.set_ylabel('PMF (kJ/mol)', fontsize=13)
ax.set_title(f'EXP-29: SH3–p41 PMF — {classification}', fontsize=14)
ax.legend(fontsize=12)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-29_pmf.png", dpi=300)
plt.close(fig)
```

### Figure 2: SMD Work Distribution

```python
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(work_values, bins=25, color='coral', edgecolor='black', density=True)
ax.axvline(x=np.mean(work_values), color='red', linestyle='--', linewidth=2,
            label=f'⟨W⟩={np.mean(work_values):.1f} kJ/mol')
ax.axvline(x=dg_smd_kj, color='blue', linestyle='--', linewidth=2,
            label=f'ΔG_Jarz={dg_smd_kcal:.1f} kcal/mol')
ax.set_xlabel('Work (kJ/mol)', fontsize=13)
ax.set_ylabel('Density', fontsize=13)
ax.set_title('EXP-29: SMD Work Distribution (50 replicates)', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-29_smd_work.png", dpi=300)
plt.close(fig)
```

### Figure 3: Method Comparison Bar Chart

```python
fig, ax = plt.subplots(figsize=(10, 6))
methods = ['WHAM', 'MBAR', 'SMD/Jarz', 'Experiment']
vals = [dg_wham_kcal, dg_mbar_kcal, dg_smd_kcal, -7.99]
colors = ['steelblue', '#3498db', 'coral', '#2ecc71']
bars = ax.bar(methods, vals, color=colors, edgecolor='black')
ax.axhspan(-12.7, -3.3, alpha=0.1, color='green', label='PASS CI')
for bar, val in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 0.5,
             f'{val:.2f}', ha='center', fontsize=11, fontweight='bold', color='white')
ax.set_ylabel('ΔG_bind (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-29: Method Comparison — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-29_method_comparison.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-29: SH3–p41 Method Validation — Results Report

**Experiment ID:** EXP-29  **Feature ID:** F-29  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| ΔG WHAM | [val] kcal/mol | [−12.7, −3.3] | [P/M/F] |
| ΔG MBAR | [val] kcal/mol | — | — |
| ΔG SMD | [val] kcal/mol | [−12.7, −3.3] | [P/M/F] |
| ΔG exp | −7.99 kcal/mol | Isvoran 2018 | — |
| WHAM-MBAR agree | [val] kcal/mol | <1.5 | [P/F] |
| WHAM-SMD agree | [val] kcal/mol | <3.0 | [P/F] |

## Figures
1. PMF profile
2. SMD work distribution
3. Method comparison bar chart

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps — they do not replace them.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-29"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures", "outputs/structures",
               "outputs/equilibration", "outputs/production",
               "outputs/smd", "outputs/umbrella"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

PROJECT_ROOT = Path("/content")
EXP_DIR = DRIVE_BASE
OUTPUT_DIR = DRIVE_BASE / "outputs"
FIGURES_DIR = DRIVE_BASE / "figures"

print(f"Colab environment ready. Drive base: {DRIVE_BASE}")
```

### §5.2 GPU Platform Selection

```python
# §5.2 GPU Platform Selection
import openmm
from src.simulate.platform import select_platform

platform = select_platform("CUDA")
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
print(f"Platform: {platform.getName()}")

# EXP-29 uses dual method (SMD + US/WHAM) — both benefit from CUDA.
# Update all Simulation() calls to include platform and properties.
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

# ─── Checkpoint points for EXP-29 ───
# After minimization:       save_checkpoint(sim, "minimize")
# After NVT/NPT:            save_checkpoint(sim, "nvt"), save_checkpoint(sim, "npt")
# Production (every 10ns):  save_checkpoint(sim, f"production_{ns}ns")
# After each SMD rep:        save_checkpoint(sim, f"smd_rep_{i:03d}")
# After each US window:      save_checkpoint(sim, f"us_window_{i:03d}")
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

# ─── EXP-29 runtime estimates (A100 40GB) ───
# Production (50 ns):   ~1–2 hours
# SMD (30 reps × 3 ns): ~2–4 hours
# US (51 windows × 10 ns): ~10–20 hours
# Total: ~14–26 hours (within single Colab session)
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
    critical = ["results.json", "checkpoints/progress.json"]
    for f in critical:
        path = DRIVE_BASE / f
        status = "✓" if path.exists() else "✗ MISSING"
        print(f"  {status} {f}")
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

# ─── EXP-29–specific recovery ───
# 1. Dual-method disagreement (SMD vs WHAM > 3 kcal/mol):
#    - Not a GPU error — indicates method limitation
#    - Document both values; classify based on WHAM (more reliable)
# 2. SH3-p41 complex instability:
#    - Weak binding (~−8 kcal/mol) — complex may dissociate during production
#    - If RMSD > 5 Å, reduce production time or add restraints
# 3. Session timeout: resume from last checkpoint
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
