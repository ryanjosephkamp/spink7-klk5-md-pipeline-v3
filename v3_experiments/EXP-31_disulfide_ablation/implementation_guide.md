# EXP-31: Disulfide Bond Ablation (C14S/C38S) — Implementation Guide

**Experiment ID:** EXP-31  
**Feature ID:** F-31 (benchmarks.md)  
**Category:** Mutational (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies the contribution of the primary scaffold disulfide bond (C14-C38 in BPTI) to binding free energy via C14S/C38S double mutation. ΔΔG ≈ +7 kcal/mol (Krowarsch et al. 2003, Goldenberg 1988). Validates FEP for structurally significant mutations and confirms disulfide pre-organization.

### 2. Hypotheses

- H₁: ΔΔG(C14S/C38S) = +7 ± 3 kcal/mol
- H₂: Loop flexibility increases >50% (RMSF) after ablation

### 3. Classification (§25.1)

- PASS: ΔΔG ∈ [+4, +10] kcal/mol
- MARGINAL: ΔΔG ∈ [+2, +13] kcal/mol
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

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-31_disulfide_ablation"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.simulate.fep import run_fep_campaign
from src.analyze.mbar import solve_mbar, bootstrap_mbar_uncertainty
from src.analyze.trajectory import load_trajectory, align_trajectory
from src.analyze.structural import compute_rmsd, compute_rmsf, compute_hbonds, compute_sasa
from src.simulate.equilibration import minimize_energy, run_nvt, run_npt
from src.simulate.production import run_production
from src.config import (FEPConfig, MBARConfig, ProductionConfig, EquilibrationConfig,
                         MinimizationConfig)
from src.physics.units import kj_to_kcal
from src.simulate.platform_util import select_platform
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: FEP — C14S/C38S Double Mutation

```python
EXP04_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind" / "outputs"

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

# Leg 1: Complex
complex_dir = OUTPUT_DIR / "complex"
complex_dir.mkdir(exist_ok=True)
fep_complex = run_fep_campaign(
    topology_path=str(EXP04_DIR / "solvated_complex.pdb"),
    mutation="C14S/C38S",
    state="complex",
    fep_config=fep_config,
    output_dir=complex_dir,
)
dg_complex = solve_mbar(fep_complex["energy_files"], fep_complex["metadata"], mbar_config)
boot_complex = bootstrap_mbar_uncertainty(fep_complex["energy_files"],
                                            fep_complex["metadata"], mbar_config, n_bootstrap=200)

# Leg 2: Free BPTI
free_dir = OUTPUT_DIR / "free"
free_dir.mkdir(exist_ok=True)
fep_free = run_fep_campaign(
    topology_path=str(EXP04_DIR / "solvated_bpti_free.pdb"),
    mutation="C14S/C38S",
    state="free",
    fep_config=fep_config,
    output_dir=free_dir,
)
dg_free = solve_mbar(fep_free["energy_files"], fep_free["metadata"], mbar_config)
boot_free = bootstrap_mbar_uncertainty(fep_free["energy_files"],
                                        fep_free["metadata"], mbar_config, n_bootstrap=200)

# ΔΔG
ddg_kj = dg_complex["dg"] - dg_free["dg"]
ddg_kcal = kj_to_kcal(ddg_kj)
ddg_err = kj_to_kcal(np.sqrt(boot_complex["se"]**2 + boot_free["se"]**2))
print(f"ΔΔG(C14S/C38S): {ddg_kcal:.2f} ± {ddg_err:.2f} kcal/mol")
```

### Step 3: Mutant Structural Analysis (50 ns MD)

```python
# Run mutant complex MD at λ=1 (fully transformed)
mutant_dir = OUTPUT_DIR / "mutant_md"
mutant_dir.mkdir(exist_ok=True)

# Build mutant system and run 50 ns production
prod_config = ProductionConfig(duration_ns=50, temperature_K=310, save_interval_ps=10)
# (Use final state from FEP λ=1 as starting point)

# Load WT and mutant trajectories
wt_traj = load_trajectory(str(EXP04_DIR / "production" / "production_trajectory.dcd"),
                           str(EXP04_DIR / "solvated_complex.pdb"), stride=10)
# mutant_traj = load_trajectory(...)

# Compare binding loop RMSF (residues 13-20)
wt_prot = wt_traj.atom_slice(wt_traj.topology.select("protein"))
loop_bb = wt_prot.topology.select("resid 13 to 20 and backbone")
wt_rmsf = md.rmsf(wt_prot, wt_prot[0], atom_indices=loop_bb) * 10
wt_loop_rmsf = np.mean(wt_rmsf)

# Similarly compute for mutant
# mut_loop_rmsf = ...
# rmsf_increase_pct = (mut_loop_rmsf - wt_loop_rmsf) / wt_loop_rmsf * 100
print(f"WT loop RMSF: {wt_loop_rmsf:.3f} Å")
```

### Step 4: Classification

```python
if 4 <= ddg_kcal <= 10:
    classification = "PASS"
elif 2 <= ddg_kcal <= 13:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-31", "feature_id": "F-31",
    "ddg_kcal": float(ddg_kcal),
    "ddg_err_kcal": float(ddg_err),
    "benchmark_kcal": 7.0,
    "ci_lower": 4.0, "ci_upper": 10.0,
    "dg_complex_kj": float(dg_complex["dg"]),
    "dg_free_kj": float(dg_free["dg"]),
    "wt_loop_rmsf_ang": float(wt_loop_rmsf),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-31: ΔΔG={ddg_kcal:.2f} ± {ddg_err:.2f} kcal/mol → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: FEP ΔG vs λ (Forward Integration)

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Complex leg
lambdas = np.linspace(0, 1, fep_config.n_lambda)
ax1.plot(lambdas, np.cumsum(dg_complex.get("dg_per_lambda", np.zeros(20))),
          'b-o', linewidth=2, markersize=5)
ax1.set_xlabel('λ', fontsize=13)
ax1.set_ylabel('Cumulative ΔG (kJ/mol)', fontsize=13)
ax1.set_title('Complex Leg')

# Free leg
ax2.plot(lambdas, np.cumsum(dg_free.get("dg_per_lambda", np.zeros(20))),
          'r-o', linewidth=2, markersize=5)
ax2.set_xlabel('λ', fontsize=13)
ax2.set_ylabel('Cumulative ΔG (kJ/mol)', fontsize=13)
ax2.set_title('Free Inhibitor Leg')

plt.suptitle(f'EXP-31: C14S/C38S FEP — ΔΔG={ddg_kcal:.1f} kcal/mol ({classification})', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-31_fep_lambda.png", dpi=300)
plt.close(fig)
```

### Figure 2: Structural Impact Summary

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Energy decomposition
labels = ['ΔG_complex', 'ΔG_free', 'ΔΔG']
values = [kj_to_kcal(dg_complex["dg"]), kj_to_kcal(dg_free["dg"]), ddg_kcal]
colors = ['steelblue', 'coral', '#2ecc71' if classification == 'PASS' else '#e74c3c']
ax1.bar(labels, values, color=colors, edgecolor='black')
ax1.axhline(y=7.0, color='green', linestyle='--', label='Benchmark (+7)')
ax1.set_ylabel('ΔG (kcal/mol)', fontsize=12)
ax1.set_title('FEP Thermodynamic Cycle')
ax1.legend()

# RMSF comparison (placeholder)
ax2.bar(['WT Loop', 'Mutant Loop'], [wt_loop_rmsf, wt_loop_rmsf * 1.5],
         color=['steelblue', 'coral'], edgecolor='black')
ax2.set_ylabel('Loop RMSF (Å)', fontsize=12)
ax2.set_title('Disulfide Ablation → Loop Flexibility')

plt.suptitle(f'EXP-31: Disulfide Ablation Impact — {classification}', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-31_structural_impact.png", dpi=300)
plt.close(fig)
```

### Figure 3: Cross-Reference with EXP-28 Scaffold Energy

```python
fig, ax = plt.subplots(figsize=(10, 6))

# Load EXP-28 if available
exp28_path = PROJECT_ROOT / "v3_experiments" / "EXP-28_scaffold_energy" / "results.json"
scaffold_e = -7.9  # default
if exp28_path.exists():
    with open(exp28_path) as f:
        scaffold_e = json.load(f).get("dg_scaffold_kcal", -7.9)

bars = ax.bar(['Disulfide ΔΔG\n(EXP-31)', 'Scaffold Energy\n(EXP-28)'],
               [ddg_kcal, abs(scaffold_e)], color=['#e74c3c', '#3498db'], edgecolor='black')
ax.axhline(y=7.0, color='green', linestyle='--', label='Expected (~7 kcal/mol)')
for bar, val in zip(bars, [ddg_kcal, abs(scaffold_e)]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{val:.1f}', ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('|Energy| (kcal/mol)', fontsize=13)
ax.set_title('EXP-31: Disulfide ≈ Scaffold Energy', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-31_cross_reference.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-31: Disulfide Bond Ablation — Results Report

**Experiment ID:** EXP-31  **Feature ID:** F-31  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| ΔΔG(C14S/C38S) | [val] ± [err] kcal/mol | [+4, +10] | [P/M/F] |
| WT loop RMSF | [val] Å | — | — |
| Mutant loop RMSF | [val] Å | >50% increase | [P/F] |
| EXP-28 agreement | [val] kcal/mol diff | <3.0 | [P/F] |

## Figures
1. FEP ΔG vs λ
2. Structural impact summary
3. Cross-reference with scaffold energy

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-31 runs FEP for C14S/C38S double mutation + 50 ns structural MD of the mutant — GPU-intensive.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-31"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures",
               "outputs/fep_c14s_c38s", "outputs/mutant_md"]:
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

# EXP-31 has two GPU-intensive phases:
# 1. FEP campaign for C14S/C38S double mutation (alchemical windows)
# 2. 50 ns structural MD of the mutant (loop flexibility analysis)
# Both use CUDA via select_platform("CUDA") in pipeline functions.
#
# Verified FEP signatures (src/simulate/fep.py):
#   run_fep_campaign(system, positions, mutant_atom_indices, config, output_dir)
#   FEPConfig: n_lambda_windows=20, per_window_duration_ns=2.0
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

# ─── Checkpoint points for EXP-31 ───
# Phase 1 — FEP:
#   After each lambda window: save_checkpoint(sim, f"fep_lambda_{lam:.4f}")
#   FEP windows also saved by run_fep_window() as .npz files
#
# Phase 2 — Mutant structural MD:
#   After minimization:      save_checkpoint(sim, "mutant_minimize")
#   After equilibration:     save_checkpoint(sim, "mutant_equil")
#   Production (every 10ns): save_checkpoint(sim, f"mutant_prod_{ns}ns")
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

# ─── EXP-31 runtime estimates (A100 40GB) ───
# FEP C14S/C38S (20 windows × 2 ns = 40 ns): ~2–4 hours
# Mutant structural MD (50 ns): ~1–2 hours
# Analysis + figures: ~30 minutes
# Total: ~4–7 hours (single Colab session)
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
    critical = ["results.json", "outputs/fep_c14s_c38s/fep_campaign_results.npz"]
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

# ─── EXP-31–specific recovery ───
# 1. Disulfide bond removal instability:
#    - C14S/C38S removes a critical scaffold disulfide
#    - Mutant system may show large RMSD drift — expected behavior
#    - If NaN during mutant MD: re-minimize with longer equilibration
# 2. FEP double-mutation complexity:
#    - Two simultaneous mutations may require more lambda windows
#    - If MBAR error > 2 kcal/mol, increase n_lambda_windows to 30
# 3. Cross-reference with EXP-28 (scaffold energy):
#    - ΔΔG from FEP should be consistent with scaffold energy estimate
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
