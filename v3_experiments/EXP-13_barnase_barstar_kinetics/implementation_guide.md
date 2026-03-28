# EXP-13: Barnase–Barstar Association Kinetics — Implementation Guide

**Experiment ID:** EXP-13  
**Feature ID:** F-13 (benchmarks.md)  
**Category:** Kinetic (Semi-Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates pipeline kinetics capabilities using barnase–barstar, a classic electrostatically-steered protein–protein association system. PDB 1BRS (2.0 Å). kon,exp = 3.7 × 10⁸ M⁻¹s⁻¹ (Schreiber & Fersht 1996). Uses PMF-based diffusion theory (Szabo-Shoup-Northrup and Smoluchowski integral) to predict kon. Electrostatic steering control: charge-neutralized proteins should reduce kon by ≥2 orders of magnitude.

### 2. Benchmark Values

| Parameter | Value | Source |
|-----------|-------|--------|
| kon | 3.7 × 10⁸ M⁻¹s⁻¹ | Schreiber & Fersht 1996 |
| Kd | ~10⁻¹⁴ M | Schreiber & Fersht 1996 |
| ΔG_bind | ~-19 kcal/mol | Derived |
| PASS range | 10⁷ – 10¹⁰ M⁻¹s⁻¹ | §25.1 semi-quantitative |

### 3. Classification Criteria (§25.1)

- PASS: kon within 10⁷–10¹⁰ M⁻¹s⁻¹ AND electrostatic control ≥ 2 orders reduction
- MARGINAL: kon within 10⁶–10¹¹ but control shows < 2 orders
- FAIL: kon outside 10⁶–10¹¹ or no electrostatic steering effect

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-13_barnase_barstar_kinetics"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import (SystemConfig, MinimizationConfig, EquilibrationConfig,
                         ProductionConfig, UmbrellaConfig, WHAMConfig)
from src.prep.pdb_fetch import fetch_pdb
from src.prep.pdb_clean import clean_structure
from src.prep.protonate import assign_protonation
from src.prep.topology import build_topology
from src.prep.solvate import solvate_system
from src.simulate.minimizer import minimize_energy
from src.simulate.equilibrate import run_nvt, run_npt
from src.simulate.production import run_production
from src.simulate.umbrella import run_umbrella_campaign, generate_window_centers
from src.simulate.platform import select_platform
from src.analyze.wham import solve_wham, bootstrap_pmf_uncertainty
from src.analyze.mbar_analysis import solve_mbar
from src.analyze.trajectory import load_trajectory
from src.physics.collective_variables import com_distance
from src.physics.units import kj_to_kcal, nm_to_angstrom, kbt
import matplotlib.pyplot as plt

print("All imports successful.")
```

### Step 2: Fetch and Prepare Barnase–Barstar Complex

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

pdb_file = fetch_pdb("1BRS", data_dir)

# 1BRS contains barnase (chain A) and barstar (chain D) — check chains
import mdtraj as md
structure = md.load(str(pdb_file))
chains = [c.chain_id for c in structure.topology.chains]
print(f"Chains in 1BRS: {set(chains)}")

# Clean: keep barnase + barstar chains
complex_clean = clean_structure(pdb_file, chains_to_keep=["A", "D"],
                                 remove_heteroatoms=True, remove_waters=True, model_index=1)

# Protonate at pH 7.4
complex_protonated = assign_protonation(complex_clean, ph=7.4, force_field="AMBER", use_propka=True)
```

### Step 3: Solvate and Equilibrate

```python
from openmm.app import PME, PDBFile, Simulation
from openmm import LangevinMiddleIntegrator
import openmm

sys_config = SystemConfig()
topology, system, modeller = build_topology(complex_protonated, sys_config,
                                             nonbonded_method=PME, nonbonded_cutoff_nm=1.0)
modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)

solvated_pdb = OUTPUT_DIR / "solvated_complex.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

platform = select_platform()
integrator = LangevinMiddleIntegrator(300*openmm.unit.kelvin, 1.0/openmm.unit.picosecond,
                                      0.002*openmm.unit.picoseconds)
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

minimize_energy(sim, MinimizationConfig())
eq_dir = OUTPUT_DIR / "equilibration"
eq_dir.mkdir(exist_ok=True)
run_nvt(sim, EquilibrationConfig(), eq_dir)
run_npt(sim, EquilibrationConfig(), eq_dir)
print("Equilibration complete.")
```

### Step 4: US/WHAM for PMF

```python
umbrella_config = UmbrellaConfig()
wham_config = WHAMConfig()

us_dir = OUTPUT_DIR / "umbrella"
us_dir.mkdir(exist_ok=True)

window_centers = generate_window_centers(umbrella_config)
us_results = run_umbrella_campaign(sim, umbrella_config, us_dir,
                                    cv_type="com_distance",
                                    group1_sel="chainid 0",
                                    group2_sel="chainid 1")

pmf_xi, pmf_G, pmf_err = solve_wham(us_results, wham_config)
boot_err = bootstrap_pmf_uncertainty(us_results, wham_config, n_bootstrap=200)

np.savez(OUTPUT_DIR / "pmf_data.npz", xi=pmf_xi, G=pmf_G, err=boot_err)
print(f"PMF computed: {len(pmf_xi)} bins, min={np.min(pmf_G):.2f}, max barrier={np.max(pmf_G):.2f}")
```

### Step 5: Compute kon via Diffusion Theory

```python
# Physical constants
kB = 1.987e-3  # kcal/(mol·K)
T = 300  # K
beta = 1.0 / (kB * T)

# Relative diffusion coefficient (typical for protein-protein)
D_rel = 10e-6  # cm²/s = 10 × 10⁻⁶ cm²/s (typical protein-protein)
D_rel_nm2_ps = D_rel * 1e14 * 1e-12  # Convert to nm²/ps

xi_nm = pmf_xi  # Already in nm from pipeline
W_kcal = pmf_G  # Already in kcal/mol

# Method A: Szabo-Shoup-Northrup (spherical, effective)
# k_on = 4π D_rel R* exp(-W_min/kT) for capture radius R*
R_star = xi_nm[np.argmin(W_kcal)]  # Minimum PMF position
W_min = np.min(W_kcal)
k_on_ssn = 4 * np.pi * D_rel * R_star * 1e-7 * np.exp(-beta * W_min) * 6.022e23
# Convert units: D in cm²/s, R in cm, result in M⁻¹s⁻¹
R_star_cm = R_star * 1e-7
k_on_ssn = 4 * np.pi * D_rel * R_star_cm * 6.022e23

print(f"Method A (SSN): kon = {k_on_ssn:.2e} M⁻¹s⁻¹")

# Method B: Smoluchowski Integral
# 1/k_on = ∫[R_min to R_max] exp(W(r)/kT) / (4π r² D(r)) dr
dr = np.diff(xi_nm)
integrand = np.exp(beta * W_kcal[:-1]) / (4 * np.pi * (xi_nm[:-1] * 1e-7)**2 * D_rel)
integral = np.sum(integrand * dr * 1e-7)
k_on_smoluchowski = 6.022e23 / integral  # M⁻¹s⁻¹

print(f"Method B (Smoluchowski): kon = {k_on_smoluchowski:.2e} M⁻¹s⁻¹")

kon_avg = np.sqrt(k_on_ssn * k_on_smoluchowski)  # Geometric mean
print(f"Geometric mean: kon = {kon_avg:.2e} M⁻¹s⁻¹")
```

### Step 6: Electrostatic Steering Control

```python
# Neutralize charges on both proteins and repeat US/WHAM
# This simulates removing electrostatic steering
ctrl_dir = OUTPUT_DIR / "control_neutral"
ctrl_dir.mkdir(exist_ok=True)

# Build system with zeroed protein charges
# NOTE: This is done by modifying partial charges in the OpenMM system
from copy import deepcopy
system_ctrl = deepcopy(system)
forces = {type(f).__name__: f for f in system_ctrl.getForces()}
nb_force = forces.get('NonbondedForce')

if nb_force:
    # Get protein atom indices
    topo_md = md.load(str(solvated_pdb))
    protein_atoms = topo_md.topology.select("protein")
    for idx in protein_atoms:
        charge, sigma, epsilon = nb_force.getParticleParameters(int(idx))
        nb_force.setParticleParameters(int(idx), 0.0 * charge.unit, sigma, epsilon)

# Re-run US/WHAM with neutralized charges
integrator_ctrl = LangevinMiddleIntegrator(300*openmm.unit.kelvin, 1.0/openmm.unit.picosecond,
                                            0.002*openmm.unit.picoseconds)
sim_ctrl = Simulation(modeller.topology, system_ctrl, integrator_ctrl, platform)
sim_ctrl.context.setPositions(modeller.positions)
minimize_energy(sim_ctrl, MinimizationConfig())

us_ctrl_results = run_umbrella_campaign(sim_ctrl, umbrella_config, ctrl_dir,
                                         cv_type="com_distance",
                                         group1_sel="chainid 0",
                                         group2_sel="chainid 1")
pmf_xi_ctrl, pmf_G_ctrl, _ = solve_wham(us_ctrl_results, wham_config)

# Compute kon for neutralized system
W_ctrl = pmf_G_ctrl
integrand_ctrl = np.exp(beta * W_ctrl[:-1]) / (4 * np.pi * (pmf_xi_ctrl[:-1] * 1e-7)**2 * D_rel)
integral_ctrl = np.sum(integrand_ctrl * np.diff(pmf_xi_ctrl) * 1e-7)
k_on_ctrl = 6.022e23 / integral_ctrl

steering_ratio = np.log10(kon_avg) - np.log10(k_on_ctrl)
print(f"Neutral control: kon = {k_on_ctrl:.2e} M⁻¹s⁻¹")
print(f"Electrostatic steering effect: {steering_ratio:.1f} orders of magnitude")
```

### Step 7: Classification

```python
kon_exp = 3.7e8
log_kon = np.log10(kon_avg)

if 7 <= log_kon <= 10 and steering_ratio >= 2:
    classification = "PASS"
elif 6 <= log_kon <= 11:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-13", "feature_id": "F-13",
    "kon_ssn": float(k_on_ssn), "kon_smoluchowski": float(k_on_smoluchowski),
    "kon_avg": float(kon_avg), "kon_exp": float(kon_exp),
    "log10_kon_pred": float(log_kon), "log10_kon_exp": np.log10(kon_exp),
    "kon_neutral_ctrl": float(k_on_ctrl),
    "electrostatic_steering_orders": float(steering_ratio),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-13: log10(kon) = {log_kon:.2f}, steering = {steering_ratio:.1f} orders → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: PMF with Wild-Type vs Neutralized Control

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 7))
ax.plot(pmf_xi * 10, pmf_G, 'b-', linewidth=2, label='Wild-type (charged)')
ax.fill_between(pmf_xi * 10, pmf_G - boot_err, pmf_G + boot_err, alpha=0.2, color='blue')
ax.plot(pmf_xi_ctrl * 10, pmf_G_ctrl, 'r--', linewidth=2, label='Charge-neutralized control')
ax.set_xlabel('COM Distance (Å)', fontsize=14)
ax.set_ylabel('PMF (kcal/mol)', fontsize=14)
ax.set_title('EXP-13: Barnase–Barstar PMF — Electrostatic Steering', fontsize=15)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-13_pmf_wt_vs_neutral.png", dpi=300)
plt.close(fig)
```

### Figure 2: kon Comparison

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
labels = ['Experiment', 'SSN', 'Smoluchowski', 'Geom. Mean', 'Neutral Ctrl']
log_values = [np.log10(kon_exp), np.log10(k_on_ssn), np.log10(k_on_smoluchowski),
              log_kon, np.log10(k_on_ctrl)]
colors = ['green', 'steelblue', 'steelblue', 'navy', 'red']
ax.barh(labels, log_values, color=colors, edgecolor='black')
ax.axvspan(7, 10, alpha=0.1, color='green', label='PASS range')
ax.set_xlabel('log₁₀(kon / M⁻¹s⁻¹)', fontsize=14)
ax.set_title(f'EXP-13: Association Rate — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-13_kon_comparison.png", dpi=300)
plt.close(fig)
```

### Figure 3: Electrostatic Steering Summary

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
# Left: kon bar chart
ax1.bar(['WT', 'Neutral'], [np.log10(kon_avg), np.log10(k_on_ctrl)],
        color=['navy', 'red'], edgecolor='black')
ax1.set_ylabel('log₁₀(kon)', fontsize=12)
ax1.set_title('Electrostatic Steering Effect', fontsize=13)
ax1.annotate(f'Δ = {steering_ratio:.1f} orders', xy=(0.5, max(np.log10(kon_avg), np.log10(k_on_ctrl))*0.9),
             fontsize=14, ha='center', fontweight='bold')

# Right: classification summary
summary_text = (f"Experiment: kon = 3.7×10⁸ M⁻¹s⁻¹\n"
                f"Predicted: kon = {kon_avg:.2e} M⁻¹s⁻¹\n"
                f"Steering: {steering_ratio:.1f} orders\n"
                f"Classification: {classification}")
ax2.text(0.5, 0.5, summary_text, transform=ax2.transAxes,
         fontsize=14, verticalalignment='center', horizontalalignment='center',
         bbox=dict(boxstyle='round', facecolor='lightblue' if classification=='PASS' else 'lightyellow'))
ax2.axis('off')
ax2.set_title('Classification Summary', fontsize=13)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-13_steering_summary.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-13: Barnase–Barstar Association Kinetics — Results Report

**Experiment ID:** EXP-13  **Feature ID:** F-13  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Abstract
Predicted kon for barnase–barstar association using PMF-based diffusion theory.

## Results
| Method | kon (M⁻¹s⁻¹) | log₁₀(kon) |
|--------|---------------|-------------|
| Experiment | 3.7×10⁸ | 8.57 |
| SSN | [value] | [value] |
| Smoluchowski | [value] | [value] |
| Geometric mean | [value] | [value] |
| Neutral control | [value] | [value] |

Electrostatic steering: [X] orders of magnitude

## Figures
1. PMF overlay (WT vs neutralized)
2. kon comparison bars
3. Steering summary

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

EXP_ID = "EXP-13"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures", "outputs/structures",
               "outputs/equilibration", "outputs/production",
               "outputs/umbrella", "outputs/umbrella_neutral"]:
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

# Update all Simulation() calls to include platform and properties.
# Both WT and charge-neutralized control simulations use CUDA.
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

# ─── Checkpoint points for EXP-13 ───
# After minimization:         save_checkpoint(sim, "minimize")
# After NVT/NPT:              save_checkpoint(sim, "nvt"), save_checkpoint(sim, "npt")
# Production (every 10ns):    save_checkpoint(sim, f"production_{ns}ns")
# After each US window (WT):  save_checkpoint(sim, f"us_wt_window_{i:03d}")
# After each US window (neutral): save_checkpoint(sim, f"us_neutral_window_{i:03d}")
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

# ─── EXP-13 runtime estimates (A100 40GB) ───
# Production (50 ns):  ~1–2 hours
# US WT (51 windows × 10 ns): ~10–20 hours
# US neutral control:  ~10–20 hours
# Total: may require 2 Colab sessions — checkpoint between WT and neutral
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

# ─── EXP-13–specific recovery ───
# 1. Electrostatic steering analysis failure with neutral control:
#    - Verify charge neutralization was applied correctly
#    - Check that PMF difference (WT vs neutral) is physically reasonable
# 2. SSN/Smoluchowski calculation divergence:
#    - Ensure PMF is smooth — apply Savitzky-Golay filter if noisy
# 3. Session split: run WT US in session 1, neutral US in session 2
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
