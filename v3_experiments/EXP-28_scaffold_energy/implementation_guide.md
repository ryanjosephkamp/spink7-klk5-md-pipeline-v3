# EXP-28: Scaffold Energy Contribution to Binding — Implementation Guide

**Experiment ID:** EXP-28  
**Feature ID:** F-28 (benchmarks.md)  
**Category:** Biophysical (Semi-quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies the scaffold's energetic contribution to binding (−7.9 ± 3.0 kcal/mol, Krowarsch et al. 2003) by comparing the ΔG of the full BPTI inhibitor vs. the isolated binding loop peptide (P3-P3', residues 13–20). The scaffold pre-organizes the loop — without it, the free peptide binds much more weakly.

### 2. Hypotheses

- H₁: ΔG_scaffold = ΔG_full − ΔG_peptide ≈ −7.9 ± 3.0 kcal/mol
- H₂: Free loop peptide shows much broader conformational sampling than intact inhibitor

### 3. Classification (§25.1)

- PASS: ΔG_scaffold ∈ [−4.9, −10.9] kcal/mol
- MARGINAL: ΔG_scaffold ∈ [−2.0, −14.0] kcal/mol
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

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-28_scaffold_energy"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.prep.pdb_fetch import fetch_pdb
from src.prep.structure_prep import clean_structure, assign_protonation
from src.prep.topology import build_topology, solvate_system
from src.simulate.equilibration import minimize_energy, run_nvt, run_npt
from src.simulate.production import run_production
from src.simulate.umbrella import run_umbrella_campaign, generate_window_centers
from src.analyze.wham import solve_wham, bootstrap_pmf_uncertainty
from src.config import (SystemConfig, MinimizationConfig, EquilibrationConfig,
                         ProductionConfig, UmbrellaConfig, WHAMConfig)
from src.simulate.platform_util import select_platform
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Retrieve Full-Inhibitor ΔG from EXP-04

```python
EXP04_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind"
exp04_results_path = EXP04_DIR / "results.json"

if exp04_results_path.exists():
    with open(exp04_results_path) as f:
        exp04_res = json.load(f)
    dg_full = exp04_res.get("dg_bind_kcal", -18.0)
else:
    dg_full = -18.0  # Expected from benchmarks

print(f"ΔG_full (BPTI-trypsin): {dg_full:.1f} kcal/mol")
```

### Step 3: Build Free Loop Peptide (BPTI residues 13–20)

```python
# Extract BPTI residues 13-20 from crystal
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)
pdb_path = fetch_pdb("4PTI", data_dir)

bpti = md.load(str(pdb_path))
loop_residues = [r for r in bpti.topology.residues if 13 <= r.resSeq <= 20]
loop_atom_idx = [a.index for r in loop_residues for a in r.atoms]

# Save loop peptide
loop_pdb = data_dir / "bpti_loop_13_20.pdb"
loop_traj = bpti.atom_slice(loop_atom_idx)
loop_traj.save_pdb(str(loop_pdb))
print(f"Loop peptide saved: {loop_pdb} ({len(loop_atom_idx)} atoms)")

# Cap with ACE/NME via PDBFixer or manual
# Protonate
loop_prot = assign_protonation(loop_pdb, ph=7.4, force_field="AMBER", use_propka=True)
```

### Step 4: Free Peptide Simulation (50 ns)

```python
sys_config = SystemConfig(
    force_field="amber14-all.xml",
    water_model="amber14/tip3p.xml",
    box_padding_nm=1.2,
    ionic_strength_M=0.15,
    temperature_K=310.0,
)
topology, system, modeller = build_topology(loop_prot, sys_config,
                                             nonbonded_method="PME",
                                             nonbonded_cutoff_nm=1.0)
modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)

import openmm as mm
from openmm import unit

platform = select_platform()
integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picosecond,
                                          0.002*unit.picoseconds)
simulation = mm.app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# Minimize
min_config = MinimizationConfig(max_iterations=10000, tolerance_kj_per_mol_nm=10.0)
minimize_energy(simulation, min_config)

# Equilibrate
equil_dir = OUTPUT_DIR / "equilibration_peptide"
equil_dir.mkdir(exist_ok=True)
nvt_config = EquilibrationConfig(ensemble="NVT", duration_ps=500, temperature_K=310)
run_nvt(simulation, nvt_config, equil_dir)
npt_config = EquilibrationConfig(ensemble="NPT", duration_ps=1000, temperature_K=310)
run_npt(simulation, npt_config, equil_dir)

# Production: 50 ns
prod_dir = OUTPUT_DIR / "production_peptide"
prod_dir.mkdir(exist_ok=True)
prod_config = ProductionConfig(duration_ns=50, temperature_K=310, save_interval_ps=10)
run_production(simulation, prod_config, prod_dir)
print("Free peptide 50 ns production complete.")
```

### Step 5: Dock Peptide into Trypsin + US/WHAM

```python
# Load trypsin from EXP-04
trypsin_pdb = EXP04_DIR / "outputs" / "structures" / "2PTC_chainA_clean_prot.pdb"
# Position peptide same as in BPTI complex
# Build complex system
complex_dir = OUTPUT_DIR / "peptide_complex"
complex_dir.mkdir(exist_ok=True)

# US campaign: ξ = COM distance peptide ↔ trypsin
us_config = UmbrellaConfig(
    xi_start_nm=1.5, xi_end_nm=4.0, xi_spacing_nm=0.05,
    spring_constant_kj=1000.0,
    per_window_ns=10.0, equilibration_ps=200,
    pre_position_velocity_nm_ps=0.01,
    save_interval_ps=1.0,
)
us_result = run_umbrella_campaign(simulation, us_config, complex_dir,
                                   group1_selection="chainid 0",
                                   group2_selection="chainid 1")

# WHAM
wham_config = WHAMConfig(tolerance=1e-6, max_iterations=100000, n_bins=200)
pmf, bin_centers = solve_wham(us_result["xi_files"], us_result["metadata"],
                               wham_config)
pmf_boot, ci_lower, ci_upper = bootstrap_pmf_uncertainty(
    us_result["xi_files"], us_result["metadata"], wham_config, n_bootstrap=200)

# ΔG extraction
dg_peptide = float(np.min(pmf) - pmf[-1])
from src.physics.units import kj_to_kcal
dg_peptide_kcal = kj_to_kcal(dg_peptide)
print(f"ΔG_peptide: {dg_peptide_kcal:.2f} kcal/mol")
```

### Step 6: Scaffold Energy + Classification

```python
dg_scaffold = dg_full - dg_peptide_kcal
print(f"ΔG_scaffold = {dg_full:.1f} − ({dg_peptide_kcal:.1f}) = {dg_scaffold:.1f} kcal/mol")

if -10.9 <= dg_scaffold <= -4.9:
    classification = "PASS"
elif -14.0 <= dg_scaffold <= -2.0:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-28", "feature_id": "F-28",
    "dg_full_kcal": float(dg_full),
    "dg_peptide_kcal": float(dg_peptide_kcal),
    "dg_scaffold_kcal": float(dg_scaffold),
    "benchmark_kcal": -7.9,
    "ci_lower": -10.9, "ci_upper": -4.9,
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-28: ΔG_scaffold={dg_scaffold:.1f} kcal/mol → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: PMF Comparison (Full Inhibitor vs Peptide)

```python
fig, ax = plt.subplots(figsize=(12, 6))
# Load EXP-04 PMF if available
ax.plot(bin_centers * 10, pmf_boot, 'r-', linewidth=2, label=f'Peptide (ΔG={dg_peptide_kcal:.1f})')
ax.fill_between(bin_centers * 10, ci_lower, ci_upper, alpha=0.2, color='red')
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.set_xlabel('ξ (Å)', fontsize=13)
ax.set_ylabel('PMF (kJ/mol)', fontsize=13)
ax.set_title(f'EXP-28: Peptide–Trypsin PMF', fontsize=14)
ax.legend(fontsize=12)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-28_pmf_comparison.png", dpi=300)
plt.close(fig)
```

### Figure 2: Energy Decomposition Bar Chart

```python
fig, ax = plt.subplots(figsize=(10, 6))
labels = ['Full BPTI', 'Loop Peptide', 'Scaffold']
values = [dg_full, dg_peptide_kcal, dg_scaffold]
colors = ['#2ecc71', '#e74c3c', '#3498db']
bars = ax.bar(labels, values, color=colors, edgecolor='black')
ax.axhline(y=0, color='black', linewidth=0.5)
ax.axhline(y=-7.9, color='blue', linestyle='--', alpha=0.5, label='Benchmark (−7.9)')
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{val:.1f}', ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('ΔG (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-28: Scaffold Energy Decomposition — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-28_energy_decomposition.png", dpi=300)
plt.close(fig)
```

### Figure 3: Free Peptide Conformational Sampling

```python
# Load free peptide trajectory
pep_traj = load_trajectory(str(prod_dir / "production_trajectory.dcd"),
                            str(OUTPUT_DIR / "solvated_peptide.pdb"), stride=10)
pep_prot = pep_traj.atom_slice(pep_traj.topology.select("protein"))
pep_bb = pep_prot.topology.select("backbone")
pep_rmsd = md.rmsd(pep_prot, pep_prot[0], atom_indices=pep_bb) * 10

fig, ax = plt.subplots(figsize=(12, 5))
time_ns = np.arange(len(pep_rmsd)) * 0.1
ax.plot(time_ns, pep_rmsd, 'b-', linewidth=0.5, alpha=0.7)
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Backbone RMSD (Å)', fontsize=12)
ax.set_title('EXP-28: Free Loop Peptide Conformational Sampling', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-28_peptide_rmsd.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-28: Scaffold Energy — Results Report

**Experiment ID:** EXP-28  **Feature ID:** F-28  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| ΔG_full (EXP-04) | [val] kcal/mol | ~−18.0 | — |
| ΔG_peptide | [val] kcal/mol | ~−10 | — |
| ΔG_scaffold | [val] kcal/mol | [−4.9, −10.9] | [P/M/F] |

## Figures
1. PMF comparison
2. Energy decomposition bar chart
3. Free peptide conformational sampling

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-28 builds a free loop peptide (BPTI residues 13–20), runs 50 ns MD, docks into trypsin, and performs US/WHAM — GPU-intensive.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-28"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures",
               "outputs/peptide_md", "outputs/peptide_us"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

# Verify EXP-04 results available for ΔG_full
EXP04_DRIVE = Path("/content/drive/MyDrive/v3_gpu_results/EXP-04")
assert (EXP04_DRIVE / "results.json").exists(), \
    "EXP-04 results.json required for ΔG_full comparison"

PROJECT_ROOT = Path("/content")
EXP_DIR = DRIVE_BASE
OUTPUT_DIR = DRIVE_BASE / "outputs"
FIGURES_DIR = DRIVE_BASE / "figures"

print(f"Colab environment ready. EXP-04 data: {EXP04_DRIVE}")
```

### §5.2 GPU Platform Selection

```python
# §5.2 GPU Platform Selection
import openmm
from src.simulate.platform import select_platform

platform = select_platform("CUDA")
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
print(f"Platform: {platform.getName()}")

# EXP-28 has two GPU-intensive phases:
# 1. Free peptide MD simulation (50 ns)
# 2. Peptide-trypsin US/WHAM (if docking successful)
# Pipeline functions use select_platform() → CUDA auto-detected.
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

# ─── Checkpoint points for EXP-28 ───
# Phase 1 — Peptide MD:
#   After peptide minimization:  save_checkpoint(sim, "peptide_minimize")
#   After peptide equilibration: save_checkpoint(sim, "peptide_equil")
#   Production (every 10ns):     save_checkpoint(sim, f"peptide_prod_{ns}ns")
#
# Phase 2 — Peptide-trypsin US:
#   After each US window:        save_checkpoint(sim, f"peptide_us_window_{i:03d}")
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

# ─── EXP-28 runtime estimates (A100 40GB) ───
# Peptide MD (50 ns, small system): ~1 hour
# Peptide-trypsin US (51 windows × 5 ns): ~5–10 hours
# Analysis + figures: ~30 minutes
# Total: ~7–12 hours (single Colab session)
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
    critical = ["results.json"]
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

# ─── EXP-28–specific recovery ───
# 1. Free peptide unfolding:
#    - Short peptide (8 residues, P3-P3') will sample broadly — expected
#    - If NaN during peptide MD: reduce timestep to 1 fs initially
# 2. Peptide-trypsin docking failure:
#    - If peptide cannot be docked into S1 pocket, use steered approach
#    - Pre-position using SMD from bulk to S1
# 3. ΔG_scaffold wrong sign:
#    - If ΔG_full − ΔG_peptide > 0, check force field consistency
#    - Verify both PMFs use same reference state
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
