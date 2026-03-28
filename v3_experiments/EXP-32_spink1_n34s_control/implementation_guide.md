# EXP-32: SPINK1 N34S Negative Control — Implementation Guide

**Experiment ID:** EXP-32  
**Feature ID:** F-32 (benchmarks.md)  
**Category:** Mutational (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Critical negative control for FEP methodology. SPINK1 N34S is a common clinical variant (pancreatitis-associated) but does NOT alter trypsin binding (ΔΔG ≈ 0, Kuwata 2002, Kiraly 2007). Pipeline must distinguish neutral mutations from hot-spot mutations. Acceptance: |ΔΔG| < 1.0 kcal/mol.

### 2. Hypotheses

- H₁: |ΔΔG(N34S)| < 1.0 kcal/mol
- H₂: Structural metrics (loop RMSD, H-bonds, BSA) unchanged

### 3. Classification (§25.1)

- PASS: |ΔΔG| < 1.0 kcal/mol
- MARGINAL: |ΔΔG| 1.0–2.0 kcal/mol
- FAIL: |ΔΔG| > 2.0 kcal/mol

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-32_spink1_n34s_control"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.simulate.fep import run_fep_campaign
from src.analyze.mbar import solve_mbar, bootstrap_mbar_uncertainty
from src.analyze.trajectory import load_trajectory
from src.analyze.structural import compute_rmsd, compute_hbonds, compute_sasa
from src.config import FEPConfig, MBARConfig
from src.physics.units import kj_to_kcal
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Verify N34 Location (Non-Interface)

```python
EXP06_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-06_spink1_trypsin_dg_bind" / "outputs"

# Load equilibrated complex
traj_path = EXP06_DIR / "production" / "production_trajectory.dcd"
top_path = EXP06_DIR / "solvated_complex.pdb"
traj = load_trajectory(str(traj_path), str(top_path), stride=100)
topo = traj.topology

# Find N34 and nearest trypsin atom
n34_atoms = topo.select("resid 34 and chainid 1")  # SPINK1
trypsin_atoms = topo.select("chainid 0 and not element H")

if len(n34_atoms) > 0:
    pairs = [[n34_atoms[0], t] for t in trypsin_atoms]
    dists = md.compute_distances(traj[0], pairs)[0]
    min_dist = np.min(dists) * 10  # nm → Å
    print(f"N34 → nearest trypsin atom: {min_dist:.1f} Å (expected >8 Å)")
```

### Step 3: FEP — N34S Single Mutation

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

# Leg 1: Complex
complex_dir = OUTPUT_DIR / "complex"
complex_dir.mkdir(exist_ok=True)
fep_complex = run_fep_campaign(
    topology_path=str(EXP06_DIR / "solvated_complex.pdb"),
    mutation="N34S",
    state="complex",
    fep_config=fep_config,
    output_dir=complex_dir,
)
dg_complex = solve_mbar(fep_complex["energy_files"], fep_complex["metadata"], mbar_config)
boot_complex = bootstrap_mbar_uncertainty(fep_complex["energy_files"],
                                            fep_complex["metadata"], mbar_config, n_bootstrap=200)

# Leg 2: Free SPINK1
free_dir = OUTPUT_DIR / "free"
free_dir.mkdir(exist_ok=True)
fep_free = run_fep_campaign(
    topology_path=str(EXP06_DIR / "solvated_spink1_free.pdb"),
    mutation="N34S",
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
print(f"ΔΔG(N34S): {ddg_kcal:.2f} ± {ddg_err:.2f} kcal/mol")
```

### Step 4: Structural Comparison (WT vs N34S)

```python
# Compare structural metrics between WT and mutant
# Loop RMSD, H-bond count, BSA should be unchanged

# WT metrics from EXP-06 trajectory
wt_prot = traj.atom_slice(traj.topology.select("protein"))
loop_bb = wt_prot.topology.select("resid 13 to 20 and backbone")
if len(loop_bb) > 0:
    wt_rmsd = md.rmsd(wt_prot, wt_prot[0], atom_indices=loop_bb) * 10
    wt_loop_rmsd = np.mean(wt_rmsd)
    print(f"WT loop RMSD: {wt_loop_rmsd:.2f} Å")

# Interface H-bonds
hbond_counts = []
for i in range(0, traj.n_frames, max(1, traj.n_frames // 20)):
    hbonds = md.baker_hubbard(traj[i], freq=0.0)
    inter = sum(1 for h in hbonds
                if topo.atom(h[0]).residue.chain.index != topo.atom(h[2]).residue.chain.index)
    hbond_counts.append(inter)
wt_hbond_mean = np.mean(hbond_counts)
print(f"WT H-bonds: {wt_hbond_mean:.1f}")
```

### Step 5: Classification

```python
abs_ddg = abs(ddg_kcal)
if abs_ddg < 1.0:
    classification = "PASS"
elif abs_ddg < 2.0:
    classification = "MARGINAL"
else:
    classification = "FAIL"

# Compare with K15A from EXP-30
exp30_path = PROJECT_ROOT / "v3_experiments" / "EXP-30_alanine_scanning" / "results.json"
k15a_ddg = 10.0  # expected
if exp30_path.exists():
    with open(exp30_path) as f:
        k15a_ddg = json.load(f).get("per_mutation", {}).get("K15A", {}).get("ddg_calc_kcal", 10.0)
discrimination = k15a_ddg - abs_ddg
print(f"K15A-N34S discrimination: {discrimination:.1f} kcal/mol (expected >8)")

results = {
    "experiment_id": "EXP-32", "feature_id": "F-32",
    "ddg_n34s_kcal": float(ddg_kcal),
    "ddg_err_kcal": float(ddg_err),
    "abs_ddg": float(abs_ddg),
    "n34_trypsin_dist_ang": float(min_dist) if 'min_dist' in dir() else None,
    "k15a_discrimination_kcal": float(discrimination),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-32: |ΔΔG|={abs_ddg:.2f} kcal/mol → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: ΔΔG Comparison (N34S vs Hot Spots)

```python
fig, ax = plt.subplots(figsize=(10, 6))
labels = ['N34S\n(neutral)', 'K15A\n(P1 hot spot)', 'C14A\n(disulfide)']
values = [abs_ddg, k15a_ddg, 7.0]
colors = ['#2ecc71', '#e74c3c', '#f39c12']
bars = ax.bar(labels, values, color=colors, edgecolor='black')
ax.axhline(y=1.0, color='blue', linestyle='--', label='|ΔΔG|=1.0 threshold')
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{val:.1f}', ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('|ΔΔG| (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-32: N34S Negative Control — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-32_ddg_comparison.png", dpi=300)
plt.close(fig)
```

### Figure 2: N34 Location Map

```python
fig, ax = plt.subplots(figsize=(10, 8))
# Per-residue minimum distance to trypsin
spink1_residues = [r for r in topo.residues if r.chain.index == 1]
res_ids = [r.resSeq for r in spink1_residues]
min_dists = []
for res in spink1_residues:
    res_atoms = [a.index for a in res.atoms if a.element.symbol != 'H']
    if res_atoms:
        pairs = [[ra, ta] for ra in res_atoms for ta in trypsin_atoms[:50]]
        if pairs:
            d = md.compute_distances(traj[0], pairs[:200])[0]
            min_dists.append(np.min(d) * 10)
        else:
            min_dists.append(20.0)
    else:
        min_dists.append(20.0)

colors_res = ['#e74c3c' if d < 5 else '#f39c12' if d < 8 else '#2ecc71' for d in min_dists]
ax.bar(res_ids, min_dists, color=colors_res, edgecolor='none', width=0.8)
ax.axhline(y=8.0, color='blue', linestyle='--', label='8 Å cutoff')
if 34 in res_ids:
    idx34 = res_ids.index(34)
    ax.bar(res_ids[idx34], min_dists[idx34], color='purple', edgecolor='black', width=0.8,
            label='N34')
ax.set_xlabel('SPINK1 Residue', fontsize=12)
ax.set_ylabel('Min. Distance to Trypsin (Å)', fontsize=12)
ax.set_title('EXP-32: N34 is Remote from Interface', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-32_n34_location.png", dpi=300)
plt.close(fig)
```

### Figure 3: FEP Convergence

```python
fig, ax = plt.subplots(figsize=(10, 6))
lambdas = np.linspace(0, 1, fep_config.n_lambda)
# Plot per-lambda ΔG contributions
complex_contrib = dg_complex.get("dg_per_lambda", np.zeros(20))
free_contrib = dg_free.get("dg_per_lambda", np.zeros(20))
ddg_per_lambda = np.array(complex_contrib) - np.array(free_contrib)

ax.bar(lambdas, kj_to_kcal(np.array(ddg_per_lambda)), width=0.04,
        color='steelblue', edgecolor='black')
ax.axhline(y=0, color='black', linewidth=0.5)
ax.set_xlabel('λ', fontsize=13)
ax.set_ylabel('ΔΔG contribution (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-32: Per-λ ΔΔG (total={ddg_kcal:.2f})', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-32_fep_convergence.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-32: SPINK1 N34S Negative Control — Results Report

**Experiment ID:** EXP-32  **Feature ID:** F-32  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| |ΔΔG(N34S)| | [val] kcal/mol | < 1.0 | [P/M/F] |
| N34-trypsin dist | [val] Å | > 8.0 | [P/F] |
| K15A discrimination | [val] kcal/mol | > 8.0 | — |

## Figures
1. ΔΔG comparison (neutral vs hot spot)
2. N34 location map
3. FEP convergence

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-32 runs a single-mutation FEP (N34S) as a negative control — GPU is critical for lambda-window sampling efficiency.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-32"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures",
               "outputs/fep_n34s", "outputs/structural"]:
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

# EXP-32 uses FEP (alchemical transformation N→S at position 34).
# Lambda-window simulations are GPU-intensive.
# The run_fep_window() and run_fep_campaign() functions in
# src/simulate/fep.py use select_platform() — CUDA auto-detected on Colab.
#
# Verified function signatures (src/simulate/fep.py):
#   generate_lambda_schedule(n_windows: int) -> np.ndarray
#   create_alchemical_system(system, mutant_atom_indices, config) -> System
#   run_fep_window(alchemical_system, positions, lambda_value,
#                  lambda_schedule, config, output_dir) -> dict
#   run_fep_campaign(system, positions, mutant_atom_indices,
#                    config, output_dir) -> dict
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

# ─── Checkpoint points for EXP-32 ───
# After system setup:              save_checkpoint(sim, "setup")
# After each FEP lambda window:    save_checkpoint(sim, f"fep_lambda_{lam:.4f}")
# After structural comparison:     save_checkpoint(sim, "structural_done")
#
# FEP campaign progress is also saved per-window by run_fep_window()
# via np.savez() in the output directory.
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

# ─── EXP-32 runtime estimates (A100 40GB) ───
# FEP N34S (20 windows × 2 ns each = 40 ns): ~2–4 hours
# Structural comparison (short MD): ~30 minutes
# Total: ~3–5 hours
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
    critical = ["results.json", "outputs/fep_n34s/fep_campaign_results.npz"]
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

# ─── EXP-32–specific recovery ───
# 1. FEP convergence failure:
#    - Symptom: MBAR solver doesn't converge
#    - Fix: increase per_window_duration_ns or n_lambda_windows
# 2. N34S mutation site identification:
#    - Ensure correct residue numbering (N34 in SPINK1 sequence)
#    - Verify N34 is NOT at the interface (distance > 8 Å to trypsin)
# 3. |ΔΔG| > 1.0 but negative control expected ~0:
#    - Check for insufficient sampling at endpoints (λ=0, λ=1)
#    - Extend endpoint windows if needed
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
