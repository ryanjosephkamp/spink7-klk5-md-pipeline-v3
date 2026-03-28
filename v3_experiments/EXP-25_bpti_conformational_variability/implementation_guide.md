# EXP-25: BPTI Conformational Variability — Implementation Guide

**Experiment ID:** EXP-25  
**Feature ID:** F-25 (benchmarks.md)  
**Category:** Dynamic (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies BPTI conformational variability from MD: the time-averaged Cα RMSD fluctuation around the mean structure. Experimental reference: 0.40 Å (NMR ensemble, Berndt et al. 1992), 95% CI [0.25, 0.55] Å. Tests whether force field produces realistic thermal dynamics.

### 2. Hypotheses

- H₁: Cα RMSD fluctuation = 0.40 ± 0.15 Å (CI [0.25, 0.55])
- H₂: Core RMSF < 0.3 Å, loop RMSF > 0.5 Å

### 3. Classification (§25.1)

- PASS: [0.25, 0.55] Å
- MARGINAL: [0.15, 0.70] Å
- FAIL: Outside marginal range

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-25_bpti_conformational_variability"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory, align_trajectory
from src.analyze.structural import compute_rmsd, compute_rmsf
from src.analyze.statistics import block_average
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load EXP-23 Trajectory

```python
EXP23_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-23_bpti_structure" / "outputs"
traj_path = EXP23_DIR / "production" / "production_trajectory.dcd"
top_path = EXP23_DIR / "solvated_bpti.pdb"

traj = load_trajectory(str(traj_path), str(top_path), stride=1)
protein_atoms = traj.topology.select("protein")
traj_prot = traj.atom_slice(protein_atoms)
print(f"Loaded: {traj_prot.n_frames} frames")
```

### Step 3: Compute Mean Structure and RMSD Fluctuation

```python
# Align all frames to first frame
ca_idx = traj_prot.topology.select("name CA")
bb_idx = traj_prot.topology.select("backbone")
traj_aligned = traj_prot.superpose(traj_prot[0], atom_indices=bb_idx)

# Compute mean Cα coordinates
ca_xyz = traj_aligned.xyz[:, ca_idx, :]  # (n_frames, n_ca, 3)
mean_xyz = np.mean(ca_xyz, axis=0)  # (n_ca, 3)

# RMSD of each frame to mean structure
rmsd_to_mean = np.sqrt(np.mean(np.sum((ca_xyz - mean_xyz)**2, axis=2), axis=1))
rmsd_to_mean_ang = rmsd_to_mean * 10  # nm → Å

conformational_variability = np.mean(rmsd_to_mean_ang)
cv_std = np.std(rmsd_to_mean_ang)

print(f"Conformational variability: {conformational_variability:.3f} ± {cv_std:.3f} Å")
print(f"Expected: 0.40 ± 0.15 Å, CI [0.25, 0.55]")
```

### Step 4: Per-Residue RMSF

```python
# Cα RMSF
rmsf_ca = np.sqrt(np.mean(np.sum((ca_xyz - mean_xyz)**2, axis=2), axis=0)) * 10  # Å
residue_ids = [traj_prot.topology.atom(idx).residue.resSeq for idx in ca_idx]

# Classify regions
core_mask = rmsf_ca < 0.3
loop_mask = rmsf_ca > 0.5
n_core = np.sum(core_mask)
n_loop = np.sum(loop_mask)
print(f"Core residues (RMSF<0.3): {n_core}, Loop residues (RMSF>0.5): {n_loop}")
print(f"RMSF: mean={np.mean(rmsf_ca):.3f}, min={np.min(rmsf_ca):.3f}, max={np.max(rmsf_ca):.3f} Å")
```

### Step 5: Block Averaging (Convergence)

```python
n_blocks = 5
block_size = traj_prot.n_frames // n_blocks
block_cvs = []

for b in range(n_blocks):
    start = b * block_size
    end = start + block_size
    block_xyz = ca_xyz[start:end]
    block_mean = np.mean(block_xyz, axis=0)
    block_rmsd = np.sqrt(np.mean(np.sum((block_xyz - block_mean)**2, axis=2), axis=1))
    block_cvs.append(np.mean(block_rmsd) * 10)

block_cvs = np.array(block_cvs)
block_se = np.std(block_cvs) / np.sqrt(n_blocks)
block_cv_pct = np.std(block_cvs) / np.mean(block_cvs) * 100

print(f"Block averages: {block_cvs}")
print(f"Block mean: {np.mean(block_cvs):.3f} ± {block_se:.3f} Å")
print(f"CV: {block_cv_pct:.1f}% (target: <15%)")
```

### Step 6: Classification

```python
if 0.25 <= conformational_variability <= 0.55:
    classification = "PASS"
elif 0.15 <= conformational_variability <= 0.70:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-25", "feature_id": "F-25",
    "conformational_variability_ang": float(conformational_variability),
    "cv_std_ang": float(cv_std),
    "benchmark_ang": 0.40,
    "ci_lower": 0.25, "ci_upper": 0.55,
    "rmsf_mean_ang": float(np.mean(rmsf_ca)),
    "core_residues": int(n_core),
    "loop_residues": int(n_loop),
    "block_averages": block_cvs.tolist(),
    "block_cv_pct": float(block_cv_pct),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-25: CV={conformational_variability:.3f} Å → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: RMSD Fluctuation Time Series

```python
fig, ax = plt.subplots(figsize=(12, 5))
time_ns = np.arange(len(rmsd_to_mean_ang)) * 0.01
ax.plot(time_ns, rmsd_to_mean_ang, 'b-', linewidth=0.5, alpha=0.7)
ax.axhspan(0.25, 0.55, alpha=0.1, color='green', label='PASS [0.25, 0.55] Å')
ax.axhline(y=0.40, color='green', linestyle='--', alpha=0.7, label='Benchmark (0.40 Å)')
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Cα RMSD to Mean (Å)', fontsize=12)
ax.set_title(f'EXP-25: BPTI Conformational Variability — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-25_variability_timeseries.png", dpi=300)
plt.close(fig)
```

### Figure 2: Per-Residue RMSF Profile

```python
fig, ax = plt.subplots(figsize=(14, 5))
colors_rmsf = ['#2ecc71' if r < 0.3 else '#e74c3c' if r > 0.5 else '#3498db' for r in rmsf_ca]
ax.bar(residue_ids, rmsf_ca, color=colors_rmsf, edgecolor='none', width=0.8)
ax.axhline(y=0.3, color='green', linestyle='--', alpha=0.5, label='Core (<0.3 Å)')
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Loop (>0.5 Å)')
ax.set_xlabel('Residue Number', fontsize=12)
ax.set_ylabel('Cα RMSF (Å)', fontsize=12)
ax.set_title('EXP-25: Per-Residue RMSF Classification', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-25_rmsf_profile.png", dpi=300)
plt.close(fig)
```

### Figure 3: Block Convergence

```python
fig, ax = plt.subplots(figsize=(10, 6))
blocks = np.arange(1, n_blocks + 1)
ax.bar(blocks, block_cvs, color='steelblue', edgecolor='black')
ax.axhline(y=conformational_variability, color='red', linestyle='-', linewidth=2,
            label=f'Overall: {conformational_variability:.3f} Å')
ax.axhspan(0.25, 0.55, alpha=0.1, color='green')
ax.set_xlabel('Block (20 ns each)', fontsize=12)
ax.set_ylabel('Conformational Variability (Å)', fontsize=12)
ax.set_title(f'EXP-25: Block Convergence (CV={block_cv_pct:.1f}%)', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-25_block_convergence.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-25: BPTI Conformational Variability — Results Report

**Experiment ID:** EXP-25  **Feature ID:** F-25  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| Cα RMSD fluctuation | [val] Å | [0.25, 0.55] | [P/M/F] |
| Benchmark | 0.40 Å | Berndt 1992 | — |
| Core residues (RMSF<0.3) | [n] | Present | — |
| Loop residues (RMSF>0.5) | [n] | Present | — |
| Block CV | [val]% | <15% | [P/F] |

## Figures
1. Variability time series
2. Per-residue RMSF profile
3. Block convergence

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-25 loads the EXP-23 apo BPTI trajectory for conformational variability analysis.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-25"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

EXP23_DRIVE = Path("/content/drive/MyDrive/v3_gpu_results/EXP-23")
assert (EXP23_DRIVE / "outputs" / "production").exists(), \
    "EXP-23 production outputs required — run EXP-23 first"

PROJECT_ROOT = Path("/content")
EXP_DIR = DRIVE_BASE
OUTPUT_DIR = DRIVE_BASE / "outputs"
FIGURES_DIR = DRIVE_BASE / "figures"

print(f"Colab environment ready. EXP-23 data: {EXP23_DRIVE}")
```

### §5.2 GPU Platform Selection

```python
# §5.2 GPU Platform Selection
import openmm
from src.simulate.platform import select_platform

platform = select_platform("CUDA")
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
print(f"Platform: {platform.getName()}")

# EXP-25 is trajectory analysis (RMSD, RMSF, block averaging).
# GPU accelerates alignment and structural computations.
```

### §5.3 Checkpoint and Resume Integration

```python
# §5.3 Checkpoint and Resume Integration
import json, time
import numpy as np

CHECKPOINT_DIR = DRIVE_BASE / "checkpoints"

def save_analysis_checkpoint(data, phase_name):
    """Save analysis progress to Drive."""
    ckpt = CHECKPOINT_DIR / f"{phase_name}.npz"
    np.savez(ckpt, **{k: v for k, v in data.items() if isinstance(v, np.ndarray)})
    meta = {k: v for k, v in data.items() if not isinstance(v, np.ndarray)}
    meta["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(CHECKPOINT_DIR / f"{phase_name}_meta.json", 'w') as f:
        json.dump(meta, f, indent=2, default=str)
    print(f"  Analysis checkpoint saved: {phase_name}")

# ─── Checkpoint points for EXP-25 ───
# After mean structure: save_analysis_checkpoint(mean_data, "mean_structure")
# After RMSD fluctuation: save_analysis_checkpoint(rmsd_data, "rmsd_fluct")
# After per-residue RMSF: save_analysis_checkpoint(rmsf_data, "rmsf")
# After block averaging: save_analysis_checkpoint(block_data, "block_avg")
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

# ─── EXP-25 runtime estimates (A100 40GB) ───
# Trajectory loading + alignment: ~10 minutes
# RMSD/RMSF analysis: ~10–20 minutes
# Block averaging (convergence): ~10 minutes
# Total: < 1 hour
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
# ─── EXP-25–specific recovery ───
# 1. EXP-23 trajectory not found: verify EXP-23 completed on Drive
# 2. Block convergence CV > 15%:
#    - Not a GPU error — indicates insufficient sampling in EXP-23
#    - Document as limitation; may need longer production run
# 3. RMSD fluctuation outside expected range:
#    - Check trajectory alignment quality
#    - Verify correct atom selection ("name CA" for Cα RMSD)
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
