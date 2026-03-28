# EXP-24: BPTI Hydrogen/Deuterium Exchange Protection Pattern — Implementation Guide

**Experiment ID:** EXP-24  
**Feature ID:** F-24 (benchmarks.md)  
**Category:** Dynamic (Semi-quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates the pipeline's ability to reproduce the H/D exchange protection pattern of BPTI. Experimental NMR identifies 11 backbone amide protons highly protected from exchange (PF > 10⁴), corresponding to buried residues in stable secondary structures or shielded by disulfide bonds. MD-derived proxies (low RMSF, high H-bond persistence, low solvent accessibility) must correctly identify ≥8 of 11 protected amides.

### 2. Hypotheses

- H₁: MD analysis correctly identifies ≥8/11 protected amides (SASA < 10% AND H-bond occupancy > 80%)
- H₂: Correlation between MD protection proxies and experimental PF yields r > 0.5

### 3. Classification (§25.1)

- PASS: ≥8/11 correctly identified; ≤3 false positives
- MARGINAL: 6–7/11 correctly identified; ≤5 false positives
- FAIL: <6/11 or >5 false positives

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-24_bpti_hd_exchange"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory, align_trajectory
from src.analyze.structural import compute_rmsd, compute_rmsf, compute_sasa, compute_hbonds
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load EXP-23 Production Trajectory

```python
EXP23_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-23_bpti_structure" / "outputs"
traj_path = EXP23_DIR / "production" / "production_trajectory.dcd"
top_path = EXP23_DIR / "solvated_bpti.pdb"

traj = load_trajectory(str(traj_path), str(top_path), stride=1)
topo = traj.topology
protein_atoms = topo.select("protein")
traj_prot = traj.atom_slice(protein_atoms)
print(f"Loaded: {traj_prot.n_frames} frames, {traj_prot.n_atoms} protein atoms")

# Define known protected amides from experimental H/D exchange (Dempsey 2001)
# Residue indices (0-based) for the 11 protected amide NHs
protected_residues_1based = [21, 22, 23, 24, 25, 14, 49, 50, 51, 52, 53]
protected_set = set(protected_residues_1based)
print(f"Protected amides (experimental): {sorted(protected_residues_1based)}")
```

### Step 3: Per-Residue RMSF

```python
# Cα RMSF
crystal_path = EXP23_DIR / "structures" / "4PTI_clean_prot.pdb"
crystal = md.load(str(crystal_path))

ca_idx = traj_prot.topology.select("name CA")
traj_aligned = traj_prot.superpose(traj_prot[0], atom_indices=traj_prot.topology.select("backbone"))
rmsf_ca = md.rmsf(traj_aligned, traj_aligned, 0, atom_indices=ca_idx) * 10  # Å

residue_ids = [traj_prot.topology.atom(idx).residue.resSeq for idx in ca_idx]
print(f"RMSF computed for {len(rmsf_ca)} residues")
```

### Step 4: Backbone N-H Solvent Accessibility

```python
# Per-residue backbone NH SASA (time-averaged)
nh_sasa_avg = {}
for res in traj_prot.topology.residues:
    if res.name == 'PRO':
        continue
    n_atoms = [a.index for a in res.atoms if a.name == 'N']
    h_atoms = [a.index for a in res.atoms if a.name == 'H']
    if n_atoms:
        sasa_frames = md.shrake_rupley(traj_prot[::10], mode='atom')
        sasa_n = np.mean(sasa_frames[:, n_atoms[0]]) * 100  # nm² → Å²
        nh_sasa_avg[res.resSeq] = sasa_n

# Normalize: fraction of max possible SASA
max_sasa = max(nh_sasa_avg.values()) if nh_sasa_avg else 1.0
nh_sasa_frac = {k: v / max_sasa * 100 for k, v in nh_sasa_avg.items()}
print(f"NH SASA computed for {len(nh_sasa_frac)} residues")
```

### Step 5: Backbone H-Bond Occupancy

```python
# For each backbone N-H, compute H-bond occupancy to any backbone C=O
hbond_occ = {}
n_sample = min(traj_prot.n_frames, 1000)
sample_idx = np.linspace(0, traj_prot.n_frames-1, n_sample, dtype=int)

hbond_counts = {}
for res in traj_prot.topology.residues:
    hbond_counts[res.resSeq] = 0

for i in sample_idx:
    hbonds = md.baker_hubbard(traj_prot[i], freq=0.0)
    for d, _h, a in hbonds:
        d_res = traj_prot.topology.atom(d).residue
        a_res = traj_prot.topology.atom(a).residue
        d_name = traj_prot.topology.atom(d).name
        a_name = traj_prot.topology.atom(a).name
        # Backbone N-H → backbone C=O
        if d_name == 'N' and a_name == 'O':
            hbond_counts[d_res.resSeq] = hbond_counts.get(d_res.resSeq, 0) + 1

for res_id, count in hbond_counts.items():
    hbond_occ[res_id] = count / n_sample * 100  # %

print(f"H-bond occupancy computed for {len(hbond_occ)} residues")
```

### Step 6: Protection Factor Proxy Classification

```python
# Binary classification: protected if SASA < 10% AND H-bond occupancy > 80%
predicted_protected = set()
for res_id in nh_sasa_frac:
    sasa_low = nh_sasa_frac.get(res_id, 100) < 10
    hb_high = hbond_occ.get(res_id, 0) > 80
    if sasa_low and hb_high:
        predicted_protected.add(res_id)

# Classification metrics
true_positives = predicted_protected & protected_set
false_negatives = protected_set - predicted_protected
false_positives = predicted_protected - protected_set

n_tp = len(true_positives)
n_fn = len(false_negatives)
n_fp = len(false_positives)

print(f"True positives: {n_tp}/11")
print(f"False negatives: {n_fn}")
print(f"False positives: {n_fp}")
print(f"Protected (predicted): {sorted(predicted_protected)}")
```

### Step 7: Classification

```python
if n_tp >= 8 and n_fp <= 3:
    classification = "PASS"
elif n_tp >= 6 and n_fp <= 5:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-24", "feature_id": "F-24",
    "true_positives": n_tp,
    "false_negatives": n_fn,
    "false_positives": n_fp,
    "predicted_protected": sorted(list(predicted_protected)),
    "experimental_protected": sorted(protected_residues_1based),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-24: {n_tp}/11 TP, {n_fp} FP → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Protection Factor Proxy Profile

```python
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
residues = sorted(nh_sasa_frac.keys())

# RMSF
rmsf_vals = []
for r in residues:
    idx = [i for i, rid in enumerate(residue_ids) if rid == r]
    rmsf_vals.append(rmsf_ca[idx[0]] if idx else 0)
ax1.bar(residues, rmsf_vals, color='steelblue', width=0.8)
ax1.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)
ax1.set_ylabel('Cα RMSF (Å)')
ax1.set_title('EXP-24: BPTI H/D Exchange Protection Analysis')

# SASA
sasa_vals = [nh_sasa_frac.get(r, 0) for r in residues]
ax2.bar(residues, sasa_vals, color='orange', width=0.8)
ax2.axhline(y=10, color='red', linestyle='--', label='10% threshold')
ax2.set_ylabel('NH SASA (%)')
ax2.legend()

# H-bond occupancy
hb_vals = [hbond_occ.get(r, 0) for r in residues]
colors_hb = ['green' if r in protected_set else 'steelblue' for r in residues]
ax3.bar(residues, hb_vals, color=colors_hb, width=0.8)
ax3.axhline(y=80, color='red', linestyle='--', label='80% threshold')
ax3.set_ylabel('H-bond Occ. (%)')
ax3.set_xlabel('Residue Number')
ax3.legend()

plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-24_protection_profile.png", dpi=300)
plt.close(fig)
```

### Figure 2: Classification Summary

```python
fig, ax = plt.subplots(figsize=(10, 6))

all_residues = sorted(set(list(protected_set) + list(predicted_protected)))
categories = []
colors_cat = []
for r in all_residues:
    if r in true_positives:
        categories.append('TP')
        colors_cat.append('#2ecc71')
    elif r in false_negatives:
        categories.append('FN')
        colors_cat.append('#e74c3c')
    elif r in false_positives:
        categories.append('FP')
        colors_cat.append('#f39c12')

ax.bar(range(len(all_residues)), [1]*len(all_residues), color=colors_cat, edgecolor='black')
ax.set_xticks(range(len(all_residues)))
ax.set_xticklabels([str(r) for r in all_residues], rotation=45)
ax.set_ylabel('Classification')
ax.set_title(f'EXP-24: Protection Classification — {classification} ({n_tp}/11 TP, {n_fp} FP)')

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#2ecc71', label=f'True Positive ({n_tp})'),
                    Patch(facecolor='#e74c3c', label=f'False Negative ({n_fn})'),
                    Patch(facecolor='#f39c12', label=f'False Positive ({n_fp})')]
ax.legend(handles=legend_elements)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-24_classification.png", dpi=300)
plt.close(fig)
```

### Figure 3: RMSF Comparison (Protected vs Unprotected)

```python
fig, ax = plt.subplots(figsize=(10, 6))
prot_rmsf = [rmsf_ca[i] for i, rid in enumerate(residue_ids) if rid in protected_set]
unprot_rmsf = [rmsf_ca[i] for i, rid in enumerate(residue_ids) if rid not in protected_set]

bp = ax.boxplot([prot_rmsf, unprot_rmsf], labels=['Protected (exp)', 'Unprotected'],
                 patch_artist=True)
bp['boxes'][0].set_facecolor('#2ecc71')
bp['boxes'][1].set_facecolor('#e74c3c')
ax.set_ylabel('Cα RMSF (Å)', fontsize=12)
ax.set_title(f'EXP-24: RMSF by Protection Status — {classification}', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-24_rmsf_boxplot.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-24: BPTI H/D Exchange Protection — Results Report

**Experiment ID:** EXP-24  **Feature ID:** F-24  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| True positives | [val]/11 | ≥8 | [P/M/F] |
| False positives | [val] | ≤3 | [P/F] |
| Protected RMSF (median) | [val] Å | < 0.5 Å | — |
| Unprotected RMSF (median) | [val] Å | > 0.5 Å | — |

## Figures
1. Protection factor proxy profile (RMSF, SASA, H-bond)
2. Classification summary (TP/FN/FP)
3. RMSF boxplot by protection status

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-24 loads the EXP-23 apo BPTI trajectory for H/D exchange protection analysis.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-24"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

# EXP-24 depends on EXP-23 (apo BPTI simulation)
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

# EXP-24 is analysis-focused (RMSF, SASA, H-bond occupancy).
# GPU accelerates trajectory loading and SASA computation.
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

# ─── Checkpoint points for EXP-24 ───
# After RMSF calculation: save_analysis_checkpoint(rmsf_data, "rmsf")
# After SASA calculation: save_analysis_checkpoint(sasa_data, "sasa")
# After H-bond occupancy: save_analysis_checkpoint(hbond_data, "hbond_occ")
# After protection classification: save_analysis_checkpoint(pf_data, "protection")
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

# ─── EXP-24 runtime estimates (A100 40GB) ───
# Trajectory loading + RMSF: ~10–20 minutes
# Per-residue SASA: ~20–40 minutes
# H-bond occupancy: ~10–20 minutes
# Total: ~1–2 hours
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
# ─── EXP-24–specific recovery ───
# 1. EXP-23 trajectory not found: verify EXP-23 completed on Drive
# 2. SASA calculation crash on large trajectory:
#    - Process in chunks: traj[::stride] with stride=10
# 3. Protection factor classification edge cases:
#    - Residues near thresholds (SASA ~10%, H-bond ~80%)
#    - Document borderline cases in results report
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
