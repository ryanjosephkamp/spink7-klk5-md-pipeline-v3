# EXP-10: BPTI–Trypsin Activation Energy — Implementation Guide

**Experiment ID:** EXP-10  
**Feature ID:** F-10 (benchmarks.md)  
**Category:** Kinetic  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Measures activation energy for BPTI-trypsin dissociation from PMF barrier height along unbinding coordinate. Experimental Ea = 10.5 kcal/mol (Quast et al. 1974, Vincent & Bhatt 2007), 95% CI [4.6, 16.4] kcal/mol. PMF barrier height from EXP-04 umbrella sampling directly estimates this quantity.

### 2. Hypothesis

**H₁:** PMF barrier height ΔW‡ = 10.5 ± 5.9 kcal/mol (95% CI [4.6, 16.4]). **H₂:** Consistent with EXP-09 koff via Kramers'/TST.

### 3. Protocol

PMF from EXP-04 US/WHAM (shared data). Identify ξ_min (bound) and ξ‡ (TS). ΔW‡ = W(ξ‡) − W(ξ_min). Bootstrap uncertainty from 200 samples. Cross-validate WHAM vs MBAR. PASS: [4.6, 16.4]; MARGINAL: [2.0, 20.0]; FAIL: outside or no clear barrier.

### 4. Controls

Positive: flat PMF at large ξ. Consistency: EXP-09 cross-check. Negative: non-interacting species → flat PMF.

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path
from scipy.signal import argrelextrema
from scipy.stats import shapiro

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-10_bpti_trypsin_ea"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import UmbrellaConfig, WHAMConfig, MBARConfig, KCAL_TO_KJ
from src.simulate.umbrella import generate_window_centers
from src.analyze.wham import solve_wham, bootstrap_pmf_uncertainty
from src.analyze.mbar import solve_mbar, bootstrap_mbar_uncertainty
from src.physics.units import kj_to_kcal, kbt
import matplotlib.pyplot as plt

print("All imports successful.")
```

### Step 2: Load EXP-04 PMF and Bootstrap Samples

```python
EXP04_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind" / "outputs"
us_output = EXP04_DIR / "umbrella"

us_config = UmbrellaConfig()
window_centers = generate_window_centers(us_config)
xi_timeseries_list = [np.load(us_output / f"window_{i:03d}" / "xi_timeseries.npy")
                      for i in range(len(window_centers))]
spring_constants = np.full(len(window_centers), us_config.spring_constant_kj_mol_nm2)

# WHAM PMF
wham_config = WHAMConfig()
wham_results = solve_wham(xi_timeseries_list, window_centers, spring_constants, 310.0, wham_config)
pmf_kj = wham_results["pmf_kj_mol"]
xi_bins = wham_results["xi_bins_nm"]
pmf_kcal = pmf_kj / KCAL_TO_KJ

# Bootstrap
wham_bootstrap = bootstrap_pmf_uncertainty(xi_timeseries_list, window_centers,
                                            spring_constants, 310.0, wham_config)
pmf_std_kj = wham_bootstrap["pmf_std_kj_mol"]
pmf_std_kcal = pmf_std_kj / KCAL_TO_KJ

# MBAR PMF
mbar_config = MBARConfig()
mbar_results = solve_mbar(xi_timeseries_list, window_centers, spring_constants, 310.0, mbar_config)
pmf_mbar_kj = mbar_results["pmf_kj_mol"]
pmf_mbar_kcal = pmf_mbar_kj / KCAL_TO_KJ

print(f"PMF loaded: {len(xi_bins)} bins")
```

### Step 3: Extract Barrier Height

```python
# Bound state minimum
idx_min = np.argmin(pmf_kcal)
xi_min = xi_bins[idx_min]

# Transition state (first maximum after minimum)
search_region = pmf_kcal[idx_min:]
maxima = argrelextrema(search_region, np.greater, order=5)[0]
if len(maxima) > 0:
    idx_ts = maxima[0] + idx_min
else:
    # Fallback: midpoint between min and plateau
    idx_ts = idx_min + len(search_region) // 2

xi_ts = xi_bins[idx_ts]
barrier_wham = pmf_kcal[idx_ts] - pmf_kcal[idx_min]

# MBAR barrier
idx_min_mbar = np.argmin(pmf_mbar_kcal)
barrier_mbar = pmf_mbar_kcal[idx_ts] - pmf_mbar_kcal[idx_min_mbar]

# Bootstrap barrier distribution
bootstrap_barriers = []
if "bootstrap_pmfs_kj" in wham_bootstrap:
    for pmf_b in wham_bootstrap["bootstrap_pmfs_kj"]:
        pmf_b_kcal = pmf_b / KCAL_TO_KJ
        b_min = np.min(pmf_b_kcal)
        b_ts = pmf_b_kcal[idx_ts]
        bootstrap_barriers.append(b_ts - b_min)
    barrier_std = np.std(bootstrap_barriers)
else:
    barrier_std = np.sqrt(pmf_std_kcal[idx_ts]**2 + pmf_std_kcal[idx_min]**2)

print(f"Barrier (WHAM):  ΔW‡ = {barrier_wham:.2f} ± {barrier_std:.2f} kcal/mol")
print(f"Barrier (MBAR):  ΔW‡ = {barrier_mbar:.2f} kcal/mol")
print(f"ξ_min = {xi_min:.2f} nm, ξ‡ = {xi_ts:.2f} nm")

# Plateau check
plateau_mask = xi_bins > 3.5
if np.any(plateau_mask):
    plateau_slope = np.polyfit(xi_bins[plateau_mask], pmf_kcal[plateau_mask], 1)[0]
    print(f"Plateau slope: {plateau_slope:.2f} kcal/mol/nm (should be < 0.5)")
```

### Step 4: Block Averaging Convergence

```python
# Split each window into first/second half, recompute PMF
half_n = min(len(ts) for ts in xi_timeseries_list) // 2
xi_first_half = [ts[:half_n] for ts in xi_timeseries_list]
xi_second_half = [ts[half_n:2*half_n] for ts in xi_timeseries_list]

wham_first = solve_wham(xi_first_half, window_centers, spring_constants, 310.0, wham_config)
wham_second = solve_wham(xi_second_half, window_centers, spring_constants, 310.0, wham_config)

pmf_first_kcal = wham_first["pmf_kj_mol"] / KCAL_TO_KJ
pmf_second_kcal = wham_second["pmf_kj_mol"] / KCAL_TO_KJ

barrier_first = pmf_first_kcal[idx_ts] - np.min(pmf_first_kcal)
barrier_second = pmf_second_kcal[idx_ts] - np.min(pmf_second_kcal)
block_agreement = abs(barrier_first - barrier_second)
print(f"Block convergence: first={barrier_first:.2f}, second={barrier_second:.2f}, diff={block_agreement:.2f} kcal/mol")
```

### Step 5: Classification

```python
ea_exp = 10.5
if 4.6 <= barrier_wham <= 16.4:
    classification = "PASS"
elif 2.0 <= barrier_wham <= 20.0:
    classification = "MARGINAL"
else:
    classification = "FAIL"

# Check for no clear barrier
if barrier_wham < 1.0:
    classification = "FAIL"
    print("WARNING: No clear barrier detected")

# WHAM-MBAR agreement
wham_mbar_diff = abs(barrier_wham - barrier_mbar)
print(f"WHAM-MBAR difference: {wham_mbar_diff:.2f} kcal/mol")

results = {
    "experiment_id": "EXP-10", "feature_id": "F-10",
    "barrier_wham_kcal": float(barrier_wham), "barrier_mbar_kcal": float(barrier_mbar),
    "barrier_std_kcal": float(barrier_std),
    "ea_exp_kcal": float(ea_exp),
    "xi_min_nm": float(xi_min), "xi_ts_nm": float(xi_ts),
    "wham_mbar_agreement_kcal": float(wham_mbar_diff),
    "block_convergence_kcal": float(block_agreement),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-10 Classification: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: PMF with Barrier Annotation

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(xi_bins, pmf_kcal, 'b-', linewidth=2, label='WHAM')
ax.fill_between(xi_bins, pmf_kcal - 1.96*pmf_std_kcal, pmf_kcal + 1.96*pmf_std_kcal,
                alpha=0.15, color='blue')
ax.plot(xi_bins, pmf_mbar_kcal, 'r--', linewidth=1.5, label='MBAR')
ax.annotate('', xy=(xi_ts, pmf_kcal[idx_ts]), xytext=(xi_ts, pmf_kcal[idx_min]),
            arrowprops=dict(arrowstyle='<->', color='darkred', lw=2))
ax.text(xi_ts + 0.1, (pmf_kcal[idx_ts] + pmf_kcal[idx_min])/2,
        f'ΔW‡ = {barrier_wham:.1f} kcal/mol', fontsize=12, color='darkred')
ax.axhline(y=ea_exp + pmf_kcal[idx_min], color='gold', linestyle=':', linewidth=2,
           label=f'Exp. Ea = {ea_exp} kcal/mol')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-10: Activation Energy from PMF Barrier', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-10_barrier_height.png", dpi=300)
plt.close(fig)
```

### Figure 2: Bootstrap Barrier Distribution

```python
if bootstrap_barriers:
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.hist(bootstrap_barriers, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
    ax.axvline(x=ea_exp, color='gold', linestyle='--', linewidth=2,
               label=f'Experimental ({ea_exp} kcal/mol)')
    ax.axvspan(4.6, 16.4, alpha=0.1, color='green', label='95% CI')
    ax.set_xlabel('Barrier Height (kcal/mol)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('EXP-10: Bootstrap Barrier Distribution', fontsize=14)
    ax.legend()
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "EXP-10_bootstrap_barrier.png", dpi=300)
    plt.close(fig)
```

### Figure 3: Block Convergence

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(xi_bins, pmf_first_kcal - np.min(pmf_first_kcal), 'b-', alpha=0.7, label='First half')
ax.plot(xi_bins, pmf_second_kcal - np.min(pmf_second_kcal), 'r-', alpha=0.7, label='Second half')
ax.plot(xi_bins, pmf_kcal - np.min(pmf_kcal), 'k-', linewidth=2, label='Full')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-10: Block Convergence', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-10_block_convergence.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-10: BPTI-Trypsin Activation Energy — Results Report

**Experiment ID:** EXP-10  **Feature ID:** F-10  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## 1. Abstract
## 2. Introduction
## 3. Hypothesis
## 4. Methods (PMF barrier extraction, bootstrap, block averaging)
## 5. Controls
## 6. Results
### 6.1 ΔW‡ (WHAM): [value] ± [σ] kcal/mol
### 6.2 ΔW‡ (MBAR): [value] kcal/mol
### 6.3 Experimental Ea: 10.5 kcal/mol
### 6.4 WHAM-MBAR agreement: [value] kcal/mol
### 6.5 Block convergence: [value] kcal/mol
## 7. Discussion
## 8. Conclusions
## 9. Figures
## 10. References
## 11. Author Block
---
Author: Ryan Kamp
Affiliation: Dept. of Computer Science, University of Cincinnati
Contact:
Email: kamprj@mail.uc.edu
GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-10 loads the EXP-04 PMF for barrier height extraction — lightweight GPU usage.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-10"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

EXP04_DRIVE = Path("/content/drive/MyDrive/v3_gpu_results/EXP-04")
assert (EXP04_DRIVE / "outputs" / "umbrella").exists(), \
    "EXP-04 umbrella outputs required — run EXP-04 first"

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

# EXP-10 is PMF-based analysis: barrier height extraction + bootstrap.
# GPU is available but analysis is primarily NumPy/SciPy.
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

# ─── Checkpoint points for EXP-10 ───
# After PMF loading: save_analysis_checkpoint(pmf_data, "pmf_loaded")
# After barrier extraction: save_analysis_checkpoint(barrier_data, "barrier")
# After bootstrap: save_analysis_checkpoint(bootstrap_data, "bootstrap")
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

# ─── EXP-10 runtime estimates (A100 40GB) ───
# PMF loading + barrier extraction: ~5–10 minutes
# Bootstrap (200 samples): ~10–20 minutes
# Block convergence: ~5 minutes
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
# ─── EXP-10–specific recovery ───
# 1. No clear barrier in PMF:
#    - Symptom: argrelextrema finds no maxima
#    - Fix: use midpoint heuristic (already in Part 2 Step 3 fallback)
#    - Document as potential FAIL if barrier < 1 kcal/mol
# 2. Bootstrap distribution non-normal:
#    - Shapiro-Wilk test p < 0.05
#    - Report median ± IQR instead of mean ± std
# 3. WHAM-MBAR disagreement > 2 kcal/mol:
#    - Indicates poor sampling in barrier region
#    - Re-check EXP-04 US window overlap near transition state
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
