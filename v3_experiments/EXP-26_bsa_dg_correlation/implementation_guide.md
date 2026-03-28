# EXP-26: BSA–ΔG Binding Correlation — Implementation Guide

**Experiment ID:** EXP-26  
**Feature ID:** F-26 (benchmarks.md)  
**Category:** Biophysical (Semi-quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Tests whether the pipeline reproduces the empirical BSA–ΔG_bind correlation (ΔG ≈ −0.01 × BSA kcal/mol/Å²; Horton & Lewis 1992). Cross-experiment meta-analysis of EXP-01 through EXP-06 (ΔG) and EXP-16 (BSA). R² > 0.5 with slope in [−0.005, −0.020] validates internal thermodynamic consistency.

### 2. Hypotheses

- H₁: Correlation R² > 0.5 between BSA and ΔG_bind across 6 systems
- H₂: Slope in [−0.005, −0.020] kcal/mol/Å²

### 3. Classification (§25.1)

- PASS: R² > 0.5 AND slope in [−0.005, −0.020]
- MARGINAL: R² > 0.3 OR correct trend (negative slope)
- FAIL: R² < 0.1 or positive slope

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path
from scipy import stats

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-26_bsa_dg_correlation"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt

print("All imports successful.")
```

### Step 2: Collect Results from Upstream Experiments

```python
# Collect BSA and ΔG from completed experiment results
EXP_ROOT = PROJECT_ROOT / "v3_experiments"

systems = {
    "BPTI-trypsin":     {"bsa_exp": "EXP-04", "dg_exp": "EXP-04", "bsa_dir": "EXP-04_bpti_trypsin_dg_bind"},
    "SPINK7-KLK5":      {"bsa_exp": "EXP-01", "dg_exp": "EXP-01", "bsa_dir": "EXP-01_spink7_klk5_dg_bind"},
    "SPINK7-KLK12":     {"bsa_exp": "EXP-02", "dg_exp": "EXP-02", "bsa_dir": "EXP-02_spink7_klk12_dg_bind"},
    "LEKTI-KLK5":       {"bsa_exp": "EXP-03", "dg_exp": "EXP-03", "bsa_dir": "EXP-03_lekti_klk5_dg_panel"},
    "PSTI-chymo":       {"bsa_exp": "EXP-05", "dg_exp": "EXP-05", "bsa_dir": "EXP-05_psti_chymotrypsinogen_dg"},
    "SPINK1-trypsin":   {"bsa_exp": "EXP-06", "dg_exp": "EXP-06", "bsa_dir": "EXP-06_spink1_trypsin_dg_bind"},
}

# Expected values (used if results.json not yet available)
expected = {
    "BPTI-trypsin":   {"bsa": 1530, "dg": -18.0},
    "SPINK7-KLK5":    {"bsa": 1300, "dg": -9.4},
    "SPINK7-KLK12":   {"bsa": 1200, "dg": -7.9},
    "LEKTI-KLK5":     {"bsa": 1400, "dg": -11.8},
    "PSTI-chymo":     {"bsa": 1450, "dg": -14.7},
    "SPINK1-trypsin": {"bsa": 1350, "dg": -11.1},
}

data = {}
for name, info in systems.items():
    results_path = EXP_ROOT / info["bsa_dir"] / "results.json"
    if results_path.exists():
        with open(results_path) as f:
            res = json.load(f)
        bsa = res.get("bsa_mean_ang2", expected[name]["bsa"])
        dg = res.get("dg_bind_kcal", expected[name]["dg"])
    else:
        bsa = expected[name]["bsa"]
        dg = expected[name]["dg"]
    data[name] = {"bsa": bsa, "dg": dg}
    print(f"{name}: BSA={bsa:.0f} Å², ΔG={dg:.1f} kcal/mol")
```

### Step 3: Linear Regression

```python
names = list(data.keys())
bsa_values = np.array([data[n]["bsa"] for n in names])
dg_values = np.array([data[n]["dg"] for n in names])

# Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(bsa_values, dg_values)
r_squared = r_value**2

print(f"Slope: {slope:.5f} kcal/mol/Å²")
print(f"Intercept: {intercept:.2f} kcal/mol")
print(f"R²: {r_squared:.3f}")
print(f"p-value: {p_value:.4f}")
print(f"Literature slope: ~-0.01 kcal/mol/Å²")

# Spearman rank correlation
rho, p_spearman = stats.spearmanr(bsa_values, dg_values)
print(f"Spearman ρ: {rho:.3f}, p={p_spearman:.4f}")
```

### Step 4: Classification

```python
slope_ok = -0.020 <= slope <= -0.005
r2_strong = r_squared > 0.5
r2_moderate = r_squared > 0.3
negative_slope = slope < 0

if r2_strong and slope_ok:
    classification = "PASS"
elif r2_moderate or negative_slope:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-26", "feature_id": "F-26",
    "slope_kcal_per_ang2": float(slope),
    "intercept_kcal": float(intercept),
    "r_squared": float(r_squared),
    "p_value": float(p_value),
    "spearman_rho": float(rho),
    "n_systems": len(names),
    "data": {n: {"bsa": data[n]["bsa"], "dg": data[n]["dg"]} for n in names},
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-26: R²={r_squared:.3f}, slope={slope:.5f} → {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: BSA vs ΔG Scatter + Regression

```python
fig, ax = plt.subplots(figsize=(10, 7))
ax.scatter(bsa_values, dg_values, s=100, c='steelblue', edgecolors='black', zorder=5)

# Regression line
bsa_fit = np.linspace(min(bsa_values)-50, max(bsa_values)+50, 100)
dg_fit = slope * bsa_fit + intercept
ax.plot(bsa_fit, dg_fit, 'r-', linewidth=2,
         label=f'ΔG = {slope:.4f}×BSA + {intercept:.1f}\nR²={r_squared:.3f}')

# Literature slope
dg_lit = -0.01 * bsa_fit
ax.plot(bsa_fit, dg_lit, 'g--', linewidth=1.5, alpha=0.7,
         label='Literature: ΔG = −0.01×BSA')

# Label points
for i, name in enumerate(names):
    ax.annotate(name, (bsa_values[i], dg_values[i]), textcoords="offset points",
                 xytext=(5, 5), fontsize=9)

ax.set_xlabel('Buried Surface Area (Å²)', fontsize=13)
ax.set_ylabel('ΔG_bind (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-26: BSA–ΔG Correlation — {classification}', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-26_bsa_dg_scatter.png", dpi=300)
plt.close(fig)
```

### Figure 2: Rank Order Comparison

```python
fig, ax = plt.subplots(figsize=(10, 6))
sorted_by_dg = sorted(names, key=lambda n: data[n]["dg"])
dg_sorted = [data[n]["dg"] for n in sorted_by_dg]
bsa_sorted = [data[n]["bsa"] for n in sorted_by_dg]

x = np.arange(len(sorted_by_dg))
width = 0.35
bars1 = ax.bar(x - width/2, [-d for d in dg_sorted], width, label='|ΔG| (kcal/mol)',
                color='steelblue', edgecolor='black')
ax2 = ax.twinx()
bars2 = ax2.bar(x + width/2, bsa_sorted, width, label='BSA (Å²)',
                 color='coral', edgecolor='black')

ax.set_xticks(x)
ax.set_xticklabels(sorted_by_dg, rotation=30, ha='right')
ax.set_ylabel('|ΔG_bind| (kcal/mol)', fontsize=12)
ax2.set_ylabel('BSA (Å²)', fontsize=12)
ax.set_title(f'EXP-26: Rank Order (Spearman ρ={rho:.3f})', fontsize=14)

lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-26_rank_order.png", dpi=300)
plt.close(fig)
```

### Figure 3: Residuals Plot

```python
fig, ax = plt.subplots(figsize=(10, 5))
dg_predicted = slope * bsa_values + intercept
residuals = dg_values - dg_predicted

ax.scatter(bsa_values, residuals, s=80, c='steelblue', edgecolors='black', zorder=5)
ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
for i, name in enumerate(names):
    ax.annotate(name, (bsa_values[i], residuals[i]), textcoords="offset points",
                 xytext=(5, 5), fontsize=9)

ax.set_xlabel('BSA (Å²)', fontsize=12)
ax.set_ylabel('Residual ΔG (kcal/mol)', fontsize=12)
ax.set_title(f'EXP-26: Regression Residuals', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-26_residuals.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-26: BSA–ΔG Correlation — Results Report

**Experiment ID:** EXP-26  **Feature ID:** F-26  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Results
| Metric | Value | Criterion | Status |
|--------|-------|-----------|--------|
| Slope | [val] kcal/mol/Å² | [−0.005, −0.020] | [P/F] |
| R² | [val] | > 0.5 | [P/M/F] |
| Spearman ρ | [val] | > 0.5 | — |
| n systems | 6 | — | — |

## Figures
1. BSA vs ΔG scatter + regression
2. Rank order comparison
3. Residuals plot

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-26 is a meta-analysis experiment that loads results.json from EXP-01 through EXP-06 and EXP-16 — minimal GPU computation but requires all upstream experiments completed.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy

import os, sys
from pathlib import Path

EXP_ID = "EXP-26"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

# Verify upstream experiment results are available
GPU_RESULTS = Path("/content/drive/MyDrive/v3_gpu_results")
upstream = ["EXP-01", "EXP-02", "EXP-03", "EXP-04", "EXP-05", "EXP-06", "EXP-16"]
for exp in upstream:
    results_path = GPU_RESULTS / exp / "results.json"
    status = "✓" if results_path.exists() else "✗ MISSING"
    print(f"  {status} {exp}/results.json")

PROJECT_ROOT = Path("/content")
EXP_DIR = DRIVE_BASE
OUTPUT_DIR = DRIVE_BASE / "outputs"
FIGURES_DIR = DRIVE_BASE / "figures"

print(f"Colab environment ready. Drive base: {DRIVE_BASE}")
```

### §5.2 GPU Platform Selection

```python
# §5.2 GPU Platform Selection
# EXP-26 is a pure analysis/meta-analysis experiment.
# No simulation is run — GPU is not actively used for computation.
# The CUDA platform is available if needed for any trajectory re-analysis.
import openmm
from src.simulate.platform import select_platform

platform = select_platform("CUDA")
print(f"Platform: {platform.getName()} (available but not actively used)")
```

### §5.3 Checkpoint and Resume Integration

```python
# §5.3 Checkpoint and Resume Integration
import json, time
import numpy as np

CHECKPOINT_DIR = DRIVE_BASE / "checkpoints"

def save_analysis_checkpoint(data, phase_name):
    """Save analysis progress to Drive."""
    meta = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    meta.update({k: v for k, v in data.items() if not isinstance(v, np.ndarray)})
    with open(CHECKPOINT_DIR / f"{phase_name}.json", 'w') as f:
        json.dump(meta, f, indent=2, default=str)
    print(f"  Analysis checkpoint saved: {phase_name}")

# ─── Checkpoint points for EXP-26 ───
# After collecting upstream results: save_analysis_checkpoint(collected, "upstream")
# After regression: save_analysis_checkpoint(regression, "regression")
```

### §5.4 Progress Monitoring

```python
# §5.4 Progress Monitoring
# EXP-26 is lightweight (~5 minutes total).
# No GPU monitoring needed beyond initial verification.
import subprocess

def report_gpu_status():
    result = subprocess.run(
        ['nvidia-smi', '--query-gpu=name,memory.used,memory.total,utilization.gpu',
         '--format=csv,noheader'], capture_output=True, text=True)
    print(f"GPU: {result.stdout.strip()}")

# ─── EXP-26 runtime estimate ───
# Total: < 5 minutes (data collection + regression + figures)
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
# ─── EXP-26–specific recovery ───
# 1. Missing upstream results:
#    - Identify which EXP-XX/results.json files are absent
#    - Run missing upstream experiments first
# 2. Regression with < 6 data points:
#    - If some upstream experiments failed, report reduced N
#    - R² with small N may be unreliable — document limitation
# 3. Outlier detection:
#    - If one system dominates the regression, test with leave-one-out
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
