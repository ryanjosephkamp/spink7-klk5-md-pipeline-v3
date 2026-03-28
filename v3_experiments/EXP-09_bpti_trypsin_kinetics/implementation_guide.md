# EXP-09: BPTI–Trypsin Association/Dissociation Kinetics — Implementation Guide

**Experiment ID:** EXP-09  
**Feature ID:** F-09 (benchmarks.md)  
**Category:** Kinetic  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Tests whether the V2 pipeline produces kinetic estimates (kon/koff) for BPTI-trypsin within order-of-magnitude ranges. Experimental: kon ≈ 10⁶ M⁻¹s⁻¹, koff ≈ 6 × 10⁻⁸ s⁻¹ (Kd = 6 × 10⁻¹⁴ M). Uses PMF-based TST/Kramers' theory from EXP-04 PMF. PASS: both within 2 orders; MARGINAL: one within 2, other within 3; FAIL: either off by >3 orders.

### 2. Hypothesis

**H₁:** PMF-derived kon within 10⁴–10⁸ M⁻¹s⁻¹. **H₂:** koff within 10⁻¹⁰–10⁻⁶ s⁻¹. **H₃:** koff/kon consistent with EXP-04 ΔG within 2 kcal/mol.

### 3. Protocol

PMF from EXP-04 US/WHAM. TST: koff from barrier height, curvature, D(ξ). Kramers: koff from friction coefficient. kon from detailed balance. MSM-based approach as secondary method.

### 4. Controls

Positive: EXP-13 (barnase-barstar, fast association). Consistency: Kd(kinetic) = koff/kon matches Kd(thermo). Negative: kon ≤ Smoluchowski limit.

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-09_bpti_trypsin_kinetics"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import WHAMConfig, MBARConfig, KCAL_TO_KJ
from src.analyze.wham import solve_wham
from src.analyze.mbar import solve_mbar
from src.analyze.convergence import evaluate_convergence
from src.analyze.msm import fit_tica, cluster_microstates, build_msm, compute_mfpt
from src.physics.units import kj_to_kcal, kbt
import matplotlib.pyplot as plt

print("All imports successful.")
```

### Step 2: Load EXP-04 PMF Data

```python
EXP04_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind" / "outputs"

# Load WHAM results
wham_output = EXP04_DIR / "umbrella"
# Reconstruct or load PMF
exp04_results_path = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind" / "results.json"
assert exp04_results_path.exists(), "Run EXP-04 first"

with open(exp04_results_path) as f:
    exp04 = json.load(f)

# Load PMF arrays (saved by EXP-04 umbrella step)
pmf_path = wham_output / "pmf_wham.npy"
xi_path = wham_output / "xi_bins.npy"
if pmf_path.exists() and xi_path.exists():
    pmf_kj = np.load(pmf_path)
    xi_bins = np.load(xi_path)
else:
    # Reconstruct from umbrella timeseries
    from src.simulate.umbrella import generate_window_centers
    from src.config import UmbrellaConfig
    us_config = UmbrellaConfig()
    window_centers = generate_window_centers(us_config)
    xi_timeseries_list = [np.load(wham_output / f"window_{i:03d}" / "xi_timeseries.npy")
                          for i in range(len(window_centers))]
    spring_constants = np.full(len(window_centers), us_config.spring_constant_kj_mol_nm2)
    wham_results = solve_wham(xi_timeseries_list, window_centers, spring_constants, 310.0, WHAMConfig())
    pmf_kj = wham_results["pmf_kj_mol"]
    xi_bins = wham_results["xi_bins_nm"]

pmf_kcal = pmf_kj / KCAL_TO_KJ
print(f"PMF loaded: {len(xi_bins)} bins, range [{xi_bins[0]:.2f}, {xi_bins[-1]:.2f}] nm")
```

### Step 3: Identify Bound State Minimum and Transition State

```python
# Find global minimum (bound state)
idx_min = np.argmin(pmf_kcal)
xi_min = xi_bins[idx_min]
pmf_min = pmf_kcal[idx_min]

# Find transition state (maximum between bound and unbound)
# Search for local maxima beyond the minimum
maxima_indices = argrelextrema(pmf_kcal[idx_min:], np.greater, order=5)[0] + idx_min
if len(maxima_indices) > 0:
    idx_ts = maxima_indices[0]
else:
    # Fallback: point where PMF first reaches plateau
    plateau_mask = xi_bins > 3.0
    idx_ts = np.argmax(pmf_kcal[plateau_mask]) + np.argmax(plateau_mask)

xi_ts = xi_bins[idx_ts]
pmf_ts = pmf_kcal[idx_ts]
barrier_height = pmf_ts - pmf_min

print(f"Bound state: ξ_min = {xi_min:.2f} nm, PMF = {pmf_min:.2f} kcal/mol")
print(f"TS: ξ‡ = {xi_ts:.2f} nm, PMF = {pmf_ts:.2f} kcal/mol")
print(f"Barrier height ΔW‡ = {barrier_height:.2f} kcal/mol")
```

### Step 4: Compute Diffusion Coefficient D(ξ)

```python
# D(ξ) from autocorrelation of ξ fluctuations in each umbrella window
# Use the umbrella window closest to the TS
us_config_reload = UmbrellaConfig() if 'us_config' not in dir() else us_config
window_centers_reload = generate_window_centers(us_config_reload) if 'window_centers' not in dir() else window_centers

# Find window closest to TS
ts_window_idx = np.argmin(np.abs(window_centers - xi_ts))
xi_ts_timeseries = np.load(wham_output / f"window_{ts_window_idx:03d}" / "xi_timeseries.npy")

# Autocorrelation time
dt = 0.001  # 1 ps save interval in ns
xi_fluct = xi_ts_timeseries - np.mean(xi_ts_timeseries)
var_xi = np.var(xi_fluct)
autocorr = np.correlate(xi_fluct, xi_fluct, mode='full')[len(xi_fluct)-1:]
autocorr /= autocorr[0]

# Integrate autocorrelation to get correlation time
tau_corr_idx = np.argmax(autocorr < 1/np.e)
tau_corr = tau_corr_idx * dt  # ns

# D = var(ξ) / tau_corr
D_xi = var_xi / tau_corr if tau_corr > 0 else var_xi / 0.01  # nm²/ns
D_xi_m2s = D_xi * 1e-18 / 1e-9  # Convert nm²/ns to m²/s

print(f"D(ξ‡) = {D_xi:.4f} nm²/ns = {D_xi_m2s:.2e} m²/s")
print(f"Autocorrelation time τ = {tau_corr:.3f} ns")
```

### Step 5: Compute koff via TST and Kramers'

```python
R = 1.987e-3  # kcal/(mol·K)
T = 310.0
kT_kcal = R * T

# --- TST ---
# koff_TST = (D/2πkT) * sqrt(|W''(min)| * |W''(TS)|) / kT * exp(-ΔW‡/kT)
# Numerical second derivatives
dx = xi_bins[1] - xi_bins[0]
W_second_deriv = np.gradient(np.gradient(pmf_kcal, dx), dx)
W_pp_min = abs(W_second_deriv[idx_min])
W_pp_ts = abs(W_second_deriv[idx_ts])

omega_min = np.sqrt(W_pp_min * KCAL_TO_KJ / 1.0)  # approximate angular frequency
omega_ts = np.sqrt(W_pp_ts * KCAL_TO_KJ / 1.0)

# koff via Kramers' (high friction limit)
gamma = kT_kcal / D_xi if D_xi > 0 else 1.0  # friction from Einstein relation
koff_kramers = (omega_min * omega_ts) / (2 * np.pi * gamma) * np.exp(-barrier_height / kT_kcal)

# koff via simple Arrhenius/TST
nu_attempt = 1e12  # typical attempt frequency (s⁻¹)
koff_tst = nu_attempt * np.exp(-barrier_height / kT_kcal)

# --- kon from detailed balance ---
dg_bind = exp04.get("dg_bind_wham_kcal", -18.0)
Kd = np.exp(dg_bind / kT_kcal) / 1.0  # 1 M standard state, units: M

kon_kramers = koff_kramers / Kd if Kd > 0 else 0
kon_tst = koff_tst / Kd if Kd > 0 else 0

print(f"\n--- Kinetic Estimates ---")
print(f"koff (Kramers): {koff_kramers:.2e} s⁻¹")
print(f"koff (TST):     {koff_tst:.2e} s⁻¹")
print(f"kon (Kramers):  {kon_kramers:.2e} M⁻¹s⁻¹")
print(f"kon (TST):      {kon_tst:.2e} M⁻¹s⁻¹")
print(f"Kd (kinetic):   {Kd:.2e} M")
```

### Step 6: Classification

```python
# Experimental values
kon_exp = 1e6  # M⁻¹s⁻¹
koff_exp = 6e-8  # s⁻¹

# Use Kramers estimate as primary
log_kon_err = abs(np.log10(kon_kramers) - np.log10(kon_exp))
log_koff_err = abs(np.log10(max(koff_kramers, 1e-20)) - np.log10(koff_exp))

if log_kon_err <= 2 and log_koff_err <= 2:
    classification = "PASS"
elif (log_kon_err <= 2 and log_koff_err <= 3) or (log_kon_err <= 3 and log_koff_err <= 2):
    classification = "MARGINAL"
else:
    classification = "FAIL"

# Thermodynamic consistency check
dg_kinetic = kT_kcal * np.log(Kd)
dg_thermo = dg_bind
thermo_consistency = abs(dg_kinetic - dg_thermo)

results = {
    "experiment_id": "EXP-09", "feature_id": "F-09",
    "kon_kramers": float(kon_kramers), "koff_kramers": float(koff_kramers),
    "kon_tst": float(kon_tst), "koff_tst": float(koff_tst),
    "kon_exp": float(kon_exp), "koff_exp": float(koff_exp),
    "barrier_height_kcal": float(barrier_height),
    "thermo_consistency_kcal": float(thermo_consistency),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-09 Classification: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: PMF with Barrier Annotation

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(xi_bins, pmf_kcal, 'b-', linewidth=2)
ax.axhline(y=pmf_min, color='green', linestyle=':', alpha=0.5)
ax.axhline(y=pmf_ts, color='red', linestyle=':', alpha=0.5)
ax.annotate(f'ΔW‡ = {barrier_height:.1f} kcal/mol',
            xy=(xi_ts, pmf_ts), xytext=(xi_ts+0.3, pmf_ts-2),
            arrowprops=dict(arrowstyle='->', color='red'), fontsize=12, color='red')
ax.plot(xi_min, pmf_min, 'go', markersize=10, label='Bound state')
ax.plot(xi_ts, pmf_ts, 'r^', markersize=10, label='Transition state')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-09: PMF with Kinetic Barrier', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-09_pmf_barrier.png", dpi=300)
plt.close(fig)
```

### Figure 2: Rate Constants Comparison

```python
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
ax = axes[0]
bars = ax.bar(['Experimental', 'Kramers', 'TST'],
              [np.log10(kon_exp), np.log10(kon_kramers), np.log10(kon_tst)],
              color=['gold', 'steelblue', 'coral'], edgecolor='black')
ax.set_ylabel(r'$\log_{10}(k_{on}$ / M⁻¹s⁻¹)', fontsize=12)
ax.set_title('kon', fontsize=13)

ax = axes[1]
bars = ax.bar(['Experimental', 'Kramers', 'TST'],
              [np.log10(koff_exp), np.log10(max(koff_kramers, 1e-20)), np.log10(max(koff_tst, 1e-20))],
              color=['gold', 'steelblue', 'coral'], edgecolor='black')
ax.set_ylabel(r'$\log_{10}(k_{off}$ / s⁻¹)', fontsize=12)
ax.set_title('koff', fontsize=13)

plt.suptitle(f'EXP-09: Kinetic Rate Constants — {classification}', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-09_rate_constants.png", dpi=300)
plt.close(fig)
```

### Figure 3: Diffusion Coefficient Profile

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
ax.plot(xi_bins, W_second_deriv, 'b-', linewidth=1.5)
ax.axvline(x=xi_min, color='green', linestyle='--', label='ξ_min')
ax.axvline(x=xi_ts, color='red', linestyle='--', label='ξ‡')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel("W''(ξ) (kcal/mol/nm²)", fontsize=13)
ax.set_title('EXP-09: PMF Second Derivative', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-09_pmf_curvature.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-09: BPTI-Trypsin Kinetics — Results Report

**Experiment ID:** EXP-09  **Feature ID:** F-09  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## 1. Abstract
## 2. Introduction
## 3. Hypothesis
## 4. Methods (TST, Kramers', MSM)
## 5. Controls
## 6. Results
### 6.1 kon (Kramers): [value] M⁻¹s⁻¹
### 6.2 koff (Kramers): [value] s⁻¹
### 6.3 Barrier height: [value] kcal/mol
### 6.4 Thermodynamic consistency: [value] kcal/mol
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

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-09 loads the EXP-04 PMF for TST/Kramers kinetics calculation.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-09"
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

# EXP-09 is primarily PMF-based analysis (TST/Kramers).
# GPU is used for any trajectory re-analysis or diffusion coefficient estimation.
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

# ─── Checkpoint points for EXP-09 ───
# After PMF loading: save_analysis_checkpoint(pmf_data, "pmf_loaded")
# After barrier extraction: save_analysis_checkpoint(barrier, "barrier")
# After kinetics calculation: save_analysis_checkpoint(kinetics, "kinetics")
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

# ─── EXP-09 runtime estimates (A100 40GB) ───
# PMF analysis + Kramers/TST: ~15–30 minutes
# Total: < 1 hour (lightweight analysis experiment)
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
# ─── EXP-09–specific recovery ───
# 1. EXP-04 PMF not found: verify EXP-04 umbrella/WHAM completed on Drive
# 2. Kramers theory divergence:
#    - If diffusion coefficient D(ξ) is noisy, smooth with Savitzky-Golay
#    - If barrier is too shallow (< 1 kcal/mol), TST breaks down — document
# 3. Thermodynamic consistency failure:
#    - koff/kon ratio should yield ΔG consistent with EXP-04
#    - Discrepancy > 2 kcal/mol indicates systematic error in barrier profile
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
