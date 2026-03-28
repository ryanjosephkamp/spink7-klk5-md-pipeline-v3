# EXP-06: SPINK1-Trypsin Binding Free Energy — Implementation Guide

**Experiment ID:** EXP-06  
**Feature ID:** F-06 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** KAZAL FAMILY BENCHMARK  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

This experiment computes the binding free energy of SPINK1 (identical to PSTI) against trypsin (Ki ≈ 7.6 nM, ΔG ≈ −11.1 kcal/mol), serving as a calibration point on the Kazal affinity ladder. The clinically significant N34S mutation does NOT alter Ki, providing a coupled negative control (EXP-32). PDB 1TGS (1.8 Å) provides starting geometry, with trypsinogen → trypsin conversion by activation peptide handling.

### 2. Hypothesis

**H₁:** ΔG_bind within 95% CI [−15.6, −6.6] kcal/mol.

**H₂:** Correctly ranked between SPINK7-KLK5 (−9.4) and PSTI-chymotrypsinogen (−14.7).

### 3. Background

SPINK1 measured at 37°C (310 K) matching pipeline default. Kuwata et al. (2002) measured Ki = 7.6 nM (recombinant) and 8.7 nM (natural), agreeing within 0.1 kcal/mol. Temperature match eliminates a common source of systematic error.

### 4. Protocol

PDB 1TGS, AMBER ff14SB, TIP3P, pH 7.4, 1.2 nm padding, 0.15 M NaCl. US: 51 windows (1.5–4.0 nm), 10 ns/window, 1000 kJ/mol/nm². WHAM/MBAR standard parameters. σ_exp = 0.5, σ_method = 2.0, σ_combined ≈ 2.3 kcal/mol. 95% CI [−15.6, −6.6].

### 5. Controls

Positive: EXP-05 (same inhibitor, different protease). Negative: EXP-32 (N34S mutation with ΔΔG ≈ 0).

### 6. Expected Outcomes

ΔG_bind = −11.1 kcal/mol, 95% CI [−15.6, −6.6].

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment and Directory Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-06_spink1_trypsin_dg_bind"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import (
    SystemConfig, MinimizationConfig, EquilibrationConfig,
    ProductionConfig, UmbrellaConfig, WHAMConfig, MBARConfig, KCAL_TO_KJ
)
from src.prep.pdb_fetch import fetch_pdb
from src.prep.pdb_clean import clean_structure
from src.prep.protonate import assign_protonation
from src.prep.topology import build_topology
from src.prep.solvate import solvate_system
from src.simulate.minimizer import minimize_energy
from src.simulate.equilibrate import run_nvt, run_npt
from src.simulate.production import run_production
from src.simulate.umbrella import (
    run_umbrella_campaign, generate_window_centers,
    diagnose_histogram_coverage, compute_overlap_matrix
)
from src.simulate.platform import select_platform
from src.analyze.wham import solve_wham, bootstrap_pmf_uncertainty
from src.analyze.mbar import solve_mbar, bootstrap_mbar_uncertainty
from src.analyze.structural import compute_rmsd
from src.analyze.trajectory import load_trajectory
from src.physics.units import kj_to_kcal, kbt
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Structure Acquisition and Trypsinogen → Trypsin Conversion

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

# Fetch PDB 1TGS (trypsinogen-PSTI complex)
complex_pdb = fetch_pdb("1TGS", data_dir)
complex_clean = clean_structure(
    complex_pdb, chains_to_keep=None,
    remove_heteroatoms=True, remove_waters=True, model_index=1
)
assert complex_clean.exists()

# NOTE: 1TGS contains trypsinogen (zymogen). For SPINK1-trypsin,
# the activation peptide (first ~6 residues) should be removed or
# the user should use a trypsin-only structure (e.g., 1S0Q or similar).
# For pipeline purposes, trypsinogen-PSTI is an acceptable approximation
# since the SPINK1 interface does not involve the activation peptide.

traj = md.load(str(complex_clean))
print(f"1TGS loaded: {traj.topology.n_chains} chains, {traj.topology.n_residues} residues")
```

### Step 3: Protonation, Topology, Solvation

```python
from openmm.app import PME, Simulation, PDBFile
from openmm import LangevinMiddleIntegrator, XmlSerializer
import openmm

complex_protonated = assign_protonation(complex_clean, ph=7.4, force_field="AMBER", use_propka=True)

sys_config = SystemConfig()
topology, system, modeller = build_topology(
    complex_protonated, sys_config, nonbonded_method=PME, nonbonded_cutoff_nm=1.0
)
modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)
print(f"Solvated: {n_waters} waters, {n_pos} Na+, {n_neg} Cl-")

solvated_pdb = OUTPUT_DIR / "solvated_complex.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
```

### Step 4: Minimization (IV-1)

```python
platform = select_platform()
integrator = LangevinMiddleIntegrator(
    310 * openmm.unit.kelvin, 1.0 / openmm.unit.picosecond, 0.002 * openmm.unit.picoseconds
)
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

e_initial = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(openmm.unit.kilojoules_per_mole)
minimize_energy(sim, MinimizationConfig())
e_final = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(openmm.unit.kilojoules_per_mole)
assert e_final < e_initial, f"IV-1 FAIL"
print(f"IV-1 PASS: E decreased {e_initial:.1f} → {e_final:.1f} kJ/mol")
```

### Step 5: Equilibration (IV-2, IV-3)

```python
eq_config = EquilibrationConfig()
eq_output = OUTPUT_DIR / "equilibration"
eq_output.mkdir(exist_ok=True)

nvt_results = run_nvt(sim, eq_config, eq_output)
t_avg = nvt_results.get("average_temperature_k", 0)
assert abs(t_avg - 310.0) < 5.0, f"IV-2 FAIL"
print(f"IV-2 PASS: T_avg = {t_avg:.1f} K")

npt_results = run_npt(sim, eq_config, eq_output)
rho = npt_results.get("average_density_g_per_cm3", 0)
assert 0.95 <= rho <= 1.05, f"IV-3 FAIL"
print(f"IV-3 PASS: density = {rho:.3f} g/cm³")

eq_state_path = eq_output / "equilibrated_state.xml"
sim.saveState(str(eq_state_path))
system_xml_path = eq_output / "system.xml"
with open(system_xml_path, "w") as f:
    f.write(XmlSerializer.serialize(system))
```

### Step 6: Production and Pull Group Setup

```python
prod_config = ProductionConfig()
prod_output = OUTPUT_DIR / "production"
prod_output.mkdir(exist_ok=True)
prod_results = run_production(sim, prod_config, prod_output)

prod_traj = load_trajectory(prod_output / "production_trajectory.dcd", solvated_pdb, stride=100)
topology_md = prod_traj.topology
spink1_indices = topology_md.select("chainid 1").tolist()
trypsin_indices = topology_md.select("chainid 0").tolist()
```

### Step 7: Umbrella Sampling, WHAM, MBAR

```python
us_config = UmbrellaConfig()
us_output = OUTPUT_DIR / "umbrella"
us_output.mkdir(exist_ok=True)
window_centers = generate_window_centers(us_config)

us_results = run_umbrella_campaign(
    equilibrated_state_path=eq_state_path, system_xml_path=system_xml_path,
    config=us_config, pull_group_1=spink1_indices, pull_group_2=trypsin_indices,
    output_dir=us_output, pdb_path=solvated_pdb, topology_pdb_path=solvated_pdb,
)

xi_timeseries_list = [np.load(us_output / f"window_{i:03d}" / "xi_timeseries.npy") for i in range(len(us_results))]

overlap = compute_overlap_matrix(xi_timeseries_list)
min_overlap = min(overlap[i, i+1] for i in range(len(overlap)-1))
assert min_overlap >= 0.10, f"IV-8 FAIL"

wham_config = WHAMConfig()
spring_constants = np.full(len(window_centers), us_config.spring_constant_kj_mol_nm2)
wham_results = solve_wham(xi_timeseries_list, window_centers, spring_constants, 310.0, wham_config)
assert wham_results["converged"], "IV-9 FAIL"
wham_bootstrap = bootstrap_pmf_uncertainty(xi_timeseries_list, window_centers, spring_constants, 310.0, wham_config)
mbar_results = solve_mbar(xi_timeseries_list, window_centers, spring_constants, 310.0, MBARConfig())
```

### Step 8: ΔG Extraction, Classification, and EXP-32 Coupling

```python
pmf_wham = wham_results["pmf_kj_mol"]
xi_bins = wham_results["xi_bins_nm"]
pmf_std = wham_bootstrap["pmf_std_kj_mol"]
kT = kbt(310.0)
beta = 1.0 / kT
C_standard_nm3 = 1.0 / (1660.0 * 0.001)

bound_mask = xi_bins < 2.0
bulk_mask = xi_bins > 3.5
pmf_shifted = pmf_wham - pmf_wham[bulk_mask].min()
I_site = np.trapz(np.exp(-beta * pmf_shifted[bound_mask]) * xi_bins[bound_mask]**2, xi_bins[bound_mask])
I_bulk = np.trapz(np.exp(-beta * pmf_shifted[bulk_mask]) * xi_bins[bulk_mask]**2, xi_bins[bulk_mask])
dg_bind_kj = -kT * np.log(C_standard_nm3/(4*np.pi) * I_site) + kT * np.log(C_standard_nm3/(4*np.pi) * I_bulk)
dg_bind_kcal = kj_to_kcal(dg_bind_kj)

dg_exp = -11.1
sigma_exp = 0.5; sigma_comp = kj_to_kcal(np.mean(pmf_std)); sigma_method = 2.0
sigma_combined = np.sqrt(sigma_exp**2 + sigma_comp**2 + sigma_method**2)
discrepancy = abs(dg_bind_kcal - dg_exp)

z = 1.96
if sigma_comp > sigma_exp:
    classification = "INCONCLUSIVE"
elif discrepancy <= z * sigma_combined:
    classification = "PASS"
elif discrepancy <= 2 * z * sigma_combined:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-06", "feature_id": "F-06",
    "system": "SPINK1-trypsin (PDB 1TGS)",
    "dg_bind_wham_kcal": float(dg_bind_kcal), "dg_exp_kcal": float(dg_exp),
    "sigma_combined": float(sigma_combined), "classification": classification,
    "note": "Coupled with EXP-32 (N34S control); WT results stored for ΔΔG comparison"
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-06 Classification: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: PMF Profile

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
pmf_kcal = pmf_wham / KCAL_TO_KJ
std_kcal = pmf_std / KCAL_TO_KJ
ax.plot(xi_bins, pmf_kcal, 'b-', linewidth=2, label='WHAM')
ax.fill_between(xi_bins, pmf_kcal - 1.96*std_kcal, pmf_kcal + 1.96*std_kcal, alpha=0.2, color='blue')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-06: SPINK1-Trypsin PMF', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-06_pmf_profile.png", dpi=300)
plt.close(fig)
```

### Figure 2: Predicted vs Experimental

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.bar(['Experimental', 'Pipeline'], [dg_exp, dg_bind_kcal],
       color=['gold', 'steelblue'], edgecolor='black')
ax.axhspan(-15.6, -6.6, alpha=0.1, color='green', label='95% CI')
ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-06: Predicted vs Experimental — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-06_predicted_vs_experimental.png", dpi=300)
plt.close(fig)
```

### Figure 3: Kazal Affinity Ladder (All Thermodynamic)

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 6))
systems = ['SPINK7-KLK5\n(EXP-01)', 'SPINK1-Trypsin\n(EXP-06)',
           'LEKTI-KLK5\n(EXP-03)', 'PSTI-Chymo\n(EXP-05)', 'BPTI-Trypsin\n(EXP-04)']
dg_exp_vals = [-9.4, -11.1, -11.8, -14.7, -18.0]
ax.barh(systems, dg_exp_vals, color='gold', edgecolor='black', alpha=0.7, label='Experimental')
ax.set_xlabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title('EXP-06: Protease-Inhibitor Affinity Ladder', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-06_affinity_ladder.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-06: SPINK1-Trypsin Binding Free Energy — Results Report

**Experiment ID:** EXP-06  **Feature ID:** F-06  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL/INCONCLUSIVE]

## 1. Abstract
## 2. Introduction/Background
## 3. Hypothesis
## 4. Methods
### 4.1 System Preparation (PDB 1TGS, trypsinogen→trypsin)
### 4.2 Umbrella Sampling
### 4.3 WHAM/MBAR
### 4.4 Deviations
## 5. Controls (coupled with EXP-32 N34S)
## 6. Results
### 6.1 ΔG_bind (WHAM): [value] kcal/mol
### 6.2 Experimental: −11.1 ± 0.5 kcal/mol
### 6.3 Classification
### 6.4 Ranking Position
## 7. Discussion
## 8. Conclusions
## 9. Figures
## 10. References
- Kuwata et al. (2002)
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

EXP_ID = "EXP-06"
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

# EXP-06 uses SPINK1 (Kazal-type) with trypsin — same pipeline as EXP-04.
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

# ─── Checkpoint points for EXP-06 ───
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

# ─── EXP-06 runtime estimates (A100 40GB) ───
# Production (100 ns): ~2–4 hours
# SMD (50 reps × 3 ns): ~3–6 hours
# US (51 windows × 10 ns): ~10–20 hours
# Total: ~16–30 hours
# Note: coupled with EXP-32 (N34S negative control) —
#   run EXP-06 first, then EXP-32 uses same system setup
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

# ─── EXP-06–specific recovery ───
# 1. SPINK1 Kazal domain model quality:
#    - If homology model RMSD drifts > 3 Å, check protonation states
# 2. Trypsinogen→trypsin conversion artifacts:
#    - Verify activation peptide processing was correct
# 3. EXP-32 dependency: save WT equilibrated state for N34S mutation
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
