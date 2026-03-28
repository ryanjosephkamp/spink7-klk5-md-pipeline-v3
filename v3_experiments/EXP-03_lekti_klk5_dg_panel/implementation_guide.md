# EXP-03: LEKTI-KLK5 Binding Free Energy Panel — Implementation Guide

**Experiment ID:** EXP-03  
**Feature ID:** F-03 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** CROSS-VALIDATION  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

This experiment determines the binding free energy of a representative LEKTI (SPINK5) single Kazal domain–KLK5 complex using US/WHAM, benchmarked against the most precise single-domain Ki value: rLEKTI(1–6) Ki = 2.35 ± 0.22 nM (ΔG = −11.8 kcal/mol) from Borgono et al. (2007). LEKTI domains share ~35% sequence identity with SPINK7 and use the same canonical Kazal inhibitory mechanism. Testing a related Kazal-KLK5 interaction validates whether the pipeline generalizes across Kazal family members and whether systematic biases observed for SPINK7-KLK5 (EXP-01) are system-specific or method-inherent.

### 2. Hypothesis

**H₁:** The V2 pipeline's US/WHAM estimate of ΔG_bind for an isolated LEKTI Kazal domain (domain 6 or domains 1–6 fragment) against KLK5 will fall within the 95% CI [−16.2, −7.4] kcal/mol.

**H₂:** The computed ΔG_bind will be more negative (tighter binding) than the SPINK7-KLK5 prediction (EXP-01), consistent with the experimental observation that LEKTI fragments bind KLK5 more tightly (Ki = 2.35 nM vs. 132 nM, ΔΔG ≈ 2.4 kcal/mol).

### 3. Background and Rationale

#### 3.1 Scientific Context

LEKTI is a 15-domain Kazal-type inhibitor that regulates KLK5 activity in the skin. Three independent laboratories have measured LEKTI-KLK5 inhibition (Borgono et al. 2007, Deraison et al. 2007, Egelrud et al. 2005), providing multi-group cross-validation. Single-domain Ki values range from 2.35–118.7 nM, while multi-domain fragments achieve sub-picomolar apparent KD by avidity effects. For MD simulation, single-domain Ki values are the appropriate targets.

#### 3.2 What This Reveals

Cross-validation with a structurally related but distinct inhibitor-protease pair tests pipeline transferability. If the pipeline performs well on LEKTI-KLK5 but poorly on SPINK7-KLK5, the issue is SPINK7-specific (likely model quality). If both fail, the issue is methodological.

### 4. Experimental Protocol

#### 4.1 System Preparation

| Parameter | Value |
|-----------|-------|
| LEKTI domain structure | Homology model of LEKTI domain 6 (or available Kazal domain structure) |
| KLK5 structure | PDB 2PSX (chain A) |
| Complex model | Docking based on canonical Kazal-protease binding mode |
| Force field | AMBER ff14SB (`amber14-all.xml`) |
| Water model | TIP3P (`amber14/tip3p.xml`) |
| pH | 7.4 |
| Box padding | 1.2 nm |
| Ionic strength | 0.15 M NaCl |
| Box shape | Cubic |

#### 4.2 Minimization and Equilibration

Same parameters as EXP-01: Minimization 10,000 steps, tolerance 10.0 kJ/mol/nm. NVT 500 ps at 310 K, friction 1.0 ps⁻¹, timestep 2 fs. NPT 1000 ps at 310 K and 1.0 atm, barostat interval 25. Restraint 1000 kJ/mol/nm², save interval 10 ps.

#### 4.3 Umbrella Sampling

| Parameter | Value | Source |
|-----------|-------|--------|
| ξ range | 1.5–4.0 nm | `config.py: UmbrellaConfig` |
| Window spacing | 0.05 nm (51 windows) | `config.py: UmbrellaConfig.window_spacing_nm` |
| Spring constant | 1000 kJ/mol/nm² | `config.py: UmbrellaConfig.spring_constant_kj_mol_nm2` |
| Per-window duration | 10.0 ns | `config.py: UmbrellaConfig.per_window_duration_ns` |
| Pre-positioning velocity | 0.01 nm/ps | `config.py: UmbrellaConfig.pre_position_velocity_nm_per_ps` |
| Equilibration per window | 200 ps | `config.py: UmbrellaConfig.equilibration_duration_ps` |
| Save interval | 1.0 ps | `config.py: UmbrellaConfig.save_interval_ps` |

#### 4.4 WHAM and MBAR Analysis

WHAM: tolerance 10⁻⁶, max iterations 100,000, bootstraps 200, bins 200. MBAR: solver "robust", tolerance 10⁻⁷, max iterations 10,000, bootstraps 200, bins 200.

#### 4.5 Statistical Comparison

σ_exp = 0.06 kcal/mol, σ_method = 2.0 kcal/mol, σ_comp from bootstrap.

### 5. Control Conditions

**Positive Control:** EXP-01 (SPINK7-KLK5) — same protease, related inhibitor.

**Negative Controls:** All invariants IV-1 through IV-9 satisfied. LEKTI domain maintains Kazal fold (RMSD < 3 Å). ΔG should be negative.

### 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind (LEKTI D6–KLK5) | −11.8 kcal/mol | [−16.2, −7.4] |
| Ranking vs. SPINK7-KLK5 | More favorable ΔG | LEKTI binds tighter |

### 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation Implied | Severity |
|-------------|--------------|-------------------|----------|
| LEKTI domain model error | Non-canonical binding | Homology model quality | High |
| Ranking inversion | LEKTI predicted weaker than SPINK7 | Cannot rank inhibitors | High |

### 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After LEKTI model | Kazal fold; disulfides present | Valid topology |
| After docking | RSL in KLK5 active site | Canonical binding |
| After equilibration | IV-2 through IV-7 | All pass |
| After US/WHAM | IV-8, IV-9; smooth PMF | Converged |

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment and Directory Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-03_lekti_klk5_dg_panel"
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
from src.physics.units import kj_to_kcal, nm_to_angstrom, kbt
from src.visualization.plot_pmf import plot_pmf
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

**If this step fails:** Verify virtual environment and dependencies.

---

### Step 2: Structure Acquisition

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

# Fetch KLK5 structure
klk5_pdb = fetch_pdb("2PSX", data_dir)
klk5_clean = clean_structure(klk5_pdb, chains_to_keep=["A"],
                              remove_heteroatoms=True, remove_waters=True, model_index=1)

# LEKTI domain 6 — requires homology model or AlphaFold
# LEKTI/SPINK5 UniProt: Q9NQ38
from src.prep.pdb_fetch import fetch_alphafold
try:
    lekti_pdb = fetch_alphafold("Q9NQ38", data_dir)
    # Extract domain 6 region (approximately residues 300-360 of full LEKTI)
    lekti_clean = clean_structure(lekti_pdb, chains_to_keep=["A"],
                                   remove_heteroatoms=True, remove_waters=True, model_index=1)
    print(f"LEKTI AlphaFold model: {lekti_clean}")
except Exception as e:
    print(f"AlphaFold fetch failed: {e}")
    lekti_clean = data_dir / "lekti_domain6.pdb"
    print(f"Place LEKTI domain 6 model at: {lekti_clean}")

assert klk5_clean.exists(), "KLK5 cleaned PDB not found"
```

**If this step fails:** If no LEKTI model is available, attempt Swiss-Model or similar. If unsuccessful, classify as INCONCLUSIVE.

---

### Step 3: Complex Assembly, Protonation, Solvation, Minimization

```python
from openmm.app import PME, Simulation, PDBFile
from openmm import LangevinMiddleIntegrator, XmlSerializer
import openmm

docked_complex = data_dir / "lekti_klk5_complex.pdb"
if not docked_complex.exists():
    print(f"WARNING: Place docked LEKTI-KLK5 complex at: {docked_complex}")

complex_protonated = assign_protonation(docked_complex, ph=7.4, force_field="AMBER", use_propka=True)

sys_config = SystemConfig()
topology, system, modeller = build_topology(complex_protonated, sys_config, nonbonded_method=PME, nonbonded_cutoff_nm=1.0)
modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)

solvated_pdb = OUTPUT_DIR / "solvated_complex.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

platform = select_platform()
integrator = LangevinMiddleIntegrator(310*openmm.unit.kelvin, 1.0/openmm.unit.picosecond, 0.002*openmm.unit.picoseconds)
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

e_initial = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(openmm.unit.kilojoules_per_mole)
minimize_energy(sim, MinimizationConfig())
e_final = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(openmm.unit.kilojoules_per_mole)
assert e_final < e_initial, f"IV-1 FAIL"
print(f"IV-1 PASS: E decreased by {e_initial - e_final:.2f} kJ/mol")
```

---

### Step 4: Equilibration

```python
eq_config = EquilibrationConfig()
eq_output = OUTPUT_DIR / "equilibration"
eq_output.mkdir(exist_ok=True)

nvt_results = run_nvt(sim, eq_config, eq_output)
t_avg = nvt_results.get("average_temperature_k", 0)
assert abs(t_avg - 310.0) < 5.0, f"IV-2 FAIL"

npt_results = run_npt(sim, eq_config, eq_output)
rho = npt_results.get("average_density_g_per_cm3", 0)
assert 0.95 <= rho <= 1.05, f"IV-3 FAIL"

eq_state_path = eq_output / "equilibrated_state.xml"
sim.saveState(str(eq_state_path))
system_xml_path = eq_output / "system.xml"
with open(system_xml_path, "w") as f:
    f.write(XmlSerializer.serialize(system))
```

---

### Step 5: Production MD and Pull Group Identification

```python
prod_config = ProductionConfig()
prod_output = OUTPUT_DIR / "production"
prod_output.mkdir(exist_ok=True)
prod_results = run_production(sim, prod_config, prod_output)

prod_traj = load_trajectory(prod_output / "production_trajectory.dcd", solvated_pdb, stride=100)
topology_md = prod_traj.topology
chain_a_indices = topology_md.select("chainid 0").tolist()  # LEKTI
chain_b_indices = topology_md.select("chainid 1").tolist()  # KLK5
```

---

### Step 6: Umbrella Sampling and PMF Reconstruction

```python
us_config = UmbrellaConfig()
us_output = OUTPUT_DIR / "umbrella"
us_output.mkdir(exist_ok=True)
window_centers = generate_window_centers(us_config)

us_results = run_umbrella_campaign(
    equilibrated_state_path=eq_state_path, system_xml_path=system_xml_path,
    config=us_config, pull_group_1=chain_a_indices, pull_group_2=chain_b_indices,
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

---

### Step 7: ΔG Extraction and Classification

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

dg_exp = -11.8
sigma_exp = 0.06
sigma_comp = kj_to_kcal(np.mean(pmf_std))
sigma_method = 2.0
sigma_combined = np.sqrt(sigma_exp**2 + sigma_comp**2 + sigma_method**2)
discrepancy = abs(dg_bind_kcal - dg_exp)

if sigma_comp > sigma_exp:
    classification = "INCONCLUSIVE"
elif discrepancy <= 1.96 * sigma_combined:
    classification = "PASS"
elif discrepancy <= 2 * 1.96 * sigma_combined:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-03", "feature_id": "F-03",
    "dg_bind_wham_kcal": float(dg_bind_kcal), "dg_exp_kcal": float(dg_exp),
    "sigma_combined": float(sigma_combined), "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"Classification: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Experimental Design — Kazal Family Comparison Schematic

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
systems = ['SPINK7-KLK5\n(EXP-01)', 'LEKTI D6-KLK5\n(EXP-03)', 'SPINK1-Trypsin\n(EXP-06)']
dg_values = [-9.4, -11.8, -11.2]
colors = ['steelblue', 'coral', 'mediumseagreen']
ax.bar(systems, dg_values, color=colors, edgecolor='black', alpha=0.85)
ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title('EXP-03: Kazal Family ΔG Comparison (Experimental Values)', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-03_kazal_comparison_schematic.png", dpi=300)
plt.close(fig)
```

### Figure 2: PMF Profile

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
pmf_kcal = pmf_wham / KCAL_TO_KJ
std_kcal = pmf_std / KCAL_TO_KJ
ax.plot(xi_bins, pmf_kcal, 'b-', linewidth=2, label='WHAM')
ax.fill_between(xi_bins, pmf_kcal - 1.96*std_kcal, pmf_kcal + 1.96*std_kcal, alpha=0.2, color='blue')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-03: LEKTI-KLK5 PMF Profile', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-03_pmf_profile.png", dpi=300)
plt.close(fig)
```

### Figure 3: Predicted vs Experimental

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.bar(['Experimental', 'Pipeline (WHAM)'], [dg_exp, dg_bind_kcal],
       color=['gold', 'steelblue'], edgecolor='black')
ax.axhspan(-16.2, -7.4, alpha=0.1, color='green', label='95% CI')
ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title('EXP-03: Predicted vs Experimental', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-03_predicted_vs_experimental.png", dpi=300)
plt.close(fig)
```

### Figure 4: Umbrella Window Overlap

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 5))
for xi_ts in xi_timeseries_list:
    ax.hist(xi_ts, bins=50, alpha=0.3, density=True)
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('Probability density', fontsize=13)
ax.set_title('EXP-03: Umbrella Window Overlap', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-03_umbrella_overlap.png", dpi=300)
plt.close(fig)
```

### Figure 5: Ranking Verification — LEKTI vs SPINK7

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
exp01_path = PROJECT_ROOT / "v3_experiments" / "EXP-01_spink7_klk5_dg_bind" / "results.json"
if exp01_path.exists():
    with open(exp01_path) as f:
        exp01 = json.load(f)
    dg_spink7 = exp01["dg_bind_wham_kcal"]
else:
    dg_spink7 = -9.4  # placeholder
ax.bar(['SPINK7-KLK5', 'LEKTI-KLK5'], [dg_spink7, dg_bind_kcal],
       color=['steelblue', 'coral'], edgecolor='black')
ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title('EXP-03: Inhibitor Ranking Verification', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-03_ranking_verification.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-03: LEKTI-KLK5 Binding Free Energy — Results Report

**Experiment ID:** EXP-03  **Feature ID:** F-03  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL/INCONCLUSIVE]

## 1. Abstract
## 2. Introduction/Background
## 3. Hypothesis
## 4. Methods
### 4.1 System Preparation
### 4.2 Enhanced Sampling
### 4.3 Analysis
### 4.4 Deviations
## 5. Controls
## 6. Results
### 6.1 Pipeline Predictions
### 6.2 Experimental Benchmark
### 6.3 Classification
### 6.4 Ranking vs EXP-01
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

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
