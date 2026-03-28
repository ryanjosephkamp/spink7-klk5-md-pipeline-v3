# EXP-02: SPINK7-KLK12 Binding Free Energy — Implementation Guide

**Experiment ID:** EXP-02  
**Feature ID:** F-02 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** SELECTIVITY TARGET  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

> **Note:** The following is the full content of `EXP-02_spink7_klk12_dg_bind/experimental_design.md`, embedded per §24.3 requirement 1.

---

### 1. Abstract

This experiment determines the binding free energy (ΔG_bind) of the SPINK7-KLK12 complex and the selectivity ratio ΔΔG_selectivity = ΔG(KLK12) − ΔG(KLK5) using the V2 pipeline. SPINK7 inhibits KLK12 with 11.5-fold weaker affinity than KLK5 (Ki = 1500 nM vs. 132 nM, ΔΔG = +1.5 kcal/mol). This relative measurement is inherently more accurate as a benchmark than absolute ΔG values because systematic computational errors partially cancel in the difference. The experiment tests whether the pipeline can discriminate between cognate and non-cognate protease targets.

### 2. Hypothesis

**H₁:** The V2 pipeline will predict that SPINK7 binds KLK12 more weakly than KLK5, with ΔΔG_selectivity (KLK12 − KLK5) within the 95% CI [−1.2, +4.2] kcal/mol.

**H₂:** The absolute ΔG_bind for SPINK7-KLK12 will fall within the range [−12.4, −3.4] kcal/mol (from −7.9 ± 1.96 × 2.3 kcal/mol).

These hypotheses are falsifiable: if the pipeline predicts KLK12 binding to be tighter than KLK5 (ΔΔG < 0 by >2.7 kcal/mol), both hypotheses are rejected.

### 3. Background and Rationale

#### 3.1 Scientific Context

Azouz et al. (2020) demonstrated that SPINK7 is not a promiscuous serine protease inhibitor. While it potently inhibits KLK5 (Ki = 132 nM), it shows 11.5-fold weaker inhibition of KLK12 (Ki = 1500 nM) and no detectable inhibition of KLK7, KLK11, or KLK13. The ΔΔG_selectivity = +1.5 kcal/mol captures a physically meaningful difference — likely arising from subsite incompatibilities between the SPINK7 RSL and KLK12 active site subsites compared to KLK5.

#### 3.2 Pipeline Capability

The pipeline can compute ΔG_bind for SPINK7-KLK12 using the same SMD/Jarzynski and US/WHAM methodology as EXP-01, with a homology model or docking-derived starting structure of the SPINK7-KLK12 complex. The selectivity ΔΔG is computed as the difference: ΔΔG = ΔG(SPINK7-KLK12) − ΔG(SPINK7-KLK5).

#### 3.3 What This Reveals

Selectivity validation tests the pipeline's ability to capture subtle energetic differences between closely related protease targets sharing the same inhibitor. A PASS validates the pipeline's sensitivity to subsite-level interactions; a FAIL (predicting equal or reversed selectivity) indicates the pipeline cannot resolve ~1.5 kcal/mol differences — a significant limitation for drug discovery applications.

### 4. Experimental Protocol

#### 4.1 System Preparation

##### 4.1.1 Structure Acquisition

| Parameter | Value |
|-----------|-------|
| SPINK7 structure | PDB 2LEO (NMR ensemble, chain A, model 1) |
| KLK12 structure | Homology model based on KLK5 (PDB 2PSX) or AlphaFold prediction |
| Complex model | ClusPro or equivalent docking of SPINK7 onto KLK12 |
| Note | No experimental co-crystal structure exists for SPINK7-KLK12 |

##### 4.1.2 Protonation and Solvation

| Parameter | Value | Source |
|-----------|-------|--------|
| Force field | AMBER ff14SB (`amber14-all.xml`) | `config.py: SystemConfig.force_field` |
| Water model | TIP3P (`amber14/tip3p.xml`) | `config.py: SystemConfig.water_model` |
| pH | 7.4 | `config.py: SystemConfig.ph` |
| Box shape | Cubic | `config.py: SystemConfig.box_shape` |
| Box padding | 1.2 nm | `config.py: SystemConfig.box_padding_nm` |
| Ionic strength | 0.15 M NaCl | `config.py: SystemConfig.ionic_strength_molar` |

##### 4.1.3 Energy Minimization

| Parameter | Value | Source |
|-----------|-------|--------|
| Max iterations | 10,000 | `config.py: MinimizationConfig.max_iterations` |
| Tolerance | 10.0 kJ/mol/nm | `config.py: MinimizationConfig.tolerance_kj_mol_nm` |

##### 4.1.4 Equilibration

| Parameter | Value | Source |
|-----------|-------|--------|
| NVT duration | 500 ps | `config.py: EquilibrationConfig.nvt_duration_ps` |
| NPT duration | 1000 ps | `config.py: EquilibrationConfig.npt_duration_ps` |
| Temperature | 310 K | `config.py: EquilibrationConfig.temperature_k` |
| Friction coefficient | 1.0 ps⁻¹ | `config.py: EquilibrationConfig.friction_per_ps` |
| Timestep | 0.002 ps | `config.py: EquilibrationConfig.timestep_ps` |
| Pressure | 1.0 atm | `config.py: EquilibrationConfig.pressure_atm` |
| Barostat interval | 25 steps | `config.py: EquilibrationConfig.barostat_interval` |
| Restraint force constant | 1000 kJ/mol/nm² | `config.py: EquilibrationConfig.restraint_k_kj_mol_nm2` |
| Save interval | 10 ps | `config.py: EquilibrationConfig.save_interval_ps` |

#### 4.2 Umbrella Sampling

| Parameter | Value | Source |
|-----------|-------|--------|
| ξ range | 1.5–4.0 nm | `config.py: UmbrellaConfig` |
| Window spacing | 0.05 nm (51 windows) | `config.py: UmbrellaConfig.window_spacing_nm` |
| Spring constant | 1000 kJ/mol/nm² | `config.py: UmbrellaConfig.spring_constant_kj_mol_nm2` |
| Per-window duration | 10.0 ns | `config.py: UmbrellaConfig.per_window_duration_ns` |
| Pre-positioning velocity | 0.01 nm/ps | `config.py: UmbrellaConfig.pre_position_velocity_nm_per_ps` |
| Equilibration per window | 200 ps | `config.py: UmbrellaConfig.equilibration_duration_ps` |
| Detect equilibration | Enabled | `config.py: UmbrellaConfig.detect_equilibration` |
| Save interval | 1.0 ps | `config.py: UmbrellaConfig.save_interval_ps` |
| Total simulation time | 51 × 10 ns = 510 ns |

#### 4.3 SMD

| Parameter | Value | Source |
|-----------|-------|--------|
| Spring constant | 1000 kJ/mol/nm² | `config.py: SMDConfig.spring_constant_kj_mol_nm2` |
| Pulling velocity | 0.001 nm/ps | `config.py: SMDConfig.pulling_velocity_nm_per_ps` |
| Pull distance | 3.0 nm | `config.py: SMDConfig.pull_distance_nm` |
| Replicates | 50 | `config.py: SMDConfig.n_replicates` |
| Save interval | 1.0 ps | `config.py: SMDConfig.save_interval_ps` |

#### 4.4 WHAM and MBAR Analysis

Parameters identical to EXP-01: WHAM tolerance 10⁻⁶, max iterations 100,000, bootstraps 200, bins 200. MBAR solver "robust", tolerance 10⁻⁷, max iterations 10,000, bootstraps 200, bins 200.

#### 4.5 Selectivity Calculation

$$\Delta\Delta G_{\text{selectivity}} = \Delta G_{\text{bind}}^{\text{KLK12}} - \Delta G_{\text{bind}}^{\text{KLK5}}$$

Uncertainty propagation:

$$\sigma_{\Delta\Delta G} = \sqrt{\sigma_{\Delta G_{\text{KLK12}}}^2 + \sigma_{\Delta G_{\text{KLK5}}}^2}$$

### 5. Control Conditions

#### 5.1 Positive Control

**EXP-01 (SPINK7-KLK5):** The companion experiment. Both must be completed to compute ΔΔG_selectivity.

#### 5.2 Negative Control / Sanity Checks

1. KLK12 structure should be stable during equilibration (backbone RMSD < 5 Å).
2. The SPINK7-KLK12 complex should maintain a canonical binding geometry.
3. ΔG_bind(KLK12) should be negative.

### 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind (SPINK7-KLK12) | −7.9 kcal/mol | [−12.4, −3.4] |
| ΔΔG_selectivity | +1.5 kcal/mol | [−1.2, +4.2] |
| Direction of selectivity | KLK12 weaker than KLK5 | Must be positive |

### 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation Implied | Severity |
|-------------|--------------|-------------------|----------|
| **KLK12 homology model error** | Complex dissociates or non-canonical binding | Model quality insufficient | High |
| **Inability to resolve 1.5 kcal/mol** | ΔΔG ≈ 0 or wrong sign | Insufficient precision | High |
| **Systematic cancellation failure** | ΔΔG grossly wrong | Error cancellation invalid | Medium |
| **KLK12 active site geometry** | Different binding mode | Docking failed | Medium |

### 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After KLK12 model acquisition | Structural alignment to KLK5 | Cα RMSD < 3 Å for active site |
| After docking | SPINK7 RSL in KLK12 active site | P1 within 5 Å of catalytic Ser |
| After minimization | IV-1 satisfied | Energy decreased |
| After equilibration | IV-2, IV-3, IV-4, IV-6, IV-7 satisfied | All pass |
| After production MD | Complex remains intact | No dissociation |
| After umbrella sampling | IV-8 satisfied | ≥10% overlap |
| After WHAM | IV-9 satisfied | Converged |
| Selectivity comparison | Both EXP-01 and EXP-02 done | ΔΔG computable |

---

## Part 2 — Step-by-Step Implementation Instructions

> **Execution context:** All code runs in the project virtual environment with `medium_project_2/` as the working directory. All imports reference the V2 pipeline as-is (§21). Parameter overrides are via runtime arguments only. This experiment depends on EXP-01 results for the selectivity (ΔΔG) comparison.

### Step 1: Environment and Directory Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path
from dataclasses import replace

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-02_spink7_klk12_dg_bind"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import (
    SystemConfig, MinimizationConfig, EquilibrationConfig,
    ProductionConfig, SMDConfig, UmbrellaConfig, WHAMConfig,
    MBARConfig, BOLTZMANN_KJ, KCAL_TO_KJ
)
from src.prep.pdb_fetch import fetch_pdb, fetch_alphafold
from src.prep.pdb_clean import clean_structure
from src.prep.protonate import assign_protonation
from src.prep.topology import build_topology
from src.prep.solvate import solvate_system
from src.simulate.minimizer import minimize_energy
from src.simulate.equilibrate import run_nvt, run_npt
from src.simulate.production import run_production
from src.simulate.smd import run_smd_campaign
from src.simulate.umbrella import (
    run_umbrella_campaign, generate_window_centers,
    diagnose_histogram_coverage, compute_overlap_matrix
)
from src.simulate.platform import select_platform
from src.analyze.wham import solve_wham, bootstrap_pmf_uncertainty
from src.analyze.mbar import solve_mbar, bootstrap_mbar_uncertainty
from src.analyze.jarzynski import (
    jarzynski_free_energy, diagnose_dissipation, evaluate_convergence
)
from src.analyze.structural import compute_rmsd
from src.analyze.trajectory import load_trajectory
from src.physics.units import kj_to_kcal, nm_to_angstrom, kbt
from src.physics.collective_variables import com_distance
from src.visualization.plot_pmf import plot_pmf
from src.visualization.plot_timeseries import plot_rmsd_timeseries
import matplotlib.pyplot as plt

print("All imports successful.")
```

**Expected output:** All imports succeed without error.

**If this step fails:** Ensure virtual environment is activated and all pipeline dependencies are installed.

---

### Step 2: KLK12 Structure Acquisition

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

# Option A: Try AlphaFold for KLK12
# KLK12 UniProt accession: Q9UKR3
try:
    klk12_pdb = fetch_alphafold("Q9UKR3", data_dir)
    print(f"KLK12 AlphaFold model: {klk12_pdb}")
except Exception as e:
    print(f"AlphaFold fetch failed: {e}")
    # Option B: Homology model from KLK5 — requires manual construction
    klk12_pdb = data_dir / "klk12_homology_model.pdb"
    if not klk12_pdb.exists():
        print("WARNING: No KLK12 structure available. Must build homology model manually.")

# Fetch SPINK7
spink7_pdb = fetch_pdb("2LEO", data_dir)
spink7_clean = clean_structure(
    pdb_path=spink7_pdb,
    chains_to_keep=["A"],
    remove_heteroatoms=True,
    remove_waters=True,
    model_index=1,
)

# Clean KLK12
klk12_clean = clean_structure(
    pdb_path=klk12_pdb,
    chains_to_keep=["A"],
    remove_heteroatoms=True,
    remove_waters=True,
    model_index=1,
)
```

**Verification:**
```python
assert spink7_clean.exists(), "SPINK7 cleaned PDB not found"
assert klk12_clean.exists(), "KLK12 cleaned PDB not found"

# Verify KLK12 structural alignment to KLK5
import mdtraj as md
klk12_traj = md.load(str(klk12_clean))
print(f"KLK12: {klk12_traj.n_atoms} atoms, {klk12_traj.n_residues} residues")
```

**If this step fails:** If AlphaFold model is unavailable, construct a homology model using Swiss-Model or similar tool. If no KLK12 model can be obtained, classify experiment as INCONCLUSIVE (requires external structure).

---

### Step 3: Complex Assembly and Protonation

```python
# Dock SPINK7 onto KLK12 using canonical Kazal binding mode
# Place the docked complex PDB at:
docked_complex_path = data_dir / "spink7_klk12_complex.pdb"

if not docked_complex_path.exists():
    print("WARNING: Pre-docked SPINK7-KLK12 complex not found.")
    print("Required: dock SPINK7 onto KLK12 analogously to KLK5 docking.")
    print(f"Place at: {docked_complex_path}")

# Assign protonation states
complex_protonated = assign_protonation(
    pdb_path=docked_complex_path,
    ph=7.4,
    force_field="AMBER",
    use_propka=True,
)
```

**Verification:**
```python
assert complex_protonated.exists(), "Protonated complex not found"
traj = md.load(str(complex_protonated))
print(f"Complex: {traj.n_atoms} atoms, chains: {traj.n_chains}")
```

---

### Step 4: Topology Building, Solvation, and Minimization

```python
from openmm.app import PME, Simulation, PDBFile
from openmm import LangevinMiddleIntegrator, XmlSerializer
import openmm

sys_config = SystemConfig()
topology, system, modeller = build_topology(
    pdb_path=complex_protonated,
    system_config=sys_config,
    nonbonded_method=PME,
    nonbonded_cutoff_nm=1.0,
)

modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)
print(f"Solvated: {n_waters} waters, {n_pos} Na+, {n_neg} Cl-")

solvated_pdb = OUTPUT_DIR / "solvated_complex.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

# Minimization
platform = select_platform()
integrator = LangevinMiddleIntegrator(
    310 * openmm.unit.kelvin, 1.0 / openmm.unit.picosecond, 0.002 * openmm.unit.picoseconds
)
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

e_initial = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
    openmm.unit.kilojoules_per_mole
)
min_results = minimize_energy(sim, MinimizationConfig())
e_final = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
    openmm.unit.kilojoules_per_mole
)

assert e_final < e_initial, f"IV-1 FAIL: E_min ({e_final}) >= E_initial ({e_initial})"
print(f"IV-1 PASS: Energy decreased by {e_initial - e_final:.2f} kJ/mol")
```

---

### Step 5: Equilibration (NVT + NPT)

```python
eq_config = EquilibrationConfig()
eq_output = OUTPUT_DIR / "equilibration"
eq_output.mkdir(exist_ok=True)

nvt_results = run_nvt(sim, eq_config, eq_output)
t_avg = nvt_results.get("average_temperature_k", 0)
assert abs(t_avg - 310.0) < 5.0, f"IV-2 FAIL: |T_avg - 310| = {abs(t_avg - 310.0):.2f} K"
print(f"IV-2 PASS: T_avg = {t_avg:.2f} K")

npt_results = run_npt(sim, eq_config, eq_output)
rho = npt_results.get("average_density_g_per_cm3", 0)
assert 0.95 <= rho <= 1.05, f"IV-3 FAIL: density = {rho:.4f}"
print(f"IV-3 PASS: density = {rho:.4f} g/cm³")

eq_state_path = eq_output / "equilibrated_state.xml"
sim.saveState(str(eq_state_path))
system_xml_path = eq_output / "system.xml"
with open(system_xml_path, "w") as f:
    f.write(XmlSerializer.serialize(system))
```

**Verification — IV-4 (structural stability):**
```python
eq_traj = load_trajectory(eq_output / "npt_trajectory.dcd", solvated_pdb)
ref_traj = md.load(str(solvated_pdb))
rmsd_eq = compute_rmsd(eq_traj, ref_traj, atom_selection="backbone")
assert np.max(rmsd_eq) < 0.5, f"IV-4 FAIL: max RMSD = {nm_to_angstrom(np.max(rmsd_eq)):.2f} Å"
print(f"IV-4 PASS: max RMSD = {nm_to_angstrom(np.max(rmsd_eq)):.2f} Å")
```

---

### Step 6: Production MD (100 ns)

```python
prod_config = ProductionConfig()
prod_output = OUTPUT_DIR / "production"
prod_output.mkdir(exist_ok=True)

prod_results = run_production(sim, prod_config, prod_output)
energy_drift = prod_results.get("energy_drift_kj_mol_ns_atom", 0)
assert abs(energy_drift) < 0.1, f"IV-5 FAIL: drift = {energy_drift}"
print(f"IV-5 PASS: energy drift = {energy_drift:.6f}")
```

---

### Step 7: Identify Pull Groups and Run SMD

```python
prod_traj = load_trajectory(prod_output / "production_trajectory.dcd", solvated_pdb, stride=100)
topology_md = prod_traj.topology

chain_a_indices = topology_md.select("chainid 0")  # SPINK7
chain_b_indices = topology_md.select("chainid 1")  # KLK12

pull_group_1 = chain_a_indices.tolist()
pull_group_2 = chain_b_indices.tolist()
print(f"Pull group 1 (SPINK7): {len(pull_group_1)} atoms")
print(f"Pull group 2 (KLK12): {len(pull_group_2)} atoms")

smd_config = SMDConfig()
smd_output = OUTPUT_DIR / "smd"
smd_output.mkdir(exist_ok=True)

smd_results = run_smd_campaign(
    equilibrated_state_path=eq_state_path,
    system_xml_path=system_xml_path,
    config=smd_config,
    pull_group_1=pull_group_1,
    pull_group_2=pull_group_2,
    output_dir=smd_output,
    topology_pdb_path=solvated_pdb,
)

work_values = np.array([r["total_work_kj_mol"] for r in smd_results])
assert np.all(np.isfinite(work_values)), "Non-finite work values"
print(f"SMD complete: {len(smd_results)} replicates, mean work = {work_values.mean():.2f} kJ/mol")
```

---

### Step 8: SMD Analysis

```python
temperature_k = 310.0
jarz_results = jarzynski_free_energy(work_values, temperature_k)
dg_jarz_kj = jarz_results["delta_g_kj_mol"]
dg_jarz_kcal = kj_to_kcal(dg_jarz_kj)
print(f"Jarzynski ΔG: {dg_jarz_kcal:.2f} kcal/mol")

dissipation = diagnose_dissipation(work_values, temperature_k)
conv_results = evaluate_convergence(work_values, temperature_k, n_subsets=10)
```

---

### Step 9: Umbrella Sampling Campaign

```python
us_config = UmbrellaConfig()
us_output = OUTPUT_DIR / "umbrella"
us_output.mkdir(exist_ok=True)

window_centers = generate_window_centers(us_config)
us_results = run_umbrella_campaign(
    equilibrated_state_path=eq_state_path,
    system_xml_path=system_xml_path,
    config=us_config,
    pull_group_1=pull_group_1,
    pull_group_2=pull_group_2,
    output_dir=us_output,
    pdb_path=solvated_pdb,
    topology_pdb_path=solvated_pdb,
)

xi_timeseries_list = []
for i in range(len(us_results)):
    xi_ts = np.load(us_output / f"window_{i:03d}" / "xi_timeseries.npy")
    xi_timeseries_list.append(xi_ts)

overlap = compute_overlap_matrix(xi_timeseries_list)
min_overlap = min(overlap[i, i+1] for i in range(len(overlap)-1))
assert min_overlap >= 0.10, f"IV-8 FAIL: min overlap = {min_overlap:.4f}"
print(f"IV-8 PASS: min adjacent overlap = {min_overlap:.4f}")
```

---

### Step 10: WHAM and MBAR Analysis

```python
wham_config = WHAMConfig()
mbar_config = MBARConfig()
spring_constants = np.full(len(window_centers), us_config.spring_constant_kj_mol_nm2)

wham_results = solve_wham(xi_timeseries_list, window_centers, spring_constants, temperature_k, wham_config)
assert wham_results["converged"], "IV-9 FAIL: WHAM did not converge"
print("IV-9 PASS: WHAM converged")

wham_bootstrap = bootstrap_pmf_uncertainty(xi_timeseries_list, window_centers, spring_constants, temperature_k, wham_config)

mbar_results = solve_mbar(xi_timeseries_list, window_centers, spring_constants, temperature_k, mbar_config)
mbar_bootstrap = bootstrap_mbar_uncertainty(xi_timeseries_list, window_centers, spring_constants, temperature_k, mbar_config)

pmf_wham = wham_results["pmf_kj_mol"]
pmf_mbar = mbar_results["pmf_kj_mol"]
xi_bins = wham_results["xi_bins_nm"]
pmf_std = wham_bootstrap["pmf_std_kj_mol"]
pmf_mbar_std = mbar_bootstrap["pmf_std_kj_mol"]
```

---

### Step 11: Binding Free Energy Extraction and Selectivity

```python
kT = kbt(temperature_k)
beta = 1.0 / kT
C_standard_nm3 = 1.0 / (1660.0 * 0.001)

bound_mask = xi_bins < 2.0
bulk_mask = xi_bins > 3.5

pmf_shifted = pmf_wham - pmf_wham[bulk_mask].min()
I_site = np.trapz(np.exp(-beta * pmf_shifted[bound_mask]) * xi_bins[bound_mask]**2, xi_bins[bound_mask])
I_bulk = np.trapz(np.exp(-beta * pmf_shifted[bulk_mask]) * xi_bins[bulk_mask]**2, xi_bins[bulk_mask])

dg_bind_klk12_kj = -kT * np.log(C_standard_nm3 / (4*np.pi) * I_site) + kT * np.log(C_standard_nm3 / (4*np.pi) * I_bulk)
dg_bind_klk12_kcal = kj_to_kcal(dg_bind_klk12_kj)
print(f"ΔG_bind (SPINK7-KLK12, WHAM): {dg_bind_klk12_kcal:.2f} kcal/mol")

# Load EXP-01 results for selectivity comparison
exp01_results_path = PROJECT_ROOT / "v3_experiments" / "EXP-01_spink7_klk5_dg_bind" / "results.json"
if exp01_results_path.exists():
    with open(exp01_results_path) as f:
        exp01_results = json.load(f)
    dg_bind_klk5_kcal = exp01_results["dg_bind_wham_kcal"]
    
    ddg_selectivity = dg_bind_klk12_kcal - dg_bind_klk5_kcal
    print(f"ΔΔG_selectivity = {ddg_selectivity:.2f} kcal/mol")
    print(f"Expected: +1.5 kcal/mol, 95% CI: [-1.2, +4.2]")
    
    if -1.2 <= ddg_selectivity <= 4.2:
        selectivity_class = "PASS"
    elif ddg_selectivity > 0:
        selectivity_class = "MARGINAL"
    else:
        selectivity_class = "FAIL"
    print(f"Selectivity classification: {selectivity_class}")
else:
    print("WARNING: EXP-01 results not yet available. Run EXP-01 first for selectivity comparison.")
    ddg_selectivity = None
```

---

### Step 12: Statistical Classification

```python
dg_exp = -7.9
sigma_exp = 0.6
sigma_comp = kj_to_kcal(np.mean(pmf_std))
sigma_method = 2.0
sigma_combined = np.sqrt(sigma_exp**2 + sigma_comp**2 + sigma_method**2)

discrepancy = abs(dg_bind_klk12_kcal - dg_exp)
z = 1.96

if sigma_comp > sigma_exp:
    classification = "INCONCLUSIVE"
elif discrepancy <= z * sigma_combined:
    classification = "PASS"
elif discrepancy <= 2 * z * sigma_combined:
    classification = "MARGINAL"
else:
    classification = "FAIL"

print(f"\n=== CLASSIFICATION ===")
print(f"Pipeline prediction (WHAM): {dg_bind_klk12_kcal:.2f} kcal/mol")
print(f"Experimental benchmark: {dg_exp:.2f} ± {sigma_exp:.2f} kcal/mol")
print(f"Classification: {classification}")

results = {
    "experiment_id": "EXP-02",
    "feature_id": "F-02",
    "dg_bind_klk12_wham_kcal": float(dg_bind_klk12_kcal),
    "dg_jarz_kcal": float(dg_jarz_kcal),
    "dg_exp_kcal": float(dg_exp),
    "ddg_selectivity_kcal": float(ddg_selectivity) if ddg_selectivity is not None else None,
    "sigma_combined": float(sigma_combined),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Selectivity Schematic

**Description:** Thermodynamic diagram comparing SPINK7-KLK5 and SPINK7-KLK12 binding energetics, showing the ΔΔG_selectivity.

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Energy level diagram
y_klk5 = -9.4
y_klk12 = -7.9
y_unbound = 0.0

ax.hlines(y_unbound, 0.5, 3.5, colors='black', linewidth=2)
ax.hlines(y_klk5, 0.5, 1.5, colors='blue', linewidth=3, label='SPINK7-KLK5')
ax.hlines(y_klk12, 2.5, 3.5, colors='red', linewidth=3, label='SPINK7-KLK12')

ax.annotate('', xy=(1.0, y_klk5), xytext=(1.0, y_unbound),
            arrowprops=dict(arrowstyle='<->', color='blue', lw=2))
ax.text(0.6, (y_klk5 + y_unbound)/2, f'ΔG = {y_klk5:.1f}\nkcal/mol', fontsize=11, color='blue')

ax.annotate('', xy=(3.0, y_klk12), xytext=(3.0, y_unbound),
            arrowprops=dict(arrowstyle='<->', color='red', lw=2))
ax.text(3.1, (y_klk12 + y_unbound)/2, f'ΔG = {y_klk12:.1f}\nkcal/mol', fontsize=11, color='red')

ax.annotate('', xy=(2.0, y_klk12), xytext=(2.0, y_klk5),
            arrowprops=dict(arrowstyle='<->', color='green', lw=2))
ax.text(2.05, (y_klk5 + y_klk12)/2, f'ΔΔG = +1.5\nkcal/mol', fontsize=11, color='green', fontweight='bold')

ax.set_ylabel('Free Energy (kcal/mol)', fontsize=13)
ax.set_title('EXP-02: SPINK7 Selectivity — KLK5 vs KLK12', fontsize=14)
ax.legend(fontsize=12)
ax.set_xlim(0, 4)
ax.set_xticks([])
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-02_selectivity_schematic.png", dpi=300)
plt.close(fig)
```

### Figure 2: PMF Comparison (KLK12)

**Description:** WHAM and MBAR PMF profiles for SPINK7-KLK12 with uncertainty bands.

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
pmf_wham_kcal = pmf_wham / KCAL_TO_KJ
pmf_std_kcal = pmf_std / KCAL_TO_KJ
pmf_mbar_kcal = pmf_mbar / KCAL_TO_KJ

ax.plot(xi_bins, pmf_wham_kcal, 'b-', linewidth=2, label='WHAM')
ax.fill_between(xi_bins, pmf_wham_kcal - 1.96*pmf_std_kcal,
                pmf_wham_kcal + 1.96*pmf_std_kcal, alpha=0.2, color='blue')
ax.plot(xi_bins, pmf_mbar_kcal, 'r--', linewidth=2, label='MBAR')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-02: SPINK7-KLK12 PMF Profile', fontsize=14)
ax.legend(fontsize=12)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-02_pmf_profile.png", dpi=300)
plt.close(fig)
```

### Figure 3: Predicted vs Experimental (Bar Chart)

**Description:** Comparison of predicted ΔG values for both KLK5 and KLK12 targets against experimental benchmarks.

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
methods = ['KLK5 (exp)', 'KLK5 (pred)', 'KLK12 (exp)', 'KLK12 (pred)']
values = [-9.4, dg_bind_klk5_kcal if ddg_selectivity is not None else 0,
          -7.9, dg_bind_klk12_kcal]
colors = ['gold', 'steelblue', 'lightsalmon', 'salmon']

ax.bar(methods, values, color=colors, edgecolor='black', alpha=0.85)
ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title('EXP-02: Selectivity — Predicted vs Experimental', fontsize=14)
plt.xticks(rotation=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-02_predicted_vs_experimental.png", dpi=300)
plt.close(fig)
```

### Figure 4: SMD Work Distribution

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 5))
work_kcal = work_values / KCAL_TO_KJ
ax.hist(work_kcal, bins=20, edgecolor='black', alpha=0.7, color='steelblue')
ax.axvline(dg_jarz_kcal, color='red', linestyle='--', linewidth=2,
           label=f'Jarzynski ΔG = {dg_jarz_kcal:.1f} kcal/mol')
ax.set_xlabel('Work (kcal/mol)', fontsize=13)
ax.set_ylabel('Count', fontsize=13)
ax.set_title('EXP-02: SMD Work Distribution (KLK12)', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-02_smd_work_distribution.png", dpi=300)
plt.close(fig)
```

### Figure 5: Umbrella Window Overlap

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 5))
for xi_ts in xi_timeseries_list:
    ax.hist(xi_ts, bins=50, alpha=0.3, density=True)
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('Probability density', fontsize=13)
ax.set_title('EXP-02: Umbrella Window Overlap (KLK12)', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-02_umbrella_overlap.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-02: SPINK7-KLK12 Binding Free Energy — Results Report

**Experiment ID:** EXP-02
**Feature ID:** F-02
**Date:** [execution date]
**Classification:** [PASS / MARGINAL / FAIL / INCONCLUSIVE]

## 1. Abstract
[One-paragraph summary including ΔG_bind and ΔΔG_selectivity results]

## 2. Introduction/Background
[SPINK7-KLK12 context; Azouz et al. 2020 selectivity data]

## 3. Hypothesis
[Restate H₁, H₂]

## 4. Methods
### 4.1 System Preparation
[KLK12 model source, docking, protonation, solvation — all explicit]
### 4.2 Enhanced Sampling
[SMD and US parameters as executed]
### 4.3 Analysis
[WHAM, MBAR, Jarzynski parameters]
### 4.4 Deviations
[Any deviations or "None"]

## 5. Controls
[EXP-01 results; sanity checks]

## 6. Results
### 6.1 Absolute ΔG_bind
| Method | ΔG_bind (kcal/mol) | 95% CI |
|--------|-------------------|---------|
| US/WHAM | X.XX ± Y.YY | [A, B] |
| SMD/Jarzynski | X.XX ± Y.YY | [A, B] |

### 6.2 Selectivity
ΔΔG_selectivity = X.XX ± Y.YY kcal/mol
Expected: +1.5 kcal/mol, 95% CI [-1.2, +4.2]

### 6.3 Classification
**Classification: [PASS / MARGINAL / FAIL / INCONCLUSIVE]**

## 7. Discussion
[Interpretation; selectivity implications]

## 8. Conclusions
[Summary]

## 9. Figures
[≥ 5 embedded figures]

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
