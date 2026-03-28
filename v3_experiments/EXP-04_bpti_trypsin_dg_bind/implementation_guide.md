# EXP-04: BPTI-Trypsin Binding Free Energy — Implementation Guide

**Experiment ID:** EXP-04  
**Feature ID:** F-04 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** GOLD STANDARD REFERENCE  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

This experiment determines the binding free energy of the BPTI-trypsin complex, the gold-standard protease-inhibitor system (Kd = 6 × 10⁻¹⁴ M, ΔG ≈ −18 kcal/mol), using US/WHAM and SMD/Jarzynski methods. Unlike EXP-01, this system benefits from a high-resolution co-crystal structure (PDB 2PTC, 1.9 Å), removing model quality as a confounding variable. This is the definitive test of the pipeline's free energy methodology for protease-inhibitor systems — success here validates the method; failure indicates fundamental methodological limitations independent of model quality.

### 2. Hypothesis

**H₁:** The V2 pipeline's US/WHAM estimate of ΔG_bind for BPTI-trypsin will fall within the 95% CI [−24.7, −11.3] kcal/mol.

**H₂:** The wide CI reflects the extreme difficulty of computationally reproducing ultra-tight binding (−18 kcal/mol). The pipeline may achieve MARGINAL classification due to the fundamental challenge of sampling the deep free energy basin required for femtomolar affinity.

### 3. Background and Rationale

#### 3.1 Scientific Context

BPTI (58-residue Kunitz-type inhibitor) and trypsin form the most extensively characterized protease-inhibitor complex. Two independent labs confirm Kd ≈ 6 × 10⁻¹⁴ M (Vincent & Lazdunski 1972, Castro & Anderson 1996). The co-crystal structure (PDB 2PTC at 1.9 Å) provides an ideal starting point. While BPTI is Kunitz-type (not Kazal), it shares the canonical serine protease inhibition mechanism with SPINK7.

#### 3.2 What This Reveals

As the gold standard with no model uncertainty, this experiment isolates method accuracy from model quality. A PASS validates the computational methodology; a FAIL attributable to sampling or force field is a fundamental pipeline limitation.

### 4. Experimental Protocol

#### 4.1 System Preparation

| Parameter | Value |
|-----------|-------|
| Complex structure | PDB 2PTC (BPTI-trypsin co-crystal, 1.9 Å) |
| Force field | AMBER ff14SB (`amber14-all.xml`) |
| Water model | TIP3P (`amber14/tip3p.xml`) |
| pH | 7.4 |
| Box padding | 1.2 nm |
| Ionic strength | 0.15 M NaCl |
| Box shape | Cubic |

#### 4.2 Minimization and Equilibration

Minimization 10,000 steps, tolerance 10.0 kJ/mol/nm. NVT 500 ps at 310 K, friction 1.0 ps⁻¹, timestep 2 fs. NPT 1000 ps at 310 K and 1.0 atm, barostat interval 25. Restraint 1000 kJ/mol/nm², save interval 10 ps.

#### 4.3 Production MD

Duration 100 ns, temperature 310 K, timestep 0.002 ps, save interval 10 ps.

#### 4.4 SMD

Spring constant 1000 kJ/mol/nm², pulling velocity 0.001 nm/ps, pull distance 3.0 nm, 50 replicates. Total 150 ns.

#### 4.5 Umbrella Sampling

ξ range 1.5–4.0 nm (51 windows at 0.05 nm spacing), spring constant 1000 kJ/mol/nm², per-window duration 10.0 ns, pre-positioning velocity 0.01 nm/ps, equilibration per window 200 ps, save interval 1.0 ps. Total 510 ns.

#### 4.6 WHAM and MBAR

WHAM: tolerance 10⁻⁶, max iterations 100,000, bootstraps 200, bins 200. MBAR: solver "robust", tolerance 10⁻⁷, max iterations 10,000, bootstraps 200, bins 200.

#### 4.7 Statistical Comparison

σ_exp = 0.5 kcal/mol, σ_method = 3.0 kcal/mol (ultra-tight binding), σ_combined ≈ 3.4 kcal/mol. 95% CI [−24.7, −11.3].

### 5. Control Conditions

**Positive Control:** EXP-29 (SH3-p41) — method validation where three computational approaches agree within 0.2 kcal/mol of experiment.

**Negative Controls:** All invariants IV-1 through IV-10 satisfied. BPTI backbone RMSD < 1.5 Å from 2PTC crystal. K15 (P1) remains in S1 pocket. Three disulfide bonds (Cys5-Cys55, Cys14-Cys38, Cys30-Cys51) intact.

### 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind (US/WHAM) | −18.0 kcal/mol | [−24.7, −11.3] |
| PMF minimum depth | Deep, well-defined | > 10 kcal/mol below bulk |

Classification: PASS within [−24.7, −11.3], MARGINAL within [−31.4, −4.6], FAIL outside.

### 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation Implied | Severity |
|-------------|--------------|-------------------|----------|
| ΔG insufficiently negative | ΔG > −11 kcal/mol | Cannot capture ultra-tight binding | High |
| Incomplete dissociation in SMD | Plateau not reached in 3 nm pull | Pull distance insufficient | Medium |
| WHAM convergence issues | Deep bound-state histogram undersampling | Need more windows | Medium |

### 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After structure loading | 2PTC correctly parsed; both chains present | BPTI + trypsin resolved |
| After minimization | IV-1 satisfied | Energy decreased |
| After equilibration | All invariants pass; K15 in S1 pocket | Stable complex |
| After production MD | RMSD < 2 Å; complex intact | Structural stability |
| After SMD | Work distributions unimodal (IV-10) | Clean pulling |
| After US | Histogram overlap ≥ 10% (IV-8) | Sufficient sampling |
| After WHAM | Convergence < 10⁻⁶ (IV-9) | WHAM converged |

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment and Directory Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import (
    SystemConfig, MinimizationConfig, EquilibrationConfig,
    ProductionConfig, SMDConfig, UmbrellaConfig, WHAMConfig, MBARConfig, KCAL_TO_KJ
)
from src.prep.pdb_fetch import fetch_pdb
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
from src.analyze.jarzynski import jarzynski_free_energy, diagnose_dissipation
from src.analyze.convergence import evaluate_convergence
from src.analyze.structural import compute_rmsd, compute_rmsf
from src.analyze.trajectory import load_trajectory, align_trajectory
from src.physics.units import kj_to_kcal, nm_to_angstrom, kbt
from src.visualization.plot_pmf import plot_pmf
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

---

### Step 2: Structure Acquisition — PDB 2PTC Co-Crystal

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

# BPTI-trypsin co-crystal — NO docking needed, co-crystal structure
complex_pdb = fetch_pdb("2PTC", data_dir)
complex_clean = clean_structure(
    complex_pdb, chains_to_keep=["E", "I"],  # E=trypsin, I=BPTI
    remove_heteroatoms=True, remove_waters=True, model_index=1
)
assert complex_clean.exists(), "2PTC cleaned PDB not found"

# Verify both chains present
traj = md.load(str(complex_clean))
chains = list(traj.topology.chains)
assert len(chains) >= 2, f"Expected 2 chains, found {len(chains)}"
n_residues = traj.topology.n_residues
print(f"2PTC loaded: {len(chains)} chains, {n_residues} residues, {traj.topology.n_atoms} atoms")

# Verify BPTI disulfide bonds (Cys5-Cys55, Cys14-Cys38, Cys30-Cys51)
bpti_chain = chains[1]  # I chain = BPTI
print("BPTI chain residues:", bpti_chain.n_residues)
```

**Expected output:** 2PTC loaded with 2 chains (~280 trypsin residues + 58 BPTI residues).

**If this step fails:** Verify PDB ID and chain identifiers. 2PTC is well-established and should always be available.

---

### Step 3: Protonation and System Setup

```python
from openmm.app import PME, Simulation, PDBFile
from openmm import LangevinMiddleIntegrator, XmlSerializer
import openmm

complex_protonated = assign_protonation(complex_clean, ph=7.4, force_field="AMBER", use_propka=True)

sys_config = SystemConfig()
topology, system, modeller = build_topology(
    complex_protonated, sys_config,
    nonbonded_method=PME, nonbonded_cutoff_nm=1.0
)
modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)
print(f"Solvated: {n_waters} waters, {n_pos} Na+, {n_neg} Cl-")

solvated_pdb = OUTPUT_DIR / "solvated_complex.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
```

---

### Step 4: Minimization (IV-1)

```python
platform = select_platform()
integrator = LangevinMiddleIntegrator(
    310 * openmm.unit.kelvin,
    1.0 / openmm.unit.picosecond,
    0.002 * openmm.unit.picoseconds
)
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

e_initial = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
    openmm.unit.kilojoules_per_mole
)
min_config = MinimizationConfig()
minimize_energy(sim, min_config)
e_final = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
    openmm.unit.kilojoules_per_mole
)

assert e_final < e_initial, f"IV-1 FAIL: E_final={e_final:.1f} >= E_initial={e_initial:.1f}"
print(f"IV-1 PASS: E decreased {e_initial:.1f} → {e_final:.1f} kJ/mol (Δ = {e_initial - e_final:.1f})")
```

---

### Step 5: Equilibration (IV-2, IV-3, IV-4)

```python
eq_config = EquilibrationConfig()
eq_output = OUTPUT_DIR / "equilibration"
eq_output.mkdir(exist_ok=True)

# NVT — temperature equilibration
nvt_results = run_nvt(sim, eq_config, eq_output)
t_avg = nvt_results.get("average_temperature_k", 0)
assert abs(t_avg - 310.0) < 5.0, f"IV-2 FAIL: T_avg = {t_avg:.1f} K"
print(f"IV-2 PASS: T_avg = {t_avg:.1f} K (target 310 K)")

# NPT — pressure/density equilibration
npt_results = run_npt(sim, eq_config, eq_output)
rho = npt_results.get("average_density_g_per_cm3", 0)
assert 0.95 <= rho <= 1.05, f"IV-3 FAIL: density = {rho:.3f} g/cm³"
print(f"IV-3 PASS: density = {rho:.3f} g/cm³")

# Save equilibrated state
eq_state_path = eq_output / "equilibrated_state.xml"
sim.saveState(str(eq_state_path))
system_xml_path = eq_output / "system.xml"
with open(system_xml_path, "w") as f:
    f.write(XmlSerializer.serialize(system))

# IV-4: RMSD < 3 Å from starting structure
eq_traj = load_trajectory(eq_output / "npt_trajectory.dcd", solvated_pdb, stride=10)
rmsd = compute_rmsd(eq_traj, eq_traj[0], atom_sel="backbone")
assert rmsd[-1] < 0.3, f"IV-4 FAIL: final RMSD = {rmsd[-1]*10:.1f} Å"
print(f"IV-4 PASS: final backbone RMSD = {rmsd[-1]*10:.1f} Å")
```

---

### Step 6: Production MD (100 ns)

```python
prod_config = ProductionConfig()
prod_output = OUTPUT_DIR / "production"
prod_output.mkdir(exist_ok=True)
prod_results = run_production(sim, prod_config, prod_output)

# IV-5: Complex stability
prod_traj = load_trajectory(
    prod_output / "production_trajectory.dcd", solvated_pdb, stride=100
)
rmsd_prod = compute_rmsd(prod_traj, prod_traj[0], atom_sel="backbone")
assert np.max(rmsd_prod) < 0.3, f"IV-5 FAIL: max RMSD = {np.max(rmsd_prod)*10:.1f} Å"
print(f"IV-5 PASS: max backbone RMSD = {np.max(rmsd_prod)*10:.1f} Å")

# Verify K15 (P1 residue) in S1 pocket — check distance
topology_md = prod_traj.topology
bpti_ca = topology_md.select("chainid 1 and name CA")
trypsin_ca = topology_md.select("chainid 0 and name CA")
```

---

### Step 7: Pull Group Identification

```python
chain_bpti_indices = topology_md.select("chainid 1").tolist()   # BPTI
chain_trypsin_indices = topology_md.select("chainid 0").tolist()  # Trypsin

print(f"Pull group 1 (BPTI): {len(chain_bpti_indices)} atoms")
print(f"Pull group 2 (Trypsin): {len(chain_trypsin_indices)} atoms")
```

---

### Step 8: SMD Campaign (50 replicates)

```python
smd_config = SMDConfig()
smd_output = OUTPUT_DIR / "smd"
smd_output.mkdir(exist_ok=True)

smd_results = run_smd_campaign(
    equilibrated_state_path=eq_state_path,
    system_xml_path=system_xml_path,
    config=smd_config,
    pull_group_1=chain_bpti_indices,
    pull_group_2=chain_trypsin_indices,
    output_dir=smd_output,
    pdb_path=solvated_pdb,
    topology_pdb_path=solvated_pdb,
)

work_values = np.array([r["work_kj_mol"] for r in smd_results])
print(f"SMD complete: {len(work_values)} replicates")
print(f"Work values: mean={np.mean(work_values):.1f}, std={np.std(work_values):.1f} kJ/mol")

# IV-10: Dissipation diagnostic
dissipation = diagnose_dissipation(work_values, 310.0, threshold=3.0)
print(f"Dissipation: {dissipation}")

# Jarzynski free energy
jarz_results = jarzynski_free_energy(work_values, 310.0)
dg_jarz_kj = jarz_results["delta_g_kj_mol"]
dg_jarz_kcal = kj_to_kcal(dg_jarz_kj)
print(f"ΔG (Jarzynski) = {dg_jarz_kcal:.2f} kcal/mol")
```

---

### Step 9: Umbrella Sampling Campaign

```python
us_config = UmbrellaConfig()
us_output = OUTPUT_DIR / "umbrella"
us_output.mkdir(exist_ok=True)
window_centers = generate_window_centers(us_config)
print(f"Umbrella windows: {len(window_centers)} at centers {window_centers[0]:.2f}–{window_centers[-1]:.2f} nm")

us_results = run_umbrella_campaign(
    equilibrated_state_path=eq_state_path,
    system_xml_path=system_xml_path,
    config=us_config,
    pull_group_1=chain_bpti_indices,
    pull_group_2=chain_trypsin_indices,
    output_dir=us_output,
    pdb_path=solvated_pdb,
    topology_pdb_path=solvated_pdb,
)

xi_timeseries_list = [
    np.load(us_output / f"window_{i:03d}" / "xi_timeseries.npy")
    for i in range(len(us_results))
]

# IV-8: Overlap matrix
overlap = compute_overlap_matrix(xi_timeseries_list)
min_overlap = min(overlap[i, i + 1] for i in range(len(overlap) - 1))
assert min_overlap >= 0.10, f"IV-8 FAIL: min overlap = {min_overlap:.3f}"
print(f"IV-8 PASS: min adjacent overlap = {min_overlap:.3f}")
```

---

### Step 10: WHAM and MBAR Analysis

```python
wham_config = WHAMConfig()
spring_constants = np.full(len(window_centers), us_config.spring_constant_kj_mol_nm2)

wham_results = solve_wham(
    xi_timeseries_list, window_centers, spring_constants, 310.0, wham_config
)
assert wham_results["converged"], "IV-9 FAIL: WHAM did not converge"
print("IV-9 PASS: WHAM converged")

wham_bootstrap = bootstrap_pmf_uncertainty(
    xi_timeseries_list, window_centers, spring_constants, 310.0, wham_config
)

mbar_config = MBARConfig()
mbar_results = solve_mbar(
    xi_timeseries_list, window_centers, spring_constants, 310.0, mbar_config
)
mbar_bootstrap = bootstrap_mbar_uncertainty(
    xi_timeseries_list, window_centers, spring_constants, 310.0, mbar_config
)
```

---

### Step 11: ΔG Extraction and Standard-State Correction

```python
pmf_wham = wham_results["pmf_kj_mol"]
xi_bins = wham_results["xi_bins_nm"]
pmf_std = wham_bootstrap["pmf_std_kj_mol"]

kT = kbt(310.0)
beta = 1.0 / kT
C_standard_nm3 = 1.0 / (1660.0 * 0.001)  # 1 M in nm⁻³

# Identify bound and bulk regions
bound_mask = xi_bins < 2.0
bulk_mask = xi_bins > 3.5
pmf_shifted = pmf_wham - pmf_wham[bulk_mask].min()

# Volumetric integrals
I_site = np.trapz(np.exp(-beta * pmf_shifted[bound_mask]) * xi_bins[bound_mask]**2,
                   xi_bins[bound_mask])
I_bulk = np.trapz(np.exp(-beta * pmf_shifted[bulk_mask]) * xi_bins[bulk_mask]**2,
                   xi_bins[bulk_mask])

dg_bind_kj = (-kT * np.log(C_standard_nm3 / (4 * np.pi) * I_site)
               + kT * np.log(C_standard_nm3 / (4 * np.pi) * I_bulk))
dg_bind_kcal = kj_to_kcal(dg_bind_kj)
print(f"ΔG_bind (WHAM) = {dg_bind_kcal:.2f} kcal/mol")

# MBAR ΔG
pmf_mbar = mbar_results["pmf_kj_mol"]
pmf_mbar_shifted = pmf_mbar - pmf_mbar[bulk_mask].min()
I_site_mbar = np.trapz(np.exp(-beta * pmf_mbar_shifted[bound_mask]) * xi_bins[bound_mask]**2,
                        xi_bins[bound_mask])
dg_mbar_kj = (-kT * np.log(C_standard_nm3 / (4 * np.pi) * I_site_mbar)
               + kT * np.log(C_standard_nm3 / (4 * np.pi) * I_bulk))
dg_mbar_kcal = kj_to_kcal(dg_mbar_kj)
print(f"ΔG_bind (MBAR) = {dg_mbar_kcal:.2f} kcal/mol")
```

---

### Step 12: Classification (§25.1)

```python
dg_exp = -18.0
sigma_exp = 0.5
sigma_comp = kj_to_kcal(np.mean(pmf_std))
sigma_method = 3.0  # ultra-tight binding
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
    "experiment_id": "EXP-04",
    "feature_id": "F-04",
    "system": "BPTI-trypsin (PDB 2PTC)",
    "dg_bind_wham_kcal": float(dg_bind_kcal),
    "dg_bind_mbar_kcal": float(dg_mbar_kcal),
    "dg_jarzynski_kcal": float(dg_jarz_kcal),
    "dg_exp_kcal": float(dg_exp),
    "sigma_exp": float(sigma_exp),
    "sigma_comp": float(sigma_comp),
    "sigma_method": float(sigma_method),
    "sigma_combined": float(sigma_combined),
    "discrepancy": float(discrepancy),
    "classification": classification,
    "ci_95_lower": float(dg_exp - z * sigma_combined),
    "ci_95_upper": float(dg_exp + z * sigma_combined),
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\n{'='*60}")
print(f"EXP-04 BPTI-Trypsin — Classification: {classification}")
print(f"ΔG_bind (WHAM)  = {dg_bind_kcal:.2f} kcal/mol")
print(f"ΔG_bind (MBAR)  = {dg_mbar_kcal:.2f} kcal/mol")
print(f"ΔG_bind (Jarz)  = {dg_jarz_kcal:.2f} kcal/mol")
print(f"ΔG_exp          = {dg_exp:.2f} kcal/mol")
print(f"σ_combined      = {sigma_combined:.2f} kcal/mol")
print(f"95% CI          = [{dg_exp - z*sigma_combined:.1f}, {dg_exp + z*sigma_combined:.1f}]")
print(f"{'='*60}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Experimental Design — Gold Standard Schematic

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: System overview
ax = axes[0]
ax.text(0.5, 0.8, 'PDB 2PTC', ha='center', fontsize=14, fontweight='bold')
ax.text(0.5, 0.6, 'Co-crystal 1.9 Å', ha='center', fontsize=11)
ax.text(0.5, 0.4, 'Trypsin + BPTI', ha='center', fontsize=11)
ax.text(0.5, 0.2, 'K15 → S1 pocket', ha='center', fontsize=10, style='italic')
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
ax.set_title('A) System', fontsize=12)

# Panel B: Energetics
ax = axes[1]
ax.barh(['Kd', 'ΔG'], [6e-14, -18], color=['coral', 'steelblue'])
ax.set_title('B) Experimental Values', fontsize=12)
ax.set_xlabel('Value')

# Panel C: Pipeline workflow
ax = axes[2]
steps = ['2PTC\n(co-crystal)', 'Equilibrate\n(100 ns)', 'SMD\n(50 reps)', 'US/WHAM\n(51 windows)', 'ΔG\n(classify)']
y_pos = range(len(steps))
ax.barh(y_pos, [1]*5, color='lightsteelblue', edgecolor='navy')
for i, s in enumerate(steps):
    ax.text(0.5, i, s, ha='center', va='center', fontsize=9)
ax.set_xlim(0, 1); ax.set_yticks([]); ax.axis('off')
ax.set_title('C) Pipeline', fontsize=12)

plt.suptitle('EXP-04: BPTI-Trypsin Gold Standard', fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-04_design_schematic.png", dpi=300)
plt.close(fig)
```

### Figure 2: PMF Profile (WHAM + MBAR)

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
pmf_kcal = pmf_wham / KCAL_TO_KJ
std_kcal = pmf_std / KCAL_TO_KJ
pmf_mbar_kcal = pmf_mbar / KCAL_TO_KJ

ax.plot(xi_bins, pmf_kcal, 'b-', linewidth=2, label='WHAM')
ax.fill_between(xi_bins, pmf_kcal - 1.96*std_kcal, pmf_kcal + 1.96*std_kcal,
                alpha=0.15, color='blue')
ax.plot(xi_bins, pmf_mbar_kcal, 'r--', linewidth=2, label='MBAR')
ax.axhline(y=dg_exp, color='gold', linestyle=':', linewidth=2, label=f'Experimental ({dg_exp} kcal/mol)')
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-04: BPTI-Trypsin PMF Profile', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-04_pmf_profile.png", dpi=300)
plt.close(fig)
```

### Figure 3: SMD Work Distribution

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
work_kcal = work_values / KCAL_TO_KJ
ax.hist(work_kcal, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
ax.axvline(x=dg_jarz_kcal, color='red', linestyle='--', linewidth=2,
           label=f'Jarzynski ΔG = {dg_jarz_kcal:.1f} kcal/mol')
ax.axvline(x=dg_exp, color='gold', linestyle=':', linewidth=2,
           label=f'Experimental = {dg_exp} kcal/mol')
ax.set_xlabel('Work (kcal/mol)', fontsize=13)
ax.set_ylabel('Count', fontsize=13)
ax.set_title('EXP-04: SMD Work Distribution (50 replicates)', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-04_smd_work_distribution.png", dpi=300)
plt.close(fig)
```

### Figure 4: Predicted vs Experimental Comparison

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
methods = ['Experimental', 'WHAM', 'MBAR', 'Jarzynski']
values = [dg_exp, dg_bind_kcal, dg_mbar_kcal, dg_jarz_kcal]
colors = ['gold', 'steelblue', 'coral', 'mediumseagreen']
bars = ax.bar(methods, values, color=colors, edgecolor='black')
ax.axhspan(dg_exp - z*sigma_combined, dg_exp + z*sigma_combined,
           alpha=0.1, color='green', label='95% CI')
ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title(f'EXP-04: BPTI-Trypsin — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-04_predicted_vs_experimental.png", dpi=300)
plt.close(fig)
```

### Figure 5: RMSD Time Series

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
time_ns = np.arange(len(rmsd_prod)) * 0.1  # stride=100, save=10ps
ax.plot(time_ns, rmsd_prod * 10, 'b-', linewidth=0.8)
ax.axhline(y=1.5, color='red', linestyle='--', label='Threshold (1.5 Å)')
ax.set_xlabel('Time (ns)', fontsize=13)
ax.set_ylabel('Backbone RMSD (Å)', fontsize=13)
ax.set_title('EXP-04: BPTI-Trypsin Production RMSD', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-04_rmsd_timeseries.png", dpi=300)
plt.close(fig)
```

### Figure 6: Umbrella Window Overlap Histogram

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 5))
for xi_ts in xi_timeseries_list:
    ax.hist(xi_ts, bins=50, alpha=0.3, density=True)
ax.set_xlabel(r'$\xi$ (nm)', fontsize=13)
ax.set_ylabel('Probability density', fontsize=13)
ax.set_title('EXP-04: Umbrella Window Overlap', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-04_umbrella_overlap.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-04: BPTI-Trypsin Binding Free Energy — Results Report

**Experiment ID:** EXP-04  **Feature ID:** F-04  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL/INCONCLUSIVE]

## 1. Abstract
## 2. Introduction/Background
## 3. Hypothesis
## 4. Methods
### 4.1 System Preparation (PDB 2PTC co-crystal)
### 4.2 Equilibration and Production
### 4.3 SMD Campaign
### 4.4 Umbrella Sampling
### 4.5 WHAM/MBAR Analysis
### 4.6 Deviations from Protocol
## 5. Controls
### 5.1 Disulfide Bond Integrity
### 5.2 K15-S1 Contact
## 6. Results
### 6.1 ΔG_bind (WHAM): [value] ± [σ] kcal/mol
### 6.2 ΔG_bind (MBAR): [value] ± [σ] kcal/mol
### 6.3 ΔG_bind (Jarzynski): [value] ± [σ] kcal/mol
### 6.4 Experimental ΔG: −18.0 ± 0.5 kcal/mol
### 6.5 σ_combined: [value]
### 6.6 Classification: [PASS/MARGINAL/FAIL/INCONCLUSIVE]
## 7. Discussion
## 8. Conclusions
## 9. Figures
## 10. References
- Vincent & Lazdunski (1972)
- Castro & Anderson (1996)
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

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps — they do not replace them. All simulation code in Part 2 should use the GPU platform and checkpoint utilities defined below when running on Colab.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
# ─── Run this cell FIRST in every Colab session ───

# 1. Verify GPU availability
!nvidia-smi

# 2. Mount Google Drive for persistent storage
from google.colab import drive
drive.mount('/content/drive')

# 3. Install dependencies
!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

# 4. Create working directories
import os, sys
from pathlib import Path

EXP_ID = "EXP-04"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
DRIVE_BASE.mkdir(parents=True, exist_ok=True)
(DRIVE_BASE / "checkpoints").mkdir(exist_ok=True)
(DRIVE_BASE / "outputs").mkdir(exist_ok=True)
(DRIVE_BASE / "figures").mkdir(exist_ok=True)
(DRIVE_BASE / "outputs" / "structures").mkdir(exist_ok=True)
(DRIVE_BASE / "outputs" / "equilibration").mkdir(exist_ok=True)
(DRIVE_BASE / "outputs" / "production").mkdir(exist_ok=True)
(DRIVE_BASE / "outputs" / "smd").mkdir(exist_ok=True)
(DRIVE_BASE / "outputs" / "umbrella").mkdir(exist_ok=True)

# 5. Upload/clone project source
# Option A: git clone (if repo is accessible)
# !git clone https://github.com/ryanjosephkamp/medium_project_2.git /content/project
# Option B: symlink from Drive
!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

# 6. Override paths for Colab
PROJECT_ROOT = Path("/content")
EXP_DIR = DRIVE_BASE
OUTPUT_DIR = DRIVE_BASE / "outputs"
FIGURES_DIR = DRIVE_BASE / "figures"

print(f"Colab environment ready. Drive base: {DRIVE_BASE}")
```

### §5.2 GPU Platform Selection

```python
# §5.2 GPU Platform Selection
# ─── Replace select_platform() in Part 2 Step 4 with explicit CUDA ───

import openmm
from src.simulate.platform import select_platform

# Method A: Use pipeline function with explicit CUDA request
platform = select_platform("CUDA")

# Method B: Direct OpenMM CUDA selection with properties
platform = openmm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}

# Verify GPU is active
print(f"Platform: {platform.getName()}")
print(f"CUDA devices: {platform.getPropertyDefaultValue('DeviceIndex')}")

# IMPORTANT: Update ALL Simulation() calls in Part 2 to include properties:
#   Step 4 (Minimization):
#     sim = Simulation(modeller.topology, system, integrator, platform, properties)
#   SMD campaign (Step 8): run_smd_campaign uses select_platform() internally —
#     ensure src/simulate/smd.py calls select_platform("CUDA")
#   US campaign (Step 9): run_umbrella_campaign uses select_platform() internally —
#     ensure src/simulate/umbrella.py calls select_platform("CUDA")
#
# The select_platform("CUDA") call propagates through the pipeline automatically
# since platform.py probes CUDA first in auto-detect mode and Colab provides CUDA.
```

### §5.3 Checkpoint and Resume Integration

```python
# §5.3 Checkpoint and Resume Integration
# ─── Enables recovery from Colab 24-hour session timeouts ───

import json, time
from pathlib import Path
from openmm import XmlSerializer
import openmm

CHECKPOINT_DIR = DRIVE_BASE / "checkpoints"

def save_checkpoint(simulation, phase_name, step=None, extra_data=None):
    """Save simulation state to Drive for session recovery.

    Args:
        simulation: OpenMM Simulation object.
        phase_name: str — e.g., 'minimize', 'nvt', 'npt', 'production_10ns'.
        step: int or None — current step number.
        extra_data: dict or None — additional metadata to save.
    """
    ckpt_path = CHECKPOINT_DIR / f"{phase_name}.chk"
    simulation.saveCheckpoint(str(ckpt_path))

    state = simulation.context.getState(
        getPositions=True, getVelocities=True, getEnergy=True
    )
    with open(CHECKPOINT_DIR / f"{phase_name}_state.xml", 'w') as f:
        f.write(XmlSerializer.serialize(state))

    progress = {
        "phase": phase_name,
        "step": step,
        "energy_kj": state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoules_per_mole
        ),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    if extra_data:
        progress.update(extra_data)
    with open(CHECKPOINT_DIR / "progress.json", 'w') as f:
        json.dump(progress, f, indent=2)
    print(f"  Checkpoint saved: {phase_name} (step={step})")


def load_checkpoint(simulation, phase_name):
    """Resume simulation from Drive checkpoint if it exists.

    Returns True if checkpoint was loaded, False otherwise.
    Verifies energy consistency after loading.
    """
    ckpt_path = CHECKPOINT_DIR / f"{phase_name}.chk"
    if not ckpt_path.exists():
        return False

    simulation.loadCheckpoint(str(ckpt_path))
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(
        openmm.unit.kilojoules_per_mole
    )
    print(f"  Resumed from checkpoint: {phase_name} (E={energy:.1f} kJ/mol)")
    return True


# ─── Checkpoint insertion points for EXP-04 ───
# After Step 4 (Minimization):
#   save_checkpoint(sim, "minimize")
#
# After Step 5 NVT:
#   save_checkpoint(sim, "nvt")
#
# After Step 5 NPT:
#   save_checkpoint(sim, "npt")
#
# During Step 6 (Production, every 10 ns):
#   for chunk in range(10):
#       sim.step(5_000_000)  # 10 ns at 2 fs
#       save_checkpoint(sim, f"production_{(chunk+1)*10}ns", step=(chunk+1)*5_000_000)
#
# During Step 8 (SMD, after each replicate):
#   # run_smd_campaign handles internally; verify outputs exist on Drive
#   save_checkpoint(sim, f"smd_rep_{rep_idx:03d}")
#
# During Step 9 (US, after each window):
#   save_checkpoint(sim, f"us_window_{win_idx:03d}")
```

### §5.4 Progress Monitoring

```python
# §5.4 Progress Monitoring
# ─── Track GPU utilization and simulation throughput ───

import time, subprocess

def report_gpu_status():
    """Print GPU memory usage and utilization."""
    result = subprocess.run(
        ['nvidia-smi', '--query-gpu=name,memory.used,memory.total,utilization.gpu',
         '--format=csv,noheader'],
        capture_output=True, text=True
    )
    print(f"GPU: {result.stdout.strip()}")


def estimate_performance(simulation, n_steps=5000, dt_ps=0.002):
    """Measure simulation throughput in ns/day.

    Args:
        simulation: OpenMM Simulation object.
        n_steps: Number of steps for benchmarking.
        dt_ps: Timestep in picoseconds.

    Returns:
        float: Performance in ns/day.
    """
    t0 = time.time()
    simulation.step(n_steps)
    elapsed = time.time() - t0
    ns_simulated = n_steps * dt_ps / 1000.0
    ns_per_day = ns_simulated / elapsed * 86400
    print(f"  Performance: {ns_per_day:.1f} ns/day ({elapsed:.1f}s for {ns_simulated:.4f} ns)")
    report_gpu_status()
    return ns_per_day


# ─── EXP-04 runtime estimates (A100 40GB, mixed precision) ───
# Production (100 ns):  ~2–4 hours
# SMD (50 reps × 3 ns): ~3–6 hours
# US (51 windows × 10 ns): ~10–20 hours
# Total estimated wall-clock: ~16–30 hours (within single Colab session)
#
# If approaching 24-hour limit:
#   1. Save checkpoint for current phase
#   2. Reconnect to new session
#   3. Resume from last checkpoint
```

### §5.5 Results Persistence

```python
# §5.5 Results Persistence
# ─── All outputs saved to Google Drive ───

import shutil

def sync_to_drive(local_dir, drive_subdir="outputs"):
    """Copy local output files to Google Drive.

    Args:
        local_dir: Path — local directory to sync.
        drive_subdir: str — subdirectory under DRIVE_BASE.
    """
    drive_dir = DRIVE_BASE / drive_subdir
    drive_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for f in Path(local_dir).rglob("*"):
        if f.is_file():
            dest = drive_dir / f.relative_to(local_dir)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dest)
            count += 1
    print(f"  Synced {count} files: {local_dir} → {drive_dir}")


def verify_drive_sync():
    """Verify critical output files exist on Drive."""
    critical_files = [
        "results.json",
        "checkpoints/progress.json",
        "outputs/umbrella/window_000/xi_timeseries.npy",
    ]
    print("Drive sync verification:")
    for f in critical_files:
        path = DRIVE_BASE / f
        if path.exists():
            print(f"  ✓ {f} ({path.stat().st_size:,} bytes)")
        else:
            print(f"  ✗ {f} MISSING")


# ─── Persistence schedule for EXP-04 ───
# After each major phase, call sync_to_drive():
#   sync_to_drive(OUTPUT_DIR / "structures", "outputs/structures")
#   sync_to_drive(OUTPUT_DIR / "equilibration", "outputs/equilibration")
#   sync_to_drive(OUTPUT_DIR / "production", "outputs/production")
#   sync_to_drive(OUTPUT_DIR / "smd", "outputs/smd")
#   sync_to_drive(OUTPUT_DIR / "umbrella", "outputs/umbrella")
#
# Final results:
#   shutil.copy2(EXP_DIR / "results.json", DRIVE_BASE / "results.json")
```

### §5.6 GPU-Optimized Figure Generation

```python
# §5.6 GPU-Optimized Figure Generation
# ─── Save figures to both local and Drive paths ───

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for Colab
import matplotlib.pyplot as plt

def save_figure(fig, name):
    """Save figure to Drive figures directory.

    Args:
        fig: matplotlib Figure object.
        name: str — filename without extension.
    """
    fig_path = DRIVE_BASE / "figures" / f"{name}.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(fig_path), dpi=150, bbox_inches='tight')
    print(f"  Figure saved: {fig_path}")
    plt.close(fig)

# Usage: replace all fig.savefig() calls in Part 3 with:
#   save_figure(fig, "EXP-04_pmf_profile")
#   save_figure(fig, "EXP-04_smd_work_distribution")
#   save_figure(fig, "EXP-04_predicted_vs_experimental")
#   save_figure(fig, "EXP-04_rmsd_timeseries")
#   save_figure(fig, "EXP-04_umbrella_overlap")
#   save_figure(fig, "EXP-04_design_schematic")

# GPU performance figure (additional)
def plot_gpu_performance(ns_per_day_values, phase_labels):
    """Generate GPU throughput summary figure."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(phase_labels, ns_per_day_values, color='steelblue', edgecolor='black')
    ax.set_xlabel('Throughput (ns/day)')
    ax.set_title(f'EXP-04: GPU Performance Summary (A100)')
    plt.tight_layout()
    save_figure(fig, "EXP-04_gpu_performance")
```

### §5.7 Error Recovery Procedures

```python
# §5.7 Error Recovery Procedures
# ─── Handle Colab-specific failure modes ───

def gpu_safe_run(func, *args, max_retries=2, **kwargs):
    """Execute a simulation function with GPU error recovery.

    Handles:
      - CUDA out-of-memory: reduces precision to 'single'
      - NaN energy: re-minimizes from last checkpoint
      - Colab timeout: saves checkpoint before exit

    Args:
        func: callable — simulation function to execute.
        max_retries: int — number of retry attempts.

    Returns:
        Result of func(*args, **kwargs).
    """
    for attempt in range(max_retries + 1):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            error_msg = str(e)
            # Log error to Drive
            with open(DRIVE_BASE / "error_log.txt", 'a') as f:
                f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} "
                        f"[attempt {attempt+1}] {type(e).__name__}: {error_msg}\n")

            if "out of memory" in error_msg.lower() or "CUDA" in error_msg:
                print(f"  GPU OOM (attempt {attempt+1}/{max_retries+1})")
                print("  Recovery: switching to single precision")
                properties['CudaPrecision'] = 'single'
                if attempt < max_retries:
                    continue

            elif "nan" in error_msg.lower():
                print(f"  NaN detected (attempt {attempt+1}/{max_retries+1})")
                print("  Recovery: reload last checkpoint and re-minimize")
                if attempt < max_retries:
                    continue

            elif "timeout" in error_msg.lower():
                print("  Colab timeout approaching — saving checkpoint")
                # Checkpoint should already be saved by periodic saves

            raise  # Re-raise if all retries exhausted


# ─── EXP-04–specific error scenarios ───
#
# 1. SMD force explosion:
#    - Symptom: NaN in work values during run_smd_campaign
#    - Fix: Reduce pulling velocity from 0.001 to 0.0005 nm/ps
#    - Resume: Re-run failed replicates only
#
# 2. US window histogram gap:
#    - Symptom: IV-8 failure (overlap < 10%)
#    - Fix: Insert additional windows at gap location
#    - Diagnostic: diagnose_histogram_coverage(xi_timeseries_list)
#
# 3. WHAM non-convergence:
#    - Symptom: IV-9 failure (tolerance not reached)
#    - Fix: Increase max_iterations or relax tolerance
#    - Verify: Check adjacent window overlap matrix
#
# 4. Colab session disconnect during production MD:
#    - Resume: load_checkpoint(sim, "production_XXns")
#    - Continue from last 10ns checkpoint
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
