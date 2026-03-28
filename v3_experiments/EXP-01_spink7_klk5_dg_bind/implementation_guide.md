# EXP-01: SPINK7-KLK5 Binding Free Energy — Implementation Guide

**Experiment ID:** EXP-01  
**Feature ID:** F-01 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** PRIMARY PIPELINE TARGET  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

> **Note:** The following is the full content of `EXP-01_spink7_klk5_dg_bind/experimental_design.md`, embedded per §24.3 requirement 1.

---

### 1. Abstract

This experiment determines the binding free energy (ΔG_bind) of the SPINK7-KLK5 protease-antiprotease complex using the V2 pipeline's two complementary enhanced sampling strategies: Steered Molecular Dynamics (SMD) with the Jarzynski equality and Umbrella Sampling (US) with WHAM reconstruction of the Potential of Mean Force (PMF). The experimentally measured Ki = 132 nM (ΔG ≈ −9.4 kcal/mol) from Azouz et al. (2020) serves as the primary benchmark. This is the single most important validation target for the entire V3 benchmarking project, as it directly tests whether the pipeline can predict the binding thermodynamics of its primary target system. Cross-validation using the Bennett Acceptance Ratio (BAR) estimator and MBAR will provide internal consistency checks across free energy estimation methods.

### 2. Hypothesis

**H₁:** The V2 pipeline's Umbrella Sampling + WHAM estimate of ΔG_bind for the SPINK7-KLK5 complex will fall within the 95% confidence interval [−14.0, −4.8] kcal/mol of the experimentally measured value (ΔG = −9.4 kcal/mol, derived from Ki = 132 nM via ΔG = RT ln(Ki)).

**H₂:** The V2 pipeline's SMD + Jarzynski estimate of ΔG_bind will fall within the 95% CI [−14.0, −4.8] kcal/mol, with the second-order cumulant expansion providing a more precise estimate than the exponential average.

**H₃:** All three free energy estimators (Jarzynski exponential average, cumulant expansion, and BAR) applied to the SMD work distributions will agree within 2 kcal/mol, and the US/WHAM and MBAR estimates from umbrella sampling data will agree within 1 kcal/mol.

These hypotheses are falsifiable: if any estimate falls outside the stated CI, the hypothesis is rejected and a pipeline limitation is identified.

### 3. Background and Rationale

#### 3.1 Scientific Context

SPINK7 (Serine Peptidase Inhibitor, Kazal Type 7) is a ~6 kDa secreted Kazal-type inhibitor that stoichiometrically inhibits KLK5, a trypsin-like serine protease. In Eosinophilic Esophagitis (EoE), IL-13-mediated transcriptional silencing of SPINK7 unleashes KLK5 proteolytic activity, degrading Desmoglein-1 and compromising the epithelial barrier. The SPINK7-KLK5 Ki = 132 nM was measured using the Morrison tight-binding inhibitor equation from dose-response curves (Azouz et al. 2020, Fig. 1B), with three independent experiments performed in duplicate.

#### 3.2 Pipeline Capability

The V2 pipeline provides two complementary enhanced sampling pathways for computing ΔG_bind:

1. **SMD + Jarzynski:** Non-equilibrium pulling of SPINK7 away from KLK5 along the COM-COM reaction coordinate ξ, with free energy extracted via the Jarzynski equality (§5.3). The pipeline implements the exponential average, second-order cumulant expansion, and BAR estimators (`src/analyze/jarzynski.py`).

2. **Umbrella Sampling + WHAM:** Equilibrium sampling of ξ in $M$ biased windows with harmonic restraints, reconstructing the PMF via WHAM (`src/analyze/wham.py`) and MBAR (`src/analyze/mbar.py`). The binding free energy is extracted from the PMF using the standard concentration correction (§5.5).

#### 3.3 What This Reveals

This experiment tests the pipeline's core capability — predicting protein-protein binding free energies from first principles. Since no co-crystal structure of SPINK7-KLK5 exists, the starting structure relies on ClusPro docking (Azouz et al. 2020, Fig. 1C), adding model uncertainty. A PASS result validates both the computational methodology and the homology-based complex model; a FAIL result identifies specific limitations (force field, sampling, or model quality) for Part 2 remediation.

### 4. Experimental Protocol

#### 4.1 System Preparation

##### 4.1.1 Structure Acquisition

| Parameter | Value |
|-----------|-------|
| SPINK7 structure | PDB 2LEO (NMR ensemble, chain A, model 1) |
| KLK5 structure | PDB 2PSX (X-ray, chain A) |
| Complex model | ClusPro docking of 2LEO onto 2PSX (canonical binding mode from Azouz et al. 2020) |
| Missing residues | Fill using PDBFixer (`src/prep/pdb_clean.py`) |
| Disulfide bonds | Verify and enforce all three SPINK7 disulfide bonds (IV-6) |

##### 4.1.2 Protonation and Solvation

| Parameter | Value | Source |
|-----------|-------|--------|
| Force field | AMBER ff14SB (`amber14-all.xml`) | `config.py: SystemConfig.force_field` |
| Water model | TIP3P (`amber14/tip3p.xml`) | `config.py: SystemConfig.water_model` |
| pH | 7.4 | `config.py: SystemConfig.ph` |
| Box shape | Cubic | `config.py: SystemConfig.box_shape` |
| Box padding | 1.2 nm | `config.py: SystemConfig.box_padding_nm` |
| Ionic strength | 0.15 M NaCl | `config.py: SystemConfig.ionic_strength_molar` |
| Positive ion | Na⁺ | `config.py: SystemConfig.positive_ion` |
| Negative ion | Cl⁻ | `config.py: SystemConfig.negative_ion` |

##### 4.1.3 Energy Minimization

| Parameter | Value | Source |
|-----------|-------|--------|
| Max iterations | 10,000 | `config.py: MinimizationConfig.max_iterations` |
| Tolerance | 10.0 kJ/mol/nm | `config.py: MinimizationConfig.tolerance_kj_mol_nm` |
| Invariant check | IV-1: $E_{\text{min}} < E_{\text{initial}}$ | §6 |

##### 4.1.4 Equilibration

| Parameter | Value | Source |
|-----------|-------|--------|
| NVT duration | 500 ps | `config.py: EquilibrationConfig.nvt_duration_ps` |
| NPT duration | 1000 ps | `config.py: EquilibrationConfig.npt_duration_ps` |
| Temperature | 310 K | `config.py: EquilibrationConfig.temperature_k` |
| Friction coefficient | 1.0 ps⁻¹ | `config.py: EquilibrationConfig.friction_per_ps` |
| Timestep | 0.002 ps (2 fs) | `config.py: EquilibrationConfig.timestep_ps` |
| Pressure | 1.0 atm | `config.py: EquilibrationConfig.pressure_atm` |
| Barostat interval | 25 steps | `config.py: EquilibrationConfig.barostat_interval` |
| Restraint force constant | 1000 kJ/mol/nm² | `config.py: EquilibrationConfig.restraint_k_kj_mol_nm2` |
| Save interval | 10 ps | `config.py: EquilibrationConfig.save_interval_ps` |
| Invariant checks | IV-2: $\|T_{\text{avg}} - 310\| < 5$ K; IV-3: $\rho \in [0.95, 1.05]$ g/cm³; IV-4: backbone RMSD < 5 Å | §6 |

#### 4.2 Production MD (Pre-Pulling Equilibrium)

| Parameter | Value | Source |
|-----------|-------|--------|
| Duration | 100 ns | `config.py: ProductionConfig.duration_ns` |
| Temperature | 310 K | `config.py: ProductionConfig.temperature_k` |
| Friction coefficient | 1.0 ps⁻¹ | `config.py: ProductionConfig.friction_per_ps` |
| Timestep | 0.002 ps (2 fs) | `config.py: ProductionConfig.timestep_ps` |
| Pressure | 1.0 atm | `config.py: ProductionConfig.pressure_atm` |
| Save interval | 10 ps | `config.py: ProductionConfig.save_interval_ps` |
| Checkpoint interval | 100 ps | `config.py: ProductionConfig.checkpoint_interval_ps` |
| Invariant check | IV-5: energy drift < 0.1 kJ/mol/ns/atom; IV-7: no periodic image artifacts | §6 |

#### 4.3 Steered Molecular Dynamics (SMD)

| Parameter | Value | Source |
|-----------|-------|--------|
| Spring constant | 1000 kJ/mol/nm² | `config.py: SMDConfig.spring_constant_kj_mol_nm2` |
| Pulling velocity | 0.001 nm/ps | `config.py: SMDConfig.pulling_velocity_nm_per_ps` |
| Pull distance | 3.0 nm | `config.py: SMDConfig.pull_distance_nm` |
| Number of replicates | 50 | `config.py: SMDConfig.n_replicates` |
| Save interval | 1.0 ps | `config.py: SMDConfig.save_interval_ps` |
| Reaction coordinate | ξ = COM-COM distance between SPINK7 and KLK5 (§5.2) |
| Pulling duration | 3.0 nm / 0.001 nm/ps = 3000 ps = 3 ns per replicate |
| Total SMD simulation time | 50 × 3 ns = 150 ns |
| Invariant check | IV-10: unimodal work distributions | §6 |

**SMD Harmonic Pulling Potential (§5.3.1):**

$$U_{\text{SMD}}(\xi, t) = \frac{k}{2} \left[ \xi(t) - \xi_0 - v \cdot t \right]^2$$

**Jarzynski Equality (§5.3.3):**

$$\Delta G = -k_B T \ln \left[ \frac{1}{N_{\text{traj}}} \sum_{j=1}^{N_{\text{traj}}} e^{-\beta W_j} \right]$$

**Second-Order Cumulant Expansion:**

$$\Delta G \approx \langle W \rangle - \frac{\beta}{2} \sigma_W^2$$

#### 4.4 Umbrella Sampling

| Parameter | Value | Source |
|-----------|-------|--------|
| ξ minimum | 1.5 nm | `config.py: UmbrellaConfig.xi_min_nm` |
| ξ maximum | 4.0 nm | `config.py: UmbrellaConfig.xi_max_nm` |
| Window spacing | 0.05 nm | `config.py: UmbrellaConfig.window_spacing_nm` |
| Number of windows | (4.0 − 1.5) / 0.05 = 51 windows |
| Spring constant | 1000 kJ/mol/nm² | `config.py: UmbrellaConfig.spring_constant_kj_mol_nm2` |
| Per-window duration | 10.0 ns | `config.py: UmbrellaConfig.per_window_duration_ns` |
| Pre-positioning velocity | 0.01 nm/ps | `config.py: UmbrellaConfig.pre_position_velocity_nm_per_ps` |
| Pre-positioning spring constant | 1000 kJ/mol/nm² | `config.py: UmbrellaConfig.pre_position_spring_constant_kj_mol_nm2` |
| Equilibration per window | 200 ps | `config.py: UmbrellaConfig.equilibration_duration_ps` |
| Automated equilibration detection | Enabled | `config.py: UmbrellaConfig.detect_equilibration` |
| Save interval | 1.0 ps | `config.py: UmbrellaConfig.save_interval_ps` |
| Total US simulation time | 51 × 10 ns = 510 ns |
| Invariant check | IV-8: ≥10% overlap between adjacent windows | §6 |

#### 4.5 WHAM Analysis

| Parameter | Value | Source |
|-----------|-------|--------|
| Tolerance | 10⁻⁶ | `config.py: WHAMConfig.tolerance` |
| Max iterations | 100,000 | `config.py: WHAMConfig.max_iterations` |
| Bootstrap resamples | 200 | `config.py: WHAMConfig.n_bootstrap` |
| Histogram bins | 200 | `config.py: WHAMConfig.histogram_bins` |
| Invariant check | IV-9: WHAM convergence $\max_i |f_i^{(n+1)} - f_i^{(n)}| < 10^{-6}$ kJ/mol | §6 |

**WHAM Equations (§5.4.2):**

$$P^{\text{unbiased}}(\xi) = \frac{\sum_{i=1}^{M} n_i \, h_i(\xi)}{\sum_{i=1}^{M} n_i \, \exp\left[ \beta \left( f_i - U_i^{\text{bias}}(\xi) \right) \right]}$$

$$e^{-\beta f_i} = \int P^{\text{unbiased}}(\xi) \, \exp\left[ -\beta \, U_i^{\text{bias}}(\xi) \right] d\xi$$

#### 4.6 MBAR Analysis

| Parameter | Value | Source |
|-----------|-------|--------|
| Solver protocol | "robust" | `config.py: MBARConfig.solver_protocol` |
| Relative tolerance | 10⁻⁷ | `config.py: MBARConfig.relative_tolerance` |
| Max iterations | 10,000 | `config.py: MBARConfig.maximum_iterations` |
| Bootstrap resamples | 200 | `config.py: MBARConfig.n_bootstrap` |
| PMF bins | 200 | `config.py: MBARConfig.n_pmf_bins` |

#### 4.7 Binding Free Energy Extraction (§5.5)

$$\Delta G_{\text{bind}}^{\circ} = -k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{site}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right] + k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{bulk}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right]$$

where $C^{\circ} = 1/1660$ Å⁻³ (standard concentration of 1 M).

#### 4.8 Statistical Comparison

1. Compute pipeline prediction: ΔG_bind ± σ_comp from bootstrap CI.
2. Construct the combined 95% CI using the §25 framework:
   - σ_exp = 0.6 kcal/mol
   - σ_comp = from pipeline bootstrap
   - σ_method = 2.0 kcal/mol (US/WHAM for protein-protein binding, §25.4)
   - σ_combined = √(σ_exp² + σ_comp² + σ_method²)
3. Classify result: PASS / MARGINAL / FAIL / INCONCLUSIVE per §25.1.

### 5. Control Conditions

#### 5.1 Positive Control

**BPTI-trypsin (EXP-04):** The gold-standard protease-inhibitor system (Kd = 6 × 10⁻¹⁴ M, ΔG ≈ −18 kcal/mol) with a high-resolution co-crystal structure (PDB 2PTC). If the pipeline fails to reproduce the BPTI-trypsin ΔG_bind within its wider CI [−24.7, −11.3], the SPINK7-KLK5 result should be interpreted with caution regardless of its own classification.

**SH3-p41 (EXP-29):** Computational methods validation target with three independent methods converging within ~0.2 kcal/mol (experimental ΔG = −7.99 kcal/mol). This verifies the PMF methodology itself.

#### 5.2 Negative Control / Sanity Checks

1. **Energy conservation:** Verify IV-5 (energy drift < 0.1 kJ/mol/ns/atom) during production MD.
2. **Structural stability:** SPINK7-KLK5 complex backbone RMSD < 5 Å during production (IV-4).
3. **Disulfide integrity:** All three SPINK7 disulfide bonds remain intact (IV-6, $d_{S-S} < 2.5$ Å).
4. **Pulling direction:** Confirm ξ increases monotonically during SMD (no re-association).
5. **PMF shape:** The PMF should show a clear minimum at the bound state (ξ ≈ 1.5–2.0 nm) and plateau at large separation (ξ > 3.5 nm).

### 6. Expected Outcomes

#### 6.1 Primary Prediction

| Method | Expected Range | Source |
|--------|---------------|--------|
| US/WHAM ΔG_bind | −14.0 to −4.8 kcal/mol (95% CI) | benchmarks.md F-01 |
| SMD/Jarzynski ΔG_bind | −14.0 to −4.8 kcal/mol (95% CI) | benchmarks.md F-01 |

#### 6.2 Classification Criteria

| Classification | US/WHAM ΔG_bind Range |
|---------------|----------------------|
| **PASS** | [−14.0, −4.8] kcal/mol |
| **MARGINAL** | [−18.7, −0.1] kcal/mol |
| **FAIL** | Outside [−18.7, −0.1] kcal/mol |
| **INCONCLUSIVE** | σ_comp > σ_exp (0.6 kcal/mol) |

#### 6.3 Internal Consistency Checks

- WHAM and MBAR estimates should agree within 1 kcal/mol.
- Jarzynski exponential average and cumulant expansion should agree within 2 kcal/mol.
- Forward and reverse cumulative ΔG estimates should converge to plateau values.

### 7. Potential Failure Modes

| Failure Mode | Expected Manifestation | Pipeline Limitation Implied | Severity |
|-------------|----------------------|---------------------------|----------|
| **Homology model inaccuracy** | ΔG too positive (weak binding); complex dissociates during equilibration | ClusPro docking model geometry does not capture correct binding pose | Critical |
| **Insufficient SMD replicates** | Jarzynski ΔG not converged; exponential average dominated by rare low-work trajectories | 50 replicates insufficient for protein-protein unbinding | Medium |
| **Inadequate US window overlap** | WHAM convergence failure (IV-9 violation) or rugged PMF with unphysical barriers | Window spacing too coarse or per-window sampling too short | Medium |
| **Force field limitation** | Systematic over/underestimation of ΔG by >5 kcal/mol despite convergence | AMBER ff14SB fixed-charge model inadequate | High |
| **Pulling speed too fast** | SMD ΔG systematically too positive (dissipative work) | v = 0.001 nm/ps too fast for near-equilibrium pulling | Medium |
| **Finite-size effects** | Systematic bias from periodic boundary conditions | Box padding insufficient for 3 nm pull distance | Medium |
| **Disulfide bond instability** | SPINK7 structure deforms during simulation (IV-6 violation) | Force field or protonation incorrect for Cys residues | Critical |

### 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After structure preparation | Visual inspection of docked complex; RSL near KLK5 active site | P1 carbonyl within 4 Å of catalytic Ser Oγ |
| After minimization | IV-1: $E_{\text{min}} < E_{\text{initial}}$ | Energy decreased |
| After NVT equilibration | IV-2: $\|T_{\text{avg}} - 310\| < 5$ K | Temperature within bounds |
| After NPT equilibration | IV-3: $\rho \in [0.95, 1.05]$ g/cm³ | Density within bounds |
| After equilibration | IV-4: backbone RMSD < 5 Å from starting structure | Complex structurally stable |
| After equilibration | IV-6: all disulfide $d_{S-S} < 2.5$ Å | Disulfides intact |
| After equilibration | IV-7: no periodic image artifacts | Box large enough |
| After production MD (10 ns check) | Stable RMSD plateau; no drift | Complex remains bound |
| After each SMD replicate | Work values finite and positive | No numerical instabilities |
| After all SMD replicates | IV-10: work distribution unimodal | No pathway bifurcation |
| After umbrella sampling | IV-8: ≥10% histogram overlap between adjacent windows | Sufficient overlap |
| After WHAM | IV-9: convergence < 10⁻⁶ kJ/mol | WHAM converged |
| After PMF reconstruction | PMF minimum at bound state, plateau at large ξ | Physically reasonable PMF |

### 9. Resource Estimates

| Component | Time | Justification |
|-----------|------|---------------|
| System preparation | ~5 min | PDB fetch + clean + protonate + solvate |
| Minimization | ~5 min | 10,000 steps |
| Equilibration (NVT + NPT) | ~30 min | 1.5 ns total |
| Production MD | ~10 hours | 100 ns at ~10 ns/hr (GPU) |
| SMD (50 replicates × 3 ns) | ~15 hours | 150 ns total |
| Umbrella sampling (51 × 10 ns) | ~50 hours | 510 ns total |
| Analysis | ~1 hour | WHAM + MBAR + Jarzynski + bootstrap |
| **Total** | **~80 hours GPU** | — |

---

## Part 2 — Step-by-Step Implementation Instructions

> **Execution context:** All code runs in the project virtual environment with the `medium_project_2/` directory as the working directory. All imports reference the V2 pipeline as-is (§21). Parameter overrides are via runtime arguments only.

### Step 1: Environment and Directory Setup

```python
import os
import sys
import json
import numpy as np
from pathlib import Path
from dataclasses import replace

# Ensure project root is on PYTHONPATH
PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

# Create experiment output directory
EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-01_spink7_klk5_dg_bind"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Import pipeline modules
from src.config import (
    SystemConfig, MinimizationConfig, EquilibrationConfig,
    ProductionConfig, SMDConfig, UmbrellaConfig, WHAMConfig,
    MBARConfig, BOLTZMANN_KJ, KCAL_TO_KJ
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
from src.analyze.jarzynski import (
    jarzynski_free_energy, diagnose_dissipation,
    bar_free_energy, evaluate_convergence
)
from src.analyze.structural import compute_rmsd, compute_rmsf, compute_sasa
from src.analyze.trajectory import load_trajectory, align_trajectory
from src.analyze.convergence import block_average, effective_sample_size
from src.analyze.equilibration import detect_equilibration
from src.analyze.contacts import compute_interface_contacts, compute_hbonds
from src.physics.units import kj_to_kcal, nm_to_angstrom, kbt
from src.physics.finite_size import finite_size_correction
from src.visualization.plot_pmf import plot_pmf
from src.visualization.plot_timeseries import (
    plot_rmsd_timeseries, plot_energy_timeseries, plot_temperature_timeseries
)
```

**Expected output:** All imports succeed without error.

**Verification:**
```python
print("All imports successful.")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Figures directory: {FIGURES_DIR}")
```

**If this step fails:** Check that the virtual environment is activated and all dependencies in `requirements.txt` are installed. Verify `PROJECT_ROOT` path is correct.

---

### Step 2: Structure Acquisition and Preparation

```python
# 2.1 Fetch PDB structures
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

spink7_pdb = fetch_pdb("2LEO", data_dir)
klk5_pdb = fetch_pdb("2PSX", data_dir)

print(f"SPINK7 PDB: {spink7_pdb}")
print(f"KLK5 PDB: {klk5_pdb}")

# 2.2 Clean structures — keep chain A, model 1 from NMR ensemble
spink7_clean = clean_structure(
    pdb_path=spink7_pdb,
    chains_to_keep=["A"],
    remove_heteroatoms=True,
    remove_waters=True,
    model_index=1,
)
klk5_clean = clean_structure(
    pdb_path=klk5_pdb,
    chains_to_keep=["A"],
    remove_heteroatoms=True,
    remove_waters=True,
    model_index=1,
)

print(f"SPINK7 cleaned: {spink7_clean}")
print(f"KLK5 cleaned: {klk5_clean}")
```

**Expected output:** Two cleaned PDB files in the `structures/` subdirectory.

**Verification:**
```python
assert spink7_clean.exists(), "SPINK7 cleaned PDB not found"
assert klk5_clean.exists(), "KLK5 cleaned PDB not found"
```

**Error handling:** If fetch fails (network issue), retry once. If PDB IDs are invalid, halt — this is a critical failure.

---

### Step 3: Complex Assembly (ClusPro Docked Model)

> **Note:** The ClusPro-docked complex model from Azouz et al. (2020) must be prepared externally or reconstructed by docking SPINK7 onto KLK5 using the canonical Kazal-type serine protease inhibitor binding mode. If a pre-docked model is available in `data/pdb/`, use that directly.

```python
# Check for pre-existing docked complex
docked_complex_path = PROJECT_ROOT / "data" / "pdb" / "spink7_klk5_complex.pdb"

if not docked_complex_path.exists():
    # Construct by superposition of SPINK7 onto canonical binding site
    # The canonical binding mode places the RSL (reactive site loop) of SPINK7
    # into the KLK5 active site cleft, with the P1 residue (Lys/Arg) positioned
    # analogously to BPTI K15 in the trypsin active site.
    print("WARNING: Pre-docked complex not found. Manual docking required.")
    print("Place the ClusPro-docked SPINK7-KLK5 complex PDB at:")
    print(f"  {docked_complex_path}")
    # For pipeline testing, we proceed with whatever is available
else:
    print(f"Docked complex found: {docked_complex_path}")

# 3.1 Assign protonation states
complex_protonated = assign_protonation(
    pdb_path=docked_complex_path,
    ph=7.4,
    force_field="AMBER",
    use_propka=True,
)
print(f"Protonated complex: {complex_protonated}")
```

**Verification:**
```python
assert complex_protonated.exists(), "Protonated complex PDB not found"
```

**Intermediate test — P1 near catalytic Ser:** After loading the protonated structure, verify the P1 carbonyl oxygen is within 4 Å of the Ser195 Oγ (or equivalent residue number in KLK5).

```python
import mdtraj as md

traj = md.load(str(complex_protonated))
topology = traj.topology

# Identify chains
chains = list(topology.chains)
print(f"Number of chains: {len(chains)}")
for c in chains:
    print(f"  Chain {c.chain_id}: {c.n_residues} residues")
```

---

### Step 4: Topology Building and Solvation

```python
from openmm.app import PME

# 4.1 Use default SystemConfig (all values from config.py)
sys_config = SystemConfig()  # All defaults per §4

# 4.2 Build topology
topology, system, modeller = build_topology(
    pdb_path=complex_protonated,
    system_config=sys_config,
    nonbonded_method=PME,
    nonbonded_cutoff_nm=1.0,
)
print(f"Topology: {topology.getNumAtoms()} atoms")

# 4.3 Solvate
modeller, n_waters, n_pos_ions, n_neg_ions = solvate_system(
    modeller=modeller,
    system_config=sys_config,
)
print(f"Solvated system: {n_waters} waters, {n_pos_ions} Na+, {n_neg_ions} Cl-")
print(f"Total atoms: {modeller.topology.getNumAtoms()}")

# Save solvated topology for later use
from openmm.app import PDBFile
solvated_pdb = OUTPUT_DIR / "solvated_complex.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
```

**Verification:**
```python
assert n_waters > 1000, f"Too few waters: {n_waters}"
total_atoms = modeller.topology.getNumAtoms()
print(f"PASS: System has {total_atoms} atoms, {n_waters} waters")
```

**If this step fails:** Check force field XML paths, ensure PDB has standard residue names.

---

### Step 5: Energy Minimization

```python
from openmm.app import Simulation
from openmm import LangevinMiddleIntegrator
import openmm

# 5.1 Create simulation
platform = select_platform()
integrator = LangevinMiddleIntegrator(
    310 * openmm.unit.kelvin,
    1.0 / openmm.unit.picosecond,
    0.002 * openmm.unit.picoseconds,
)
sim = Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

# 5.2 Record initial energy
initial_state = sim.context.getState(getEnergy=True)
e_initial = initial_state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoules_per_mole)
print(f"Initial potential energy: {e_initial:.2f} kJ/mol")

# 5.3 Minimize
min_config = MinimizationConfig()  # defaults: max_iter=10000, tol=10.0
min_results = minimize_energy(sim, min_config)

# 5.4 Record final energy
final_state = sim.context.getState(getEnergy=True)
e_final = final_state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoules_per_mole)
print(f"Final potential energy: {e_final:.2f} kJ/mol")
```

**Verification — IV-1:**
```python
assert e_final < e_initial, f"IV-1 FAIL: E_min ({e_final}) >= E_initial ({e_initial})"
print(f"IV-1 PASS: Energy decreased by {e_initial - e_final:.2f} kJ/mol")
```

---

### Step 6: Equilibration (NVT + NPT)

```python
eq_config = EquilibrationConfig()  # All defaults per §4
eq_output = OUTPUT_DIR / "equilibration"
eq_output.mkdir(exist_ok=True)

# 6.1 NVT equilibration (500 ps)
nvt_results = run_nvt(sim, eq_config, eq_output)
print(f"NVT complete: {nvt_results}")

# Verify IV-2: Temperature
t_avg = nvt_results.get("average_temperature_k", 0)
assert abs(t_avg - 310.0) < 5.0, f"IV-2 FAIL: |T_avg - 310| = {abs(t_avg - 310.0):.2f} K"
print(f"IV-2 PASS: T_avg = {t_avg:.2f} K")

# 6.2 NPT equilibration (1000 ps)
npt_results = run_npt(sim, eq_config, eq_output)
print(f"NPT complete: {npt_results}")

# Verify IV-3: Density
rho = npt_results.get("average_density_g_per_cm3", 0)
assert 0.95 <= rho <= 1.05, f"IV-3 FAIL: density = {rho:.4f} g/cm³"
print(f"IV-3 PASS: density = {rho:.4f} g/cm³")

# Save equilibrated state
eq_state_path = eq_output / "equilibrated_state.xml"
sim.saveState(str(eq_state_path))

# Save system XML
system_xml_path = eq_output / "system.xml"
with open(system_xml_path, "w") as f:
    f.write(openmm.XmlSerializer.serialize(system))

print(f"Equilibrated state saved: {eq_state_path}")
print(f"System XML saved: {system_xml_path}")
```

**Verification — IV-4 (structural stability):**
```python
# Load equilibrated trajectory and check backbone RMSD
eq_traj = load_trajectory(
    trajectory_path=eq_output / "npt_trajectory.dcd",
    topology_path=solvated_pdb,
)
ref_traj = md.load(str(solvated_pdb))
rmsd_eq = compute_rmsd(eq_traj, ref_traj, atom_selection="backbone")
max_rmsd = np.max(rmsd_eq)
print(f"Max backbone RMSD during equilibration: {nm_to_angstrom(max_rmsd):.2f} Å")
assert max_rmsd < 0.5, f"IV-4 FAIL: max RMSD = {nm_to_angstrom(max_rmsd):.2f} Å > 5 Å"
print("IV-4 PASS: Backbone RMSD < 5 Å")
```

---

### Step 7: Production MD (100 ns)

```python
prod_config = ProductionConfig()  # defaults: 100 ns
prod_output = OUTPUT_DIR / "production"
prod_output.mkdir(exist_ok=True)

prod_results = run_production(sim, prod_config, prod_output)
print(f"Production MD complete: {prod_results}")
```

**Verification — IV-5 (energy drift):**
```python
energy_drift = prod_results.get("energy_drift_kj_mol_ns_atom", 0)
assert abs(energy_drift) < 0.1, f"IV-5 FAIL: energy drift = {energy_drift:.4f} kJ/mol/ns/atom"
print(f"IV-5 PASS: energy drift = {energy_drift:.6f} kJ/mol/ns/atom")
```

**10 ns checkpoint verification:**
```python
prod_traj = load_trajectory(
    trajectory_path=prod_output / "production_trajectory.dcd",
    topology_path=solvated_pdb,
    stride=100,  # every 1 ns for quick check
)
rmsd_prod = compute_rmsd(prod_traj, ref_traj, atom_selection="backbone")
print(f"Production RMSD range: {nm_to_angstrom(rmsd_prod.min()):.2f} - {nm_to_angstrom(rmsd_prod.max()):.2f} Å")

# Verify RMSD at 10 ns
rmsd_10ns = rmsd_prod[min(10, len(rmsd_prod)-1)]
assert rmsd_10ns < 0.5, f"Complex unstable at 10 ns: RMSD = {nm_to_angstrom(rmsd_10ns):.2f} Å"
print(f"10 ns checkpoint PASS: RMSD = {nm_to_angstrom(rmsd_10ns):.2f} Å")
```

---

### Step 8: Identify Pull Groups

```python
# Identify SPINK7 and KLK5 atom indices from topology
prod_traj_full = load_trajectory(
    trajectory_path=prod_output / "production_trajectory.dcd",
    topology_path=solvated_pdb,
)
topology_md = prod_traj_full.topology

# Select pull groups by chain
chains_list = list(topology_md.chains)
chain_a_indices = topology_md.select("chainid 0")  # SPINK7
chain_b_indices = topology_md.select("chainid 1")  # KLK5

pull_group_1 = chain_a_indices.tolist()  # SPINK7
pull_group_2 = chain_b_indices.tolist()  # KLK5

print(f"Pull group 1 (SPINK7): {len(pull_group_1)} atoms")
print(f"Pull group 2 (KLK5): {len(pull_group_2)} atoms")

# Compute initial COM distance
from src.physics.collective_variables import com_distance
positions = prod_traj_full.xyz[-1]  # last frame
masses = np.array([a.element.mass for a in topology_md.atoms])

xi_initial = com_distance(
    positions, masses,
    np.array(pull_group_1), np.array(pull_group_2),
)
print(f"Initial COM-COM distance: {xi_initial:.4f} nm ({nm_to_angstrom(xi_initial):.2f} Å)")
```

---

### Step 9: SMD Campaign (50 Replicates)

```python
smd_config = SMDConfig()  # defaults: k=1000, v=0.001, d=3.0, n=50
smd_output = OUTPUT_DIR / "smd"
smd_output.mkdir(exist_ok=True)

# Compute pull direction (away from KLK5 COM)
com_spink7 = np.mean(positions[pull_group_1], axis=0)
com_klk5 = np.mean(positions[pull_group_2], axis=0)
pull_direction = com_spink7 - com_klk5
pull_direction = pull_direction / np.linalg.norm(pull_direction)

smd_results = run_smd_campaign(
    equilibrated_state_path=eq_state_path,
    system_xml_path=system_xml_path,
    config=smd_config,
    pull_group_1=pull_group_1,
    pull_group_2=pull_group_2,
    output_dir=smd_output,
    topology_pdb_path=solvated_pdb,
)

print(f"SMD campaign complete: {len(smd_results)} replicates")
```

**Verification — work values:**
```python
work_values = np.array([r["total_work_kj_mol"] for r in smd_results])
print(f"Work values: mean = {work_values.mean():.2f} kJ/mol, "
      f"std = {work_values.std():.2f} kJ/mol")

# Check all work values are finite and positive
assert np.all(np.isfinite(work_values)), "FAIL: Non-finite work values detected"
assert np.all(work_values > 0), "FAIL: Negative work values detected"
print("PASS: All work values finite and positive")
```

---

### Step 10: SMD Analysis (Jarzynski + BAR)

```python
temperature_k = 310.0

# 10.1 Jarzynski free energy estimate
jarz_results = jarzynski_free_energy(work_values, temperature_k)
dg_jarz_kj = jarz_results["delta_g_kj_mol"]
dg_jarz_kcal = kj_to_kcal(dg_jarz_kj)
print(f"Jarzynski ΔG (exponential average): {dg_jarz_kcal:.2f} kcal/mol")

# 10.2 Cumulant expansion
dg_cumulant_kj = jarz_results.get("delta_g_cumulant_kj_mol", None)
if dg_cumulant_kj is not None:
    dg_cumulant_kcal = kj_to_kcal(dg_cumulant_kj)
    print(f"Jarzynski ΔG (cumulant expansion): {dg_cumulant_kcal:.2f} kcal/mol")

# 10.3 Dissipation diagnostics
dissipation = diagnose_dissipation(work_values, temperature_k)
print(f"Dissipation: {dissipation}")

# 10.4 Convergence analysis
conv_results = evaluate_convergence(work_values, temperature_k, n_subsets=10)
print("Convergence analysis complete")

# 10.5 Check IV-10 (unimodal work distribution)
from scipy.stats import kurtosis
work_kurt = kurtosis(work_values)
print(f"Work distribution kurtosis: {work_kurt:.2f}")
# High kurtosis (>6) suggests multimodality
if work_kurt > 6.0:
    print("WARNING: IV-10 — work distribution may be multimodal")
else:
    print("IV-10 PASS: Work distribution appears unimodal")
```

---

### Step 11: Umbrella Sampling Campaign (51 Windows)

```python
us_config = UmbrellaConfig()  # defaults: 51 windows, 10 ns each
us_output = OUTPUT_DIR / "umbrella"
us_output.mkdir(exist_ok=True)

window_centers = generate_window_centers(us_config)
print(f"Umbrella windows: {len(window_centers)} windows from "
      f"{window_centers[0]:.3f} to {window_centers[-1]:.3f} nm")

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

print(f"Umbrella sampling complete: {len(us_results)} windows")
```

**Verification — IV-8 (histogram overlap):**
```python
# Collect ξ timeseries from each window
xi_timeseries_list = []
for i, result in enumerate(us_results):
    xi_ts = np.load(us_output / f"window_{i:03d}" / "xi_timeseries.npy")
    xi_timeseries_list.append(xi_ts)

# Diagnose histogram coverage
coverage = diagnose_histogram_coverage(xi_timeseries_list, window_centers)
print(f"Coverage analysis: {coverage}")

# Compute overlap matrix
overlap = compute_overlap_matrix(xi_timeseries_list)
min_adjacent_overlap = min(overlap[i, i+1] for i in range(len(overlap)-1))
print(f"Minimum adjacent overlap: {min_adjacent_overlap:.4f}")
assert min_adjacent_overlap >= 0.10, f"IV-8 FAIL: min overlap = {min_adjacent_overlap:.4f} < 0.10"
print("IV-8 PASS: All adjacent windows have ≥10% overlap")
```

---

### Step 12: WHAM Analysis

```python
wham_config = WHAMConfig()  # defaults: tol=1e-6, max_iter=100000, n_boot=200
spring_constants = np.full(len(window_centers), us_config.spring_constant_kj_mol_nm2)

# 12.1 Solve WHAM
wham_results = solve_wham(
    xi_timeseries_list=xi_timeseries_list,
    window_centers=window_centers,
    spring_constants=spring_constants,
    temperature_k=temperature_k,
    config=wham_config,
)

pmf_wham = wham_results["pmf_kj_mol"]
xi_bins = wham_results["xi_bins_nm"]
converged = wham_results["converged"]
print(f"WHAM converged: {converged}")

# Verify IV-9
assert converged, "IV-9 FAIL: WHAM did not converge"
print("IV-9 PASS: WHAM converged")

# 12.2 Bootstrap uncertainty
wham_bootstrap = bootstrap_pmf_uncertainty(
    xi_timeseries_list=xi_timeseries_list,
    window_centers=window_centers,
    spring_constants=spring_constants,
    temperature_k=temperature_k,
    config=wham_config,
)
pmf_std = wham_bootstrap["pmf_std_kj_mol"]
```

---

### Step 13: MBAR Analysis

```python
mbar_config = MBARConfig()  # defaults: robust, tol=1e-7

mbar_results = solve_mbar(
    xi_timeseries_list=xi_timeseries_list,
    window_centers=window_centers,
    spring_constants=spring_constants,
    temperature_k=temperature_k,
    config=mbar_config,
)

pmf_mbar = mbar_results["pmf_kj_mol"]
print("MBAR analysis complete")

# Bootstrap
mbar_bootstrap = bootstrap_mbar_uncertainty(
    xi_timeseries_list=xi_timeseries_list,
    window_centers=window_centers,
    spring_constants=spring_constants,
    temperature_k=temperature_k,
    config=mbar_config,
)
pmf_mbar_std = mbar_bootstrap["pmf_std_kj_mol"]
```

---

### Step 14: Binding Free Energy Extraction

```python
# Extract ΔG_bind from PMF using standard concentration correction (§5.5)
kT = kbt(temperature_k)  # kJ/mol
beta = 1.0 / kT

# Standard concentration: C° = 1 M = 1/1660 Å⁻³ = 1/(1660 * 0.001) nm⁻³
C_standard_nm3 = 1.0 / (1660.0 * 0.001)  # nm⁻³

# Define bound and bulk regions
bound_mask = xi_bins < 2.0  # bound state: ξ < 2.0 nm
bulk_mask = xi_bins > 3.5   # bulk/unbound: ξ > 3.5 nm

# Compute integrals
# ΔG_bind = -kT ln[ (C°/4π) ∫_site exp(-β*G) ξ² dξ ]
#           + kT ln[ (C°/4π) ∫_bulk exp(-β*G) ξ² dξ ]

dxi = xi_bins[1] - xi_bins[0]

# WHAM-based binding free energy
pmf_wham_shifted = pmf_wham - pmf_wham[bulk_mask].min()  # shift bulk to zero
integrand_site = np.exp(-beta * pmf_wham_shifted[bound_mask]) * xi_bins[bound_mask]**2
integrand_bulk = np.exp(-beta * pmf_wham_shifted[bulk_mask]) * xi_bins[bulk_mask]**2

I_site = np.trapz(integrand_site, xi_bins[bound_mask])
I_bulk = np.trapz(integrand_bulk, xi_bins[bulk_mask])

dg_bind_wham_kj = -kT * np.log(C_standard_nm3 / (4 * np.pi) * I_site) + \
                   kT * np.log(C_standard_nm3 / (4 * np.pi) * I_bulk)
dg_bind_wham_kcal = kj_to_kcal(dg_bind_wham_kj)

print(f"ΔG_bind (WHAM): {dg_bind_wham_kcal:.2f} kcal/mol")

# MBAR-based binding free energy
pmf_mbar_shifted = pmf_mbar - pmf_mbar[bulk_mask].min()
integrand_site_m = np.exp(-beta * pmf_mbar_shifted[bound_mask]) * xi_bins[bound_mask]**2
integrand_bulk_m = np.exp(-beta * pmf_mbar_shifted[bulk_mask]) * xi_bins[bulk_mask]**2

I_site_m = np.trapz(integrand_site_m, xi_bins[bound_mask])
I_bulk_m = np.trapz(integrand_bulk_m, xi_bins[bulk_mask])

dg_bind_mbar_kj = -kT * np.log(C_standard_nm3 / (4 * np.pi) * I_site_m) + \
                   kT * np.log(C_standard_nm3 / (4 * np.pi) * I_bulk_m)
dg_bind_mbar_kcal = kj_to_kcal(dg_bind_mbar_kj)

print(f"ΔG_bind (MBAR): {dg_bind_mbar_kcal:.2f} kcal/mol")

# Internal consistency check: WHAM vs MBAR
wham_mbar_diff = abs(dg_bind_wham_kcal - dg_bind_mbar_kcal)
print(f"|WHAM - MBAR| = {wham_mbar_diff:.2f} kcal/mol")
if wham_mbar_diff < 1.0:
    print("PASS: WHAM and MBAR agree within 1 kcal/mol")
else:
    print("WARNING: WHAM and MBAR differ by more than 1 kcal/mol")
```

---

### Step 15: Statistical Classification

```python
# Experimental benchmark
dg_exp = -9.4  # kcal/mol (from Ki = 132 nM)
sigma_exp = 0.6  # kcal/mol

# Computational uncertainty (from bootstrap)
# Use WHAM bootstrap std of ΔG_bind
sigma_comp = kj_to_kcal(np.mean(pmf_std))  # approximate

# Method systematic error
sigma_method = 2.0  # kcal/mol (US/WHAM for protein-protein, §25.4)

# Combined uncertainty
sigma_combined = np.sqrt(sigma_exp**2 + sigma_comp**2 + sigma_method**2)

# Classification (§25.1)
z = 1.96  # z_{0.975}
discrepancy = abs(dg_bind_wham_kcal - dg_exp)

if discrepancy <= z * sigma_combined:
    classification = "PASS"
elif discrepancy <= 2 * z * sigma_combined:
    classification = "MARGINAL"
else:
    classification = "FAIL"

# Check for INCONCLUSIVE
if sigma_comp > sigma_exp:
    classification = "INCONCLUSIVE"

print(f"\n=== CLASSIFICATION ===")
print(f"Pipeline prediction (WHAM): {dg_bind_wham_kcal:.2f} kcal/mol")
print(f"Experimental benchmark: {dg_exp:.2f} ± {sigma_exp:.2f} kcal/mol")
print(f"σ_combined: {sigma_combined:.2f} kcal/mol")
print(f"Discrepancy: {discrepancy:.2f} kcal/mol")
print(f"Classification: {classification}")

# Save results
results = {
    "experiment_id": "EXP-01",
    "feature_id": "F-01",
    "dg_bind_wham_kcal": float(dg_bind_wham_kcal),
    "dg_bind_mbar_kcal": float(dg_bind_mbar_kcal),
    "dg_jarz_exp_avg_kcal": float(dg_jarz_kcal),
    "dg_exp_kcal": float(dg_exp),
    "sigma_exp": float(sigma_exp),
    "sigma_comp": float(sigma_comp),
    "sigma_method": float(sigma_method),
    "sigma_combined": float(sigma_combined),
    "discrepancy": float(discrepancy),
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {EXP_DIR / 'results.json'}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Experimental Design Schematic

**Description:** Thermodynamic cycle showing the SMD pulling pathway and umbrella sampling windows along the COM-COM reaction coordinate ξ.

```python
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Draw schematic of pulling pathway
xi_range = np.linspace(1.5, 4.5, 200)
# Idealized PMF shape
pmf_schematic = 15 * np.exp(-2 * (xi_range - 1.8)**2) - 10 * np.exp(-0.5 * (xi_range - 1.8)**2)
pmf_schematic -= pmf_schematic[-1]

ax.plot(xi_range, pmf_schematic, 'k-', linewidth=2, label='PMF (schematic)')
ax.fill_between(xi_range, pmf_schematic, alpha=0.1, color='blue')

# Mark umbrella windows
for center in np.arange(1.5, 4.05, 0.05):
    ax.axvline(center, color='red', alpha=0.15, linewidth=0.5)

# Annotations
ax.annotate('Bound state', xy=(1.8, pmf_schematic.min()), fontsize=11,
            ha='center', va='top', color='green', fontweight='bold')
ax.annotate('Unbound\n(bulk)', xy=(4.0, 0), fontsize=11,
            ha='center', va='bottom', color='gray')
ax.set_xlabel(r'COM-COM distance $\xi$ (nm)', fontsize=13)
ax.set_ylabel(r'PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-01: SPINK7-KLK5 Binding Free Energy — Experimental Design', fontsize=14)
ax.legend(fontsize=11)
ax.set_xlim(1.3, 4.7)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-01_experimental_design_schematic.png", dpi=300)
plt.close(fig)
print("Figure 1 saved.")
```

### Figure 2: PMF Profile (WHAM + MBAR)

**Description:** Reconstructed PMF from WHAM and MBAR with bootstrap uncertainty bands, showing the bound-state minimum and unbound plateau.

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

pmf_wham_kcal = pmf_wham / KCAL_TO_KJ
pmf_std_kcal = pmf_std / KCAL_TO_KJ
pmf_mbar_kcal = pmf_mbar / KCAL_TO_KJ
pmf_mbar_std_kcal = pmf_mbar_std / KCAL_TO_KJ

ax.plot(xi_bins, pmf_wham_kcal, 'b-', linewidth=2, label='WHAM')
ax.fill_between(xi_bins, pmf_wham_kcal - 1.96*pmf_std_kcal,
                pmf_wham_kcal + 1.96*pmf_std_kcal, alpha=0.2, color='blue')

ax.plot(xi_bins, pmf_mbar_kcal, 'r--', linewidth=2, label='MBAR')
ax.fill_between(xi_bins, pmf_mbar_kcal - 1.96*pmf_mbar_std_kcal,
                pmf_mbar_kcal + 1.96*pmf_mbar_std_kcal, alpha=0.2, color='red')

ax.set_xlabel(r'COM-COM distance $\xi$ (nm)', fontsize=13)
ax.set_ylabel(r'PMF (kcal/mol)', fontsize=13)
ax.set_title('EXP-01: PMF Profile — WHAM vs. MBAR', fontsize=14)
ax.legend(fontsize=12)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-01_pmf_profile.png", dpi=300)
plt.close(fig)
print("Figure 2 saved.")
```

### Figure 3: SMD Work Distribution

**Description:** Histogram of work values from 50 SMD replicates, with Jarzynski ΔG estimate marked.

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

work_kcal = work_values / KCAL_TO_KJ

ax.hist(work_kcal, bins=20, edgecolor='black', alpha=0.7, color='steelblue',
        label=f'N = {len(work_kcal)} replicates')
ax.axvline(dg_jarz_kcal, color='red', linestyle='--', linewidth=2,
           label=f'Jarzynski ΔG = {dg_jarz_kcal:.1f} kcal/mol')
ax.axvline(np.mean(work_kcal), color='green', linestyle='-', linewidth=2,
           label=f'⟨W⟩ = {np.mean(work_kcal):.1f} kcal/mol')

ax.set_xlabel('Work (kcal/mol)', fontsize=13)
ax.set_ylabel('Count', fontsize=13)
ax.set_title('EXP-01: SMD Work Distribution', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-01_smd_work_distribution.png", dpi=300)
plt.close(fig)
print("Figure 3 saved.")
```

### Figure 4: Predicted vs. Experimental Comparison

**Description:** Bar chart comparing pipeline ΔG predictions (WHAM, MBAR, Jarzynski) with the experimental value, including error bars and CI bounds.

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

methods = ['Experimental', 'US/WHAM', 'US/MBAR', 'SMD/Jarzynski']
values = [dg_exp, dg_bind_wham_kcal, dg_bind_mbar_kcal, dg_jarz_kcal]
errors = [sigma_exp, kj_to_kcal(np.mean(pmf_std)),
          kj_to_kcal(np.mean(pmf_mbar_std)), work_kcal.std() / np.sqrt(len(work_kcal))]
colors = ['gold', 'steelblue', 'salmon', 'mediumseagreen']

bars = ax.bar(methods, values, yerr=[1.96*e for e in errors], capsize=8,
              color=colors, edgecolor='black', alpha=0.85)

# Draw CI bounds
ax.axhspan(-14.0, -4.8, alpha=0.1, color='green', label='95% CI (PASS)')
ax.axhline(-9.4, color='black', linestyle=':', linewidth=1.5, label='Experimental ΔG')

ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)', fontsize=13)
ax.set_title('EXP-01: Predicted vs. Experimental ΔG_bind', fontsize=14)
ax.legend(fontsize=10)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-01_predicted_vs_experimental.png", dpi=300)
plt.close(fig)
print("Figure 4 saved.")
```

### Figure 5: Umbrella Window Histogram Overlap

**Description:** Stacked histogram showing the sampling distribution for each umbrella window, demonstrating sufficient overlap for WHAM convergence.

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 5))

for i, xi_ts in enumerate(xi_timeseries_list):
    ax.hist(xi_ts, bins=50, alpha=0.3, density=True, label=None)

ax.set_xlabel(r'COM-COM distance $\xi$ (nm)', fontsize=13)
ax.set_ylabel('Probability density', fontsize=13)
ax.set_title('EXP-01: Umbrella Window Histogram Overlap', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-01_umbrella_overlap.png", dpi=300)
plt.close(fig)
print("Figure 5 saved.")
```

### Figure 6: Jarzynski Convergence

**Description:** Cumulative Jarzynski ΔG estimate as a function of the number of SMD trajectories, showing convergence behavior.

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

n_traj_array = conv_results["n_trajectories"]
dg_cumulative = conv_results["delta_g_cumulative_kj_mol"]
dg_cumulative_kcal = dg_cumulative / KCAL_TO_KJ

ax.plot(n_traj_array, dg_cumulative_kcal, 'b-o', markersize=4, linewidth=1.5)
ax.axhline(dg_exp, color='red', linestyle='--', label=f'Experimental: {dg_exp:.1f} kcal/mol')
ax.fill_between(n_traj_array, -14.0, -4.8, alpha=0.1, color='green', label='95% CI')

ax.set_xlabel('Number of SMD trajectories', fontsize=13)
ax.set_ylabel(r'$\Delta G$ (kcal/mol)', fontsize=13)
ax.set_title('EXP-01: Jarzynski Convergence', fontsize=14)
ax.legend(fontsize=11)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-01_jarzynski_convergence.png", dpi=300)
plt.close(fig)
print("Figure 6 saved.")
```

### Figure 7: Production MD RMSD Time Series

**Description:** Backbone RMSD over the 100 ns production trajectory, verifying structural stability of the complex.

```python
time_ns = np.arange(len(rmsd_prod)) * (prod_config.save_interval_ps / 1000.0)

fig = plot_rmsd_timeseries(
    time_ns=time_ns,
    rmsd_nm=rmsd_prod,
    output_path=FIGURES_DIR / "EXP-01_production_rmsd.png",
)
plt.close(fig)
print("Figure 7 saved.")
```

---

## Part 4 — Results Documentation Template

The results report (`results_report.md`) must follow this structure per §28.1:

```markdown
# EXP-01: SPINK7-KLK5 Binding Free Energy — Results Report

**Experiment ID:** EXP-01
**Feature ID:** F-01
**Date:** [execution date]
**Classification:** [PASS / MARGINAL / FAIL / INCONCLUSIVE]

## 1. Abstract
[One-paragraph summary: system, methods, key numerical result, classification]

## 2. Introduction/Background
[Scientific context for SPINK7-KLK5 ΔG_bind; cite Azouz et al. 2020]

## 3. Hypothesis
[Restate H₁, H₂, H₃ from experimental design]

## 4. Methods
### 4.1 System Preparation
[PDB IDs, cleaning steps, protonation, solvation parameters — all explicit]
### 4.2 Production MD
[Duration, temperature, all parameters as executed]
### 4.3 Steered Molecular Dynamics
[k, v, d, N_rep — as executed]
### 4.4 Umbrella Sampling
[Windows, duration, spring constant — as executed]
### 4.5 Analysis
[WHAM, MBAR, Jarzynski parameters — as executed]
### 4.6 Deviations from Experimental Design
[Any deviations with justification, or "None"]

## 5. Controls
[Results of positive control (EXP-04) and negative controls/sanity checks]

## 6. Results
### 6.1 Pipeline Predictions
| Method | ΔG_bind (kcal/mol) | 95% CI |
|--------|-------------------|---------|
| US/WHAM | X.XX ± Y.YY | [A, B] |
| US/MBAR | X.XX ± Y.YY | [A, B] |
| SMD/Jarzynski (exp avg) | X.XX ± Y.YY | [A, B] |
| SMD/Jarzynski (cumulant) | X.XX ± Y.YY | [A, B] |

### 6.2 Experimental Benchmark
ΔG_exp = −9.4 ± 0.6 kcal/mol (Ki = 132 nM, Azouz et al. 2020)

### 6.3 Discrepancy and Classification
σ_combined = X.XX kcal/mol
Discrepancy = X.XX kcal/mol
**Classification: [PASS / MARGINAL / FAIL / INCONCLUSIVE]**

### 6.4 Internal Consistency
|WHAM - MBAR| = X.XX kcal/mol
|Jarzynski exp avg - cumulant| = X.XX kcal/mol

## 7. Discussion
[Interpretation of results; if FAIL, root-cause analysis per §25.7]

## 8. Conclusions
[One-paragraph summary]

## 9. Figures
[≥ 7 figures embedded with captions]

## 10. References
[Full citations]

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
