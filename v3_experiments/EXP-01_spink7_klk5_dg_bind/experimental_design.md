# EXP-01: SPINK7-KLK5 Binding Free Energy

**Experiment ID:** EXP-01  
**Feature ID:** F-01 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** PRIMARY PIPELINE TARGET  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment determines the binding free energy (ΔG_bind) of the SPINK7-KLK5 protease-antiprotease complex using the V2 pipeline's two complementary enhanced sampling strategies: Steered Molecular Dynamics (SMD) with the Jarzynski equality and Umbrella Sampling (US) with WHAM reconstruction of the Potential of Mean Force (PMF). The experimentally measured Ki = 132 nM (ΔG ≈ −9.4 kcal/mol) from Azouz et al. (2020) serves as the primary benchmark. This is the single most important validation target for the entire V3 benchmarking project, as it directly tests whether the pipeline can predict the binding thermodynamics of its primary target system. Cross-validation using the Bennett Acceptance Ratio (BAR) estimator and MBAR will provide internal consistency checks across free energy estimation methods.

---

## 2. Hypothesis

**H₁:** The V2 pipeline's Umbrella Sampling + WHAM estimate of ΔG_bind for the SPINK7-KLK5 complex will fall within the 95% confidence interval [−14.0, −4.8] kcal/mol of the experimentally measured value (ΔG = −9.4 kcal/mol, derived from Ki = 132 nM via ΔG = RT ln(Ki)).

**H₂:** The V2 pipeline's SMD + Jarzynski estimate of ΔG_bind will fall within the 95% CI [−14.0, −4.8] kcal/mol, with the second-order cumulant expansion providing a more precise estimate than the exponential average.

**H₃:** All three free energy estimators (Jarzynski exponential average, cumulant expansion, and BAR) applied to the SMD work distributions will agree within 2 kcal/mol, and the US/WHAM and MBAR estimates from umbrella sampling data will agree within 1 kcal/mol.

These hypotheses are falsifiable: if any estimate falls outside the stated CI, the hypothesis is rejected and a pipeline limitation is identified.

---

## 3. Background and Rationale

### 3.1 Scientific Context

SPINK7 (Serine Peptidase Inhibitor, Kazal Type 7) is a ~6 kDa secreted Kazal-type inhibitor that stoichiometrically inhibits KLK5, a trypsin-like serine protease. In Eosinophilic Esophagitis (EoE), IL-13-mediated transcriptional silencing of SPINK7 unleashes KLK5 proteolytic activity, degrading Desmoglein-1 and compromising the epithelial barrier. The SPINK7-KLK5 Ki = 132 nM was measured using the Morrison tight-binding inhibitor equation from dose-response curves (Azouz et al. 2020, Fig. 1B), with three independent experiments performed in duplicate.

### 3.2 Pipeline Capability

The V2 pipeline provides two complementary enhanced sampling pathways for computing ΔG_bind:

1. **SMD + Jarzynski:** Non-equilibrium pulling of SPINK7 away from KLK5 along the COM-COM reaction coordinate ξ, with free energy extracted via the Jarzynski equality (§5.3). The pipeline implements the exponential average, second-order cumulant expansion, and BAR estimators (`src/analyze/jarzynski.py`).

2. **Umbrella Sampling + WHAM:** Equilibrium sampling of ξ in $M$ biased windows with harmonic restraints, reconstructing the PMF via WHAM (`src/analyze/wham.py`) and MBAR (`src/analyze/mbar.py`). The binding free energy is extracted from the PMF using the standard concentration correction (§5.5).

### 3.3 What This Reveals

This experiment tests the pipeline's core capability — predicting protein-protein binding free energies from first principles. Since no co-crystal structure of SPINK7-KLK5 exists, the starting structure relies on ClusPro docking (Azouz et al. 2020, Fig. 1C), adding model uncertainty. A PASS result validates both the computational methodology and the homology-based complex model; a FAIL result identifies specific limitations (force field, sampling, or model quality) for Part 2 remediation.

---

## 4. Experimental Protocol

### 4.1 System Preparation

#### 4.1.1 Structure Acquisition

| Parameter | Value |
|-----------|-------|
| SPINK7 structure | PDB 2LEO (NMR ensemble, chain A, model 1) |
| KLK5 structure | PDB 2PSX (X-ray, chain A) |
| Complex model | ClusPro docking of 2LEO onto 2PSX (canonical binding mode from Azouz et al. 2020) |
| Missing residues | Fill using PDBFixer (`src/prep/pdb_clean.py`) |
| Disulfide bonds | Verify and enforce all three SPINK7 disulfide bonds (IV-6) |

#### 4.1.2 Protonation and Solvation

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

#### 4.1.3 Energy Minimization

| Parameter | Value | Source |
|-----------|-------|--------|
| Max iterations | 10,000 | `config.py: MinimizationConfig.max_iterations` |
| Tolerance | 10.0 kJ/mol/nm | `config.py: MinimizationConfig.tolerance_kj_mol_nm` |
| Invariant check | IV-1: $E_{\text{min}} < E_{\text{initial}}$ | §6 |

#### 4.1.4 Equilibration

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

### 4.2 Production MD (Pre-Pulling Equilibrium)

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

### 4.3 Steered Molecular Dynamics (SMD)

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

### 4.4 Umbrella Sampling

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

### 4.5 WHAM Analysis

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

### 4.6 MBAR Analysis

| Parameter | Value | Source |
|-----------|-------|--------|
| Solver protocol | "robust" | `config.py: MBARConfig.solver_protocol` |
| Relative tolerance | 10⁻⁷ | `config.py: MBARConfig.relative_tolerance` |
| Max iterations | 10,000 | `config.py: MBARConfig.maximum_iterations` |
| Bootstrap resamples | 200 | `config.py: MBARConfig.n_bootstrap` |
| PMF bins | 200 | `config.py: MBARConfig.n_pmf_bins` |

### 4.7 Binding Free Energy Extraction (§5.5)

$$\Delta G_{\text{bind}}^{\circ} = -k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{site}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right] + k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{bulk}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right]$$

where $C^{\circ} = 1/1660$ Å⁻³ (standard concentration of 1 M).

### 4.8 Statistical Comparison

1. Compute pipeline prediction: ΔG_bind ± σ_comp from bootstrap CI.
2. Construct the combined 95% CI using the §25 framework:
   - σ_exp = 0.6 kcal/mol
   - σ_comp = from pipeline bootstrap
   - σ_method = 2.0 kcal/mol (US/WHAM for protein-protein binding, §25.4)
   - σ_combined = √(σ_exp² + σ_comp² + σ_method²)
3. Classify result: PASS / MARGINAL / FAIL / INCONCLUSIVE per §25.1.

---

## 5. Control Conditions

### 5.1 Positive Control

**BPTI-trypsin (EXP-04):** The gold-standard protease-inhibitor system (Kd = 6 × 10⁻¹⁴ M, ΔG ≈ −18 kcal/mol) with a high-resolution co-crystal structure (PDB 2PTC). If the pipeline fails to reproduce the BPTI-trypsin ΔG_bind within its wider CI [−24.7, −11.3], the SPINK7-KLK5 result should be interpreted with caution regardless of its own classification.

**SH3-p41 (EXP-29):** Computational methods validation target with three independent methods converging within ~0.2 kcal/mol (experimental ΔG = −7.99 kcal/mol). This verifies the PMF methodology itself.

### 5.2 Negative Control / Sanity Checks

1. **Energy conservation:** Verify IV-5 (energy drift < 0.1 kJ/mol/ns/atom) during production MD.
2. **Structural stability:** SPINK7-KLK5 complex backbone RMSD < 5 Å during production (IV-4).
3. **Disulfide integrity:** All three SPINK7 disulfide bonds remain intact (IV-6, $d_{S-S} < 2.5$ Å).
4. **Pulling direction:** Confirm ξ increases monotonically during SMD (no re-association).
5. **PMF shape:** The PMF should show a clear minimum at the bound state (ξ ≈ 1.5–2.0 nm) and plateau at large separation (ξ > 3.5 nm).

---

## 6. Expected Outcomes

### 6.1 Primary Prediction

| Method | Expected Range | Source |
|--------|---------------|--------|
| US/WHAM ΔG_bind | −14.0 to −4.8 kcal/mol (95% CI) | benchmarks.md F-01 |
| SMD/Jarzynski ΔG_bind | −14.0 to −4.8 kcal/mol (95% CI) | benchmarks.md F-01 |

### 6.2 Classification Criteria

| Classification | US/WHAM ΔG_bind Range |
|---------------|----------------------|
| **PASS** | [−14.0, −4.8] kcal/mol |
| **MARGINAL** | [−18.7, −0.1] kcal/mol |
| **FAIL** | Outside [−18.7, −0.1] kcal/mol |
| **INCONCLUSIVE** | σ_comp > σ_exp (0.6 kcal/mol) |

### 6.3 Internal Consistency Checks

- WHAM and MBAR estimates should agree within 1 kcal/mol.
- Jarzynski exponential average and cumulant expansion should agree within 2 kcal/mol.
- Forward and reverse cumulative ΔG estimates should converge to plateau values.

---

## 7. Potential Failure Modes

| Failure Mode | Expected Manifestation | Pipeline Limitation Implied | Severity |
|-------------|----------------------|---------------------------|----------|
| **Homology model inaccuracy** | ΔG too positive (weak binding); complex dissociates during equilibration | ClusPro docking model geometry does not capture correct binding pose | Critical — no co-crystal structure available |
| **Insufficient SMD replicates** | Jarzynski ΔG not converged; exponential average dominated by rare low-work trajectories | 50 replicates insufficient for protein-protein unbinding | Medium — can increase to 100+ |
| **Inadequate US window overlap** | WHAM convergence failure (IV-9 violation) or rugged PMF with unphysical barriers | Window spacing too coarse or per-window sampling too short | Medium — can refine spacing |
| **Force field limitation** | Systematic over/underestimation of ΔG by >5 kcal/mol despite convergence | AMBER ff14SB fixed-charge model inadequate for charge-transfer interactions at protease-inhibitor interface | High — requires Part 2 force field investigation |
| **Pulling speed too fast** | SMD ΔG systematically too positive (dissipative work) | v = 0.001 nm/ps too fast for near-equilibrium pulling | Medium — known SMD limitation |
| **Finite-size effects** | Systematic bias from periodic boundary conditions | Box padding insufficient for 3 nm pull distance | Medium — finite-size corrections available |
| **Disulfide bond instability** | SPINK7 structure deforms during simulation (IV-6 violation) | Force field or protonation state assignment incorrect for Cys residues | Critical |

---

## 8. Intermediate Verification Tests

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

---

## 9. Resource Estimates

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

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
