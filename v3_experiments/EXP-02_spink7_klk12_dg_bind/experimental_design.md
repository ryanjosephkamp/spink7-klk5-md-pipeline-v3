# EXP-02: SPINK7-KLK12 Binding Free Energy

**Experiment ID:** EXP-02  
**Feature ID:** F-02 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** SELECTIVITY TARGET  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment determines the binding free energy (ΔG_bind) of the SPINK7-KLK12 complex and the selectivity ratio ΔΔG_selectivity = ΔG(KLK12) − ΔG(KLK5) using the V2 pipeline. SPINK7 inhibits KLK12 with 11.5-fold weaker affinity than KLK5 (Ki = 1500 nM vs. 132 nM, ΔΔG = +1.5 kcal/mol). This relative measurement is inherently more accurate as a benchmark than absolute ΔG values because systematic computational errors partially cancel in the difference. The experiment tests whether the pipeline can discriminate between cognate and non-cognate protease targets.

---

## 2. Hypothesis

**H₁:** The V2 pipeline will predict that SPINK7 binds KLK12 more weakly than KLK5, with ΔΔG_selectivity (KLK12 − KLK5) within the 95% CI [−1.2, +4.2] kcal/mol.

**H₂:** The absolute ΔG_bind for SPINK7-KLK12 will fall within the range [−12.4, −3.4] kcal/mol (from −7.9 ± 1.96 × 2.3 kcal/mol).

These hypotheses are falsifiable: if the pipeline predicts KLK12 binding to be tighter than KLK5 (ΔΔG < 0 by >2.7 kcal/mol), both hypotheses are rejected.

---

## 3. Background and Rationale

### 3.1 Scientific Context

Azouz et al. (2020) demonstrated that SPINK7 is not a promiscuous serine protease inhibitor. While it potently inhibits KLK5 (Ki = 132 nM), it shows 11.5-fold weaker inhibition of KLK12 (Ki = 1500 nM) and no detectable inhibition of KLK7, KLK11, or KLK13. The ΔΔG_selectivity = +1.5 kcal/mol captures a physically meaningful difference — likely arising from subsite incompatibilities between the SPINK7 RSL and KLK12 active site subsites compared to KLK5.

### 3.2 Pipeline Capability

The pipeline can compute ΔG_bind for SPINK7-KLK12 using the same SMD/Jarzynski and US/WHAM methodology as EXP-01, with a homology model or docking-derived starting structure of the SPINK7-KLK12 complex. The selectivity ΔΔG is computed as the difference: ΔΔG = ΔG(SPINK7-KLK12) − ΔG(SPINK7-KLK5).

### 3.3 What This Reveals

Selectivity validation tests the pipeline's ability to capture subtle energetic differences between closely related protease targets sharing the same inhibitor. A PASS validates the pipeline's sensitivity to subsite-level interactions; a FAIL (predicting equal or reversed selectivity) indicates the pipeline cannot resolve ~1.5 kcal/mol differences — a significant limitation for drug discovery applications.

---

## 4. Experimental Protocol

### 4.1 System Preparation

#### 4.1.1 Structure Acquisition

| Parameter | Value |
|-----------|-------|
| SPINK7 structure | PDB 2LEO (NMR ensemble, chain A, model 1) |
| KLK12 structure | Homology model based on KLK5 (PDB 2PSX) or AlphaFold prediction |
| Complex model | ClusPro or equivalent docking of SPINK7 onto KLK12, analogous to the SPINK7-KLK5 approach |
| Note | No experimental co-crystal structure exists for SPINK7-KLK12 |

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

### 4.2 Umbrella Sampling

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

### 4.3 SMD (Same parameters as EXP-01)

| Parameter | Value | Source |
|-----------|-------|--------|
| Spring constant | 1000 kJ/mol/nm² | `config.py: SMDConfig.spring_constant_kj_mol_nm2` |
| Pulling velocity | 0.001 nm/ps | `config.py: SMDConfig.pulling_velocity_nm_per_ps` |
| Pull distance | 3.0 nm | `config.py: SMDConfig.pull_distance_nm` |
| Replicates | 50 | `config.py: SMDConfig.n_replicates` |
| Save interval | 1.0 ps | `config.py: SMDConfig.save_interval_ps` |

### 4.4 WHAM and MBAR Analysis

Parameters identical to EXP-01 §4.5 and §4.6:
- WHAM tolerance: 10⁻⁶, max iterations: 100,000, bootstraps: 200, bins: 200
- MBAR solver: "robust", tolerance: 10⁻⁷, max iterations: 10,000, bootstraps: 200, bins: 200

### 4.5 Selectivity Calculation

$$\Delta\Delta G_{\text{selectivity}} = \Delta G_{\text{bind}}^{\text{KLK12}} - \Delta G_{\text{bind}}^{\text{KLK5}}$$

Uncertainty propagation:

$$\sigma_{\Delta\Delta G} = \sqrt{\sigma_{\Delta G_{\text{KLK12}}}^2 + \sigma_{\Delta G_{\text{KLK5}}}^2}$$

Where the systematic method errors partially cancel (same force field, same protocol), reducing σ_method from ~2 kcal/mol per system to ~1 kcal/mol for the difference.

---

## 5. Control Conditions

### 5.1 Positive Control

**EXP-01 (SPINK7-KLK5):** The companion experiment. Both must be completed to compute ΔΔG_selectivity. If EXP-01 fails, the selectivity comparison from EXP-02 has limited interpretive value.

### 5.2 Negative Control / Sanity Checks

1. KLK12 structure should be stable during equilibration (backbone RMSD < 5 Å).
2. The SPINK7-KLK12 complex should maintain a canonical protease-inhibitor binding geometry (RSL in active site).
3. ΔG_bind(KLK12) should be negative (SPINK7 does bind KLK12, just more weakly).

---

## 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind (SPINK7-KLK12) | −7.9 kcal/mol | [−12.4, −3.4] |
| ΔΔG_selectivity | +1.5 kcal/mol | [−1.2, +4.2] |
| Direction of selectivity | KLK12 weaker than KLK5 | Must be positive |

### Classification Criteria

- **PASS:** ΔΔG within [−1.2, +4.2] kcal/mol and correct direction (KLK12 weaker)
- **MARGINAL:** Correct direction but ΔΔG outside 95% CI but within 2× CI
- **FAIL:** Wrong direction (KLK12 predicted tighter) or grossly wrong magnitude

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation Implied | Severity |
|-------------|--------------|-------------------|----------|
| **KLK12 homology model error** | Complex dissociates or adopts non-canonical binding mode | KLK12 model quality insufficient for quantitative ΔG | High |
| **Inability to resolve 1.5 kcal/mol ΔΔG** | ΔΔG ≈ 0 or wrong sign | Pipeline lacks sufficient precision for selectivity prediction | High |
| **Systematic cancellation failure** | ΔΔG grossly wrong despite reasonable absolute ΔG values | Error cancellation assumption invalid | Medium |
| **KLK12 active site geometry differences** | Different binding mode vs. KLK5 | Docking procedure failed to identify correct binding pose | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After KLK12 model acquisition | Structural alignment to KLK5 | Cα RMSD < 3 Å for active site residues |
| After docking | SPINK7 RSL in KLK12 active site | P1 residue within 5 Å of catalytic Ser |
| After minimization | IV-1 satisfied | Energy decreased |
| After equilibration | IV-2, IV-3, IV-4, IV-6, IV-7 satisfied | All invariants pass |
| After production MD | Complex remains intact; RMSD stable | No dissociation |
| After umbrella sampling | IV-8: histogram overlap ≥ 10% | Sufficient overlap |
| After WHAM | IV-9: convergence achieved | WHAM converged |
| Selectivity comparison | Both EXP-01 and EXP-02 completed | ΔΔG computable |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
