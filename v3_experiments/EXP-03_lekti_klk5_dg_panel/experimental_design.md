# EXP-03: LEKTI-KLK5 Binding Free Energy Panel

**Experiment ID:** EXP-03  
**Feature ID:** F-03 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** CROSS-VALIDATION  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment determines the binding free energy of a representative LEKTI (SPINK5) single Kazal domain–KLK5 complex using US/WHAM, benchmarked against the most precise single-domain Ki value: rLEKTI(1–6) Ki = 2.35 ± 0.22 nM (ΔG = −11.8 kcal/mol) from Borgono et al. (2007). LEKTI domains share ~35% sequence identity with SPINK7 and use the same canonical Kazal inhibitory mechanism. Testing a related Kazal-KLK5 interaction validates whether the pipeline generalizes across Kazal family members and whether systematic biases observed for SPINK7-KLK5 (EXP-01) are system-specific or method-inherent.

---

## 2. Hypothesis

**H₁:** The V2 pipeline's US/WHAM estimate of ΔG_bind for an isolated LEKTI Kazal domain (domain 6 or domains 1–6 fragment) against KLK5 will fall within the 95% CI [−16.2, −7.4] kcal/mol.

**H₂:** The computed ΔG_bind will be more negative (tighter binding) than the SPINK7-KLK5 prediction (EXP-01), consistent with the experimental observation that LEKTI fragments bind KLK5 more tightly (Ki = 2.35 nM vs. 132 nM, ΔΔG ≈ 2.4 kcal/mol).

---

## 3. Background and Rationale

### 3.1 Scientific Context

LEKTI is a 15-domain Kazal-type inhibitor that regulates KLK5 activity in the skin. Three independent laboratories have measured LEKTI-KLK5 inhibition (Borgono et al. 2007, Deraison et al. 2007, Egelrud et al. 2005), providing multi-group cross-validation. Single-domain Ki values range from 2.35–118.7 nM, while multi-domain fragments achieve sub-picomolar apparent KD by avidity effects. For MD simulation, single-domain Ki values are the appropriate targets.

### 3.2 What This Reveals

Cross-validation with a structurally related but distinct inhibitor-protease pair tests pipeline transferability. If the pipeline performs well on LEKTI-KLK5 but poorly on SPINK7-KLK5, the issue is SPINK7-specific (likely model quality). If both fail, the issue is methodological.

---

## 4. Experimental Protocol

### 4.1 System Preparation

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

### 4.2 Minimization and Equilibration

Same parameters as EXP-01 §4.1.3 and §4.1.4:
- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT: 500 ps at 310 K, friction 1.0 ps⁻¹, timestep 2 fs
- NPT: 1000 ps at 310 K and 1.0 atm, barostat interval 25
- Restraint: 1000 kJ/mol/nm², save interval 10 ps

### 4.3 Umbrella Sampling

| Parameter | Value | Source |
|-----------|-------|--------|
| ξ range | 1.5–4.0 nm | `config.py: UmbrellaConfig` |
| Window spacing | 0.05 nm (51 windows) | `config.py: UmbrellaConfig.window_spacing_nm` |
| Spring constant | 1000 kJ/mol/nm² | `config.py: UmbrellaConfig.spring_constant_kj_mol_nm2` |
| Per-window duration | 10.0 ns | `config.py: UmbrellaConfig.per_window_duration_ns` |
| Pre-positioning velocity | 0.01 nm/ps | `config.py: UmbrellaConfig.pre_position_velocity_nm_per_ps` |
| Equilibration per window | 200 ps | `config.py: UmbrellaConfig.equilibration_duration_ps` |
| Save interval | 1.0 ps | `config.py: UmbrellaConfig.save_interval_ps` |

### 4.4 WHAM and MBAR Analysis

- WHAM: tolerance 10⁻⁶, max iterations 100,000, bootstraps 200, bins 200
- MBAR: solver "robust", tolerance 10⁻⁷, max iterations 10,000, bootstraps 200, bins 200

### 4.5 Statistical Comparison

Using the §25 framework with:
- σ_exp = 0.06 kcal/mol (from 0.22 nM on 2.35 nM)
- σ_comp = from bootstrap
- σ_method = 2.0 kcal/mol (protein-protein US/WHAM, §25.4)

---

## 5. Control Conditions

### 5.1 Positive Control

**EXP-01 (SPINK7-KLK5):** Same protease, related inhibitor. Comparison reveals whether discrepancies are system-specific or method-inherent.

### 5.2 Negative Control

1. All physical validity invariants (IV-1 through IV-9) must be satisfied.
2. LEKTI domain should maintain its Kazal fold during MD (backbone RMSD < 3 Å).
3. ΔG should be negative (binding observed experimentally).

---

## 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind (LEKTI D6–KLK5) | −11.8 kcal/mol | [−16.2, −7.4] |
| Ranking vs. SPINK7-KLK5 | More favorable ΔG | LEKTI binds tighter |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation Implied | Severity |
|-------------|--------------|-------------------|----------|
| **LEKTI domain model error** | Non-canonical binding mode | Homology model quality | High |
| **Multi-domain avidity not captured** | ΔG weaker than single-domain Ki predicts | Single-domain simulation appropriate | Low (expected) |
| **Ranking inversion** | LEKTI predicted weaker than SPINK7 | Pipeline cannot rank related inhibitors | High |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After LEKTI model | Kazal fold integrity; disulfides present | Valid Kazal topology |
| After docking | RSL in KLK5 active site | Canonical binding mode |
| After equilibration | IV-2 through IV-7 satisfied | All pass |
| After US/WHAM | IV-8, IV-9 satisfied; smooth PMF | Converged result |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
