# EXP-12: LEKTI–KLK5 pH-Dependent Kinetics

**Experiment ID:** EXP-12  
**Feature ID:** F-12 (benchmarks.md)  
**Category:** Kinetic  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment tests whether the V2 pipeline can reproduce the pH-dependent binding kinetics of LEKTI domain 6 with KLK5. Experimentally, LEKTI-KLK5 binding is tight at neutral pH (Ki = 2.35 nM at pH 7.4) but weakens dramatically at acidic pH (reduced inhibition at pH 4.5), with ΔΔG(pH 4.5 − pH 7.5) > 2 kcal/mol (Deraison et al. 2007). The pipeline will run parallel US/WHAM calculations at different pH values (modeled via protonation state assignment) to quantify the pH-dependent shift in ΔG_bind.

---

## 2. Hypothesis

**H₁:** The difference in binding free energy between pH 4.5 and pH 7.5, ΔΔG = ΔG(pH 4.5) − ΔG(pH 7.5), will be > +2 kcal/mol (weaker binding at acidic pH).

**H₂:** The pH sensitivity arises from protonation of histidine residues at the binding interface, disrupting key electrostatic interactions.

---

## 3. Background and Rationale

The pH-dependent regulation of KLK5 by LEKTI is physiologically critical: in the upper epidermis (pH ~5.5), LEKTI inhibits KLK5 to prevent premature desquamation. As pH drops in the outer stratum corneum (pH ~4.5–5.0), LEKTI releases KLK5, enabling corneodesmosome cleavage. This pH switch is mediated by protonation of interfacial histidine residues. Reproducing this pH-dependent behavior validates the pipeline's ability to model electrostatic modulation of binding.

In classical MD, pH effects are modeled by fixing protonation states using pKa prediction tools (e.g., PROPKA or H++) at each target pH. This is a static approximation (no constant-pH MD), but should capture the dominant effect of histidine protonation.

---

## 4. Experimental Protocol

### 4.1 System Preparation at Multiple pH Values

For each pH condition (pH 4.5, 5.5, 6.5, 7.5):

1. Assign protonation states using PROPKA at the target pH.
2. Key residues to monitor: all histidines at the LEKTI-KLK5 interface.
   - At pH 7.5: His residues predominantly neutral (Nε protonated).
   - At pH 4.5: His residues protonated (both Nδ and Nε, net +1 charge).
3. Build each system independently:
   - Force field: AMBER ff14SB (`amber14-all.xml`)
   - Water model: TIP3P (`amber14/tip3p.xml`)
   - Box padding: 1.2 nm, cubic box
   - Ionic strength: 0.15 M NaCl

### 4.2 Equilibration (per pH condition)

- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT equilibration: 500 ps at 310 K, friction 1.0 ps⁻¹, restraints 1000 kJ/mol/nm²
- NPT equilibration: 1000 ps at 310 K, 1.0 atm, barostat interval 25 steps
- Timestep: 2 fs throughout

### 4.3 Umbrella Sampling (per pH condition)

- Reaction coordinate ξ: COM distance, LEKTI D6 ↔ KLK5
- ξ range: 1.5–4.0 nm, spacing 0.05 nm (51 windows)
- Spring constant: 1000 kJ/mol/nm²
- Per-window: 10 ns production (after 200 ps equilibration)
- Pre-positioning velocity: 0.01 nm/ps
- Save interval: 1.0 ps

### 4.4 PMF Reconstruction and ΔG Extraction

- WHAM: tolerance = 10⁻⁶, max_iterations = 100,000, bins = 200
- Bootstrap: n_bootstrap = 200
- ΔG_bind at each pH from PMF well depth + volume correction (C° = 1/1660 Å⁻³)
- Also verify with MBAR (solver = "robust", tolerance = 10⁻⁷)

### 4.5 ΔΔG Calculation

$$\Delta\Delta G = \Delta G_{bind}(\text{pH 4.5}) - \Delta G_{bind}(\text{pH 7.5})$$

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | ΔΔG(pH 4.5 − pH 7.5) |
|---------------|----------------------|
| **PASS** | > +2 kcal/mol (correct direction and magnitude) |
| **MARGINAL** | +1 to +2 kcal/mol (correct direction, weaker than expected) |
| **FAIL** | ≤ 0 kcal/mol (wrong direction) or > +10 kcal/mol (unphysical) |

---

## 5. Control Conditions

### 5.1 Positive Control

**EXP-03 (LEKTI-KLK5 at pH 7.4):** The ΔG_bind at pH 7.5 should be consistent with the Ki = 2.35 nM value measured at neutral pH.

### 5.2 Internal Consistency Control

**pH 6.5 and 5.5 intermediate points:** The pH dependence should be monotonic — binding weakens progressively from pH 7.5 → 6.5 → 5.5 → 4.5.

### 5.3 Negative Control

**BPTI-trypsin pH sensitivity:** BPTI-trypsin binding is relatively pH-insensitive over this range (no histidines at interface). Running pH 4.5 vs 7.5 for BPTI-trypsin should yield |ΔΔG| < 1 kcal/mol.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| ΔΔG(pH 4.5 − pH 7.5) | > +2 kcal/mol | Deraison 2007 |
| ΔG_bind(pH 7.5) | ≈ −11.8 kcal/mol | Ki = 2.35 nM |
| ΔG_bind(pH 4.5) | ≈ −9.8 kcal/mol or weaker | Estimated from ΔΔG |
| pH dependence direction | Weaker binding at lower pH | Physiological mechanism |
| Key protonation change | Interface His residues | PROPKA prediction |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Fixed protonation states** | Static pKa assignment misses coupling | No constant-pH MD | High |
| **No interface His residues** | No protonation state change at interface | LEKTI-KLK5 model issues | High |
| **PROPKA pKa error** | Wrong residues protonated | Empirical pKa prediction | Medium |
| **Insufficient ΔΔG** | Correct direction but < 2 kcal/mol | Charge redistribution not captured | Medium |
| **4× computational cost** | 4 pH conditions × 510 ns each | Resource constraint | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| PROPKA at pH 7.5 | All His neutral (or one protonated) | Consistent with NMR |
| PROPKA at pH 4.5 | All interface His protonated (+1) | pKa > 4.5 |
| Equilibrated structures | Stable RMSD at each pH | Backbone RMSD < 2.5 Å |
| Histogram overlap | All 51 windows adequately sampled per pH | >10% overlap |
| PMF quality per pH | Smooth PMFs with clear minima | Bootstrap error < 2 kcal/mol |
| Monotonic pH trend | ΔG weakens: pH 7.5 > 6.5 > 5.5 > 4.5 | Consistent trend |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
