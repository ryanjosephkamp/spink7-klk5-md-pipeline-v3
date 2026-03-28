# Consolidated Experimental Benchmarks — V3 Part 1

**Document Class:** V3 Literature Review — Post-Review Benchmarks Compilation  
**Domain:** Computational Biophysics — SPINK7-KLK5 Molecular Dynamics Pipeline  
**Date:** 2026-03-22  
**Version:** 1.0  
**Input Documents:** 28 source reports (`v3_literature/source_reports/*.md`), `literature_review_progress.csv`

---

## Table of Contents

1. [Overview](#1-overview)
2. [Feature Summary Table](#2-feature-summary-table)
3. [Category I — Thermodynamic Features](#3-category-i--thermodynamic-features)
4. [Category II — Kinetic Features](#4-category-ii--kinetic-features)
5. [Category III — Structural Features](#5-category-iii--structural-features)
6. [Category IV — Dynamic Features](#6-category-iv--dynamic-features)
7. [Category V — Biophysical Features](#7-category-v--biophysical-features)
8. [Category VI — Mutational Features](#8-category-vi--mutational-features)
9. [Summary Statistics](#9-summary-statistics)
10. [Coverage Assessment](#10-coverage-assessment)
11. [Key PDB Codes and Structural References](#11-key-pdb-codes-and-structural-references)
12. [References](#12-references)

---

## 1. Overview

This document consolidates all experimentally measurable features identified during the Step 3 literature review of the V3 Part 1 pipeline project. A total of **33 unique features** were extracted from **28 source publications**, organized across **six categories** (thermodynamic, kinetic, structural, dynamic, biophysical, mutational). Each feature includes:

- **(a)** A thorough background explanation of the feature, its physical significance, and relevance to the SPINK7-KLK5 pipeline.
- **(b)** Literature source(s) with full citations and specific table/figure/page references.
- **(c)** Experimental value(s) with exact values, units, and uncertainty where reported.
- **(d)** 95% confidence interval for pipeline prediction, with scientific justification based on the §25 framework.

Features are deduplicated — when the same measurement appears in multiple sources, entries are consolidated with all supporting citations. The most precise experimental value is designated as the primary benchmark.

### Systems Hierarchy

| System | Relationship to SPINK7-KLK5 | Role in Benchmarking |
|--------|------------------------------|---------------------|
| SPINK7-KLK5 | Primary target | Direct validation |
| SPINK7-KLK12 | Same inhibitor, different protease | Selectivity validation |
| LEKTI-KLK5 | Same protease, related Kazal inhibitor | Cross-validation |
| SPINK1-trypsin | Same Kazal family, different system | Kazal-family transferability |
| PSTI-chymotrypsinogen | Kazal-type, co-crystal available | Structural reference |
| BPTI-trypsin | Protease-inhibitor gold standard | Method validation |
| Barnase-barstar | Ultra-tight PPI reference | PPI energetics methodology |
| SH3-p41 peptide | Computational methods benchmark | PMF/alchemical calibration |

---

## 2. Feature Summary Table

| ID | Feature | Category | Primary System | Primary Source(s) | Primary Value | Pipeline Method |
|----|---------|----------|---------------|-------------------|---------------|----------------|
| F-01 | SPINK7-KLK5 ΔG_bind | Thermodynamic | SPINK7-KLK5 | S-03 | Ki = 132 nM (ΔG ≈ −9.4 kcal/mol) | SMD/Jarzynski, US/WHAM |
| F-02 | SPINK7-KLK12 ΔG_bind | Thermodynamic | SPINK7-KLK12 | S-03 | Ki = 1500 nM (ΔG ≈ −7.9 kcal/mol) | SMD/Jarzynski, US/WHAM |
| F-03 | LEKTI-KLK5 ΔG_bind panel | Thermodynamic | LEKTI-KLK5 | S-08, S-09, S-11 | Ki = 2.35–118.7 nM | US/WHAM, FEP |
| F-04 | BPTI-trypsin ΔG_bind | Thermodynamic | BPTI-trypsin | S-16, S-17 | Kd = 6 × 10⁻¹⁴ M (ΔG ≈ −18 kcal/mol) | US/WHAM, SMD |
| F-05 | PSTI-chymotrypsinogen ΔG_bind | Thermodynamic | PSTI-chymotrypsinogen | S-19 | Ki = 1.6 × 10⁻¹¹ M (ΔG ≈ −14.7 kcal/mol) | US/WHAM |
| F-06 | SPINK1-trypsin ΔG_bind | Thermodynamic | SPINK1-trypsin | S-20 | Ki ≈ 7.6 nM (ΔG ≈ −11.1 kcal/mol) | US/WHAM |
| F-07 | P1 residue energetic contribution | Thermodynamic | Canonical inhibitors | S-14, S-17 | ~70% of total ΔG_bind | Energy decomposition |
| F-08 | Interfacial H-bond energy | Thermodynamic | Canonical inhibitors | S-14 | ~1.5 kcal/mol per H-bond | H-bond energy analysis |
| F-09 | BPTI-trypsin kon/koff | Kinetic | BPTI-trypsin | S-16, S-17 | kon = 1.1 × 10⁶ M⁻¹s⁻¹; koff = 6.6 × 10⁻⁸ s⁻¹ | MSM, PMF-derived |
| F-10 | BPTI-trypsin Ea(association) | Kinetic | BPTI-trypsin | S-16 | Ea = 10.5 kcal/mol | PMF barrier height |
| F-11 | KLK5 substrate Km/kcat | Kinetic | KLK5 | S-06 | Km = 0.19–1.69 mM; kcat = 4.56–196.8 min⁻¹ | Qualitative validation |
| F-12 | LEKTI-KLK5 pH-dependent kinetics | Kinetic | LEKTI-KLK5 | S-09 | ka/kd vary 10–100× between pH 4.5 and 7.5 | Multi-pH simulations |
| F-13 | Barnase-barstar kon/koff | Kinetic | Barnase-barstar | S-22 | kon = 3.7 × 10⁸ M⁻¹s⁻¹; koff = 3.7 × 10⁻⁶ s⁻¹ | MSM reference |
| F-14 | Binding loop backbone geometry | Structural | Canonical inhibitors | S-13, S-14, S-19 | φ/ψ ranges for P3–P1 | Trajectory dihedral analysis |
| F-15 | Catalytic contact distance | Structural | Canonical complexes | S-12, S-13 | P1 C=O to Ser195 Oγ ≈ 2.7 Å | Distance monitoring |
| F-16 | Protease-inhibitor BSA | Structural | 75 PPI dataset | S-24, S-14, S-13 | 1530 ± 170 Å² (protease-inhibitor mean) | SASA calculation |
| F-17 | Interface H-bond count | Structural | Protease-inhibitor | S-24, S-14 | ~10 H-bonds (protease-inhibitor mean) | H-bond analysis |
| F-18 | Interface water molecules | Structural | Protease-inhibitor | S-24 | ~15 waters per interface | Solvation analysis |
| F-19 | Interface non-polar fraction | Structural | Protease-inhibitor | S-24 | 61% non-polar | Contact composition |
| F-20 | Binding loop structural conservation | Structural | Kazal family | S-19, S-13, S-14 | rmsd 0.25–0.35 Å across family | Cross-system RMSD |
| F-21 | KLK5 subsite specificity | Structural | KLK5 | S-07 | P1: Arg exclusive; P2: Ser; P3: Met/Phe/Tyr; P4: Gly | Binding pose analysis |
| F-22 | SPINK7-KLK5 complex geometry | Structural | SPINK7-KLK5 | S-03 | RSL ≈ 2.5 Å from catalytic Ser | Complex equilibration |
| F-23 | BPTI high-resolution structure | Structural | BPTI | S-18 | 0.94 Å X-ray (PDB 4PTI/5PTI) | Structural RMSD |
| F-24 | BPTI amide H/D exchange | Dynamic | BPTI | S-18 | 11 protected amide H after 3 months D₂O | N-H order parameters |
| F-25 | BPTI conformational variability | Dynamic | BPTI | S-18 | Form I/II rmsd = 0.40 Å (main chain) | Trajectory ensemble width |
| F-26 | BSA-ΔG empirical correlation | Biophysical | PPI dataset | S-23 | ΔG ≈ 0.025 × BSA (kcal/mol per Å²); r = 0.54 | Cross-check |
| F-27 | Interface close-packing quality | Biophysical | Protease-inhibitor | S-24 | V/Vo = 1.00 ± 0.01 | Voronoi packing analysis |
| F-28 | Inhibitor scaffold energy | Biophysical | BPTI, Kazal | S-13 | −33 kJ/mol (−7.9 kcal/mol) | Scaffold vs. loop decomposition |
| F-29 | SH3-p41 ΔG_bind (method validation) | Biophysical | SH3-p41 | S-28, S-25 | Expt: −7.99 kcal/mol; PMF: −7.8; Alch: −7.7 | PMF via REMD-US |
| F-30 | BPTI alanine scan ΔΔG panel | Mutational | BPTI-trypsin | S-17 | 15 mutants; K15A ΔΔG ≈ 10 kcal/mol | Alchemical FEP |
| F-31 | BPTI disulfide ablation | Mutational | BPTI-trypsin | S-16 | C14-C38 reduction: Kd increases ~10⁵-fold | Modified topology simulation |
| F-32 | SPINK1 N34S functional impact | Mutational | SPINK1-trypsin | S-20 | No change in Ki (negative control) | FEP (ΔΔG ≈ 0) |
| F-33 | Barnase-barstar DMC coupling | Mutational | Barnase-barstar | S-22 | ~45 coupling pairs; max 7.0 kcal/mol | FEP thermodynamic cycles |

---

## 3. Category I — Thermodynamic Features

### F-01: SPINK7-KLK5 Binding Free Energy

**Status: PRIMARY PIPELINE TARGET**

**(a) Feature Description**

The binding free energy of the SPINK7-KLK5 complex is the single most important benchmark for the V2 pipeline. SPINK7 (Serine Peptidase Inhibitor, Kazal Type 7) is a small (~6 kDa) secreted Kazal-type inhibitor that stoichiometrically inhibits KLK5, a trypsin-like serine protease responsible for desquamation regulation in the esophageal epithelium. In Eosinophilic Esophagitis (EoE), IL-13-driven down-regulation of SPINK7 unleashes KLK5 proteolytic activity, degrading Desmoglein-1 and compromising the epithelial barrier. The binding affinity, expressed as Ki and converted to ΔG via ΔG = RT ln(Ki), quantifies the thermodynamic stability of the inhibitory complex and is the primary validation target for all free energy calculations in the pipeline.

The measurement was performed using the Morrison tight-binding inhibitor equation, which is appropriate because Ki is comparable to or smaller than the enzyme concentration used in the assay. This provides a more reliable Ki estimate than the classical Michaelis-Menten approach for potent inhibitors.

**(b) Literature Source(s)**

- **Azouz et al. (2020)** — *Science Translational Medicine*, 12(545), eaaz7773. **Fig. 1B.** Dose-response curve of SPINK7 against KLK5 using fluorogenic substrate; Morrison equation fit. Three independent experiments performed in duplicate.
- **Source Report:** `source_reports/azouz_et_al_2020_report.md` (S-03, Feature 1)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Ki | 132 nM |
| ΔG_bind (298 K) | −9.4 kcal/mol (calculated: ΔG = RT ln(Ki) = 0.592 × ln(132 × 10⁻⁹)) |
| Method | Morrison tight-binding equation |
| Substrate | Fluorogenic (BOC-VPR-AMC typical for trypsin-like proteases) |
| Replicates | 3 independent experiments, each in duplicate |
| Uncertainty | Not explicitly reported; estimated ±2–3× from typical Morrison Ki variability |

**(d) 95% Confidence Interval for Pipeline Prediction**

$$\text{CI}_{95\%}(\Delta G) = -9.4 \pm 1.96 \times \sqrt{\sigma_{\text{exp}}^2 + \sigma_{\text{comp}}^2 + \sigma_{\text{method}}^2}$$

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.6 kcal/mol | Morrison Ki assay 2–3× variability → 50–350 nM range → ΔG = −10.0 to −8.8 |
| σ_comp | 1.0 kcal/mol | Expected from 50 SMD replicates or 51 umbrella windows (S-25, S-28) |
| σ_method | 2.0 kcal/mol | US/WHAM for protein-protein binding (§25.4; S-25: 1–2.5 kcal/mol RMSE) |
| **σ_combined** | **2.4 kcal/mol** | √(0.6² + 1.0² + 2.0²) = 2.37 |
| **95% CI** | **−14.0 to −4.8 kcal/mol** | −9.4 ± (1.96 × 2.37) = −9.4 ± 4.6 |

**PASS criterion:** Pipeline ΔG_bind within [−14.0, −4.8] kcal/mol  
**MARGINAL criterion:** Pipeline ΔG_bind within [−18.7, −0.1] kcal/mol  
**FAIL criterion:** Pipeline ΔG_bind outside [−18.7, −0.1] kcal/mol

---

### F-02: SPINK7-KLK12 Binding Free Energy

**Status: SELECTIVITY TARGET**

**(a) Feature Description**

SPINK7 inhibits KLK12 with 11.5-fold weaker affinity than KLK5, establishing that SPINK7 is not a promiscuous serine protease inhibitor but possesses target selectivity. The selectivity ratio (ΔΔG_selectivity = ΔG_KLK12 − ΔG_KLK5 = −7.9 − (−9.4) = +1.5 kcal/mol) provides a relative benchmark that is inherently more accurate than absolute ΔG values because systematic errors in the computational method partially cancel in the difference. SPINK7 shows no detectable inhibition of KLK7, KLK11, or KLK13.

**(b) Literature Source(s)**

- **Azouz et al. (2020)** — *Sci. Transl. Med.*, 12(545), eaaz7773. **Fig. 1B, fig. S1B.** Same experimental setup as F-01.
- **Source Report:** `source_reports/azouz_et_al_2020_report.md` (S-03, Features 2 and 3)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Ki (SPINK7-KLK12) | 1500 nM |
| ΔG_bind (298 K) | −7.9 kcal/mol |
| ΔΔG_selectivity (KLK12 − KLK5) | +1.5 kcal/mol |
| Selectivity ratio | 11.5× (132 vs 1500 nM) |
| Non-inhibited proteases | KLK7, KLK11, KLK13 (no detectable inhibition) |

**(d) 95% CI for Pipeline Prediction**

For the selectivity ratio (ΔΔG), systematic errors partially cancel:

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.8 kcal/mol | Propagated from both Ki measurements |
| σ_comp | 0.5 kcal/mol | Partial cancellation of systematic error in ΔΔG |
| σ_method | 1.0 kcal/mol | Relative FEP accuracy better than absolute (S-25) |
| **σ_combined** | **1.4 kcal/mol** | √(0.8² + 0.5² + 1.0²) |
| **95% CI for ΔΔG** | **−1.2 to +4.2 kcal/mol** | +1.5 ± (1.96 × 1.4) = +1.5 ± 2.7 |

**PASS criterion:** Pipeline predicts KLK12 binding weaker than KLK5, with ΔΔG within [−1.2, +4.2] kcal/mol

---

### F-03: LEKTI-KLK5 Binding Free Energy Panel

**Status: CROSS-VALIDATION**

**(a) Feature Description**

LEKTI (Lympho-Epithelial Kazal-Type Inhibitor), encoded by the SPINK5 gene, is a large multi-domain (15 Kazal domains) serine protease inhibitor that regulates KLK5 activity in the skin. Multiple recombinant LEKTI fragments have been characterized for KLK5 inhibition, providing a panel of related Kazal-type binding affinities that serve as cross-validation targets. SPINK7 shares 35% sequence identity with SPINK5 domain 15, and both use the canonical Kazal inhibitory mechanism. These data establish the expected affinity range for Kazal-KLK5 interactions (Ki ≈ 2–120 nM) and demonstrate that multi-domain LEKTI fragments can achieve sub-picomolar apparent KD by avidity effects (SPR).

Three independent laboratories (Borgono/Diamandis, Deraison/Hovnanian, Egelrud/Brattsand) have measured LEKTI-KLK5 inhibition, providing multi-group cross-validation.

**(b) Literature Source(s)**

- **Borgono et al. (2007)** — *J. Biol. Chem.*, 282(6), 3640–3652. **Table II.** Ki values for four recombinant LEKTI fragments against KLK5, mixed inhibition.
- **Deraison et al. (2007)** — *Mol. Biol. Cell*, 18, 3607–3619. **Table 1, Fig. 2, Fig. 3.** Ki and SPR KD for four LEKTI domain constructs; pH-dependent ka/kd.
- **Egelrud et al. (2005)** — *Br. J. Dermatol.*, 153, 1200–1203. **Fig. 1.** LD6-KLK5 IC₅₀.
- **Source Reports:** `borgono_et_al_2007_report.md` (S-08), `deraison_et_al_2007_report.md` (S-09), `egelrud_et_al_2005_report.md` (S-11)

**(c) Experimental Value(s)**

**From S-08 (Borgono et al. 2007) — Functional Inhibition:**

| LEKTI Fragment | Ki (nM) | Mechanism | ΔG_bind (298 K) |
|----------------|---------|-----------|-----------------|
| rLEKTI(1–6) | 2.35 ± 0.22 | Mixed | −11.8 kcal/mol |
| rLEKTI(6–9') | 4.68 ± 0.66 | Mixed | −11.4 kcal/mol |
| rLEKTI(9–12) | 2.75 ± 0.24 | Mixed | −11.7 kcal/mol |
| rLEKTI(12–15) | 21.80 ± 2.40 | Mixed | −10.5 kcal/mol |

**From S-09 (Deraison et al. 2007) — Functional Inhibition + SPR:**

| LEKTI Fragment | Ki (nM) | KD by SPR (M) | ΔG from KD (298 K) |
|----------------|---------|---------------|---------------------|
| D5 | 32.8 | 9.3 × 10⁻¹⁰ | −12.3 kcal/mol |
| D6 | 83.3 | 3.6 × 10⁻⁹ | −11.5 kcal/mol |
| D8–D11 | 3.7 | 1.1 × 10⁻¹² | −16.3 kcal/mol |
| D9–D15 | 118.7 | 9.5 × 10⁻¹⁰ | −12.3 kcal/mol |

**From S-11 (Egelrud et al. 2005):**

| LEKTI Fragment | IC₅₀ (nM) | Conditions |
|----------------|-----------|------------|
| LD6 (native) | ~60 | 80 mM Tris pH 8.0, 80 mM NaCl |

**Note:** The large discrepancy between Ki and SPR KD for multi-domain fragments (e.g., D8-D11: Ki = 3.7 nM vs KD = 1.1 pM) reflects avidity effects in SPR where multiple domains can engage the immobilized protease surface simultaneously. Single-domain Ki values are the most appropriate benchmarks for MD simulation of isolated Kazal domain binding.

**(d) 95% CI for Pipeline Prediction**

Using the most precise single-domain value, rLEKTI(1-6) Ki = 2.35 ± 0.22 nM (ΔG = −11.8 kcal/mol):

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.06 kcal/mol | Propagated from 0.22 nM uncertainty on 2.35 nM |
| σ_comp | 1.0 kcal/mol | Expected from US/WHAM convergence |
| σ_method | 2.0 kcal/mol | Protein-protein US/WHAM (§25.4) |
| **σ_combined** | **2.2 kcal/mol** | √(0.06² + 1.0² + 2.0²) |
| **95% CI** | **−16.2 to −7.4 kcal/mol** | −11.8 ± (1.96 × 2.2) |

---

### F-04: BPTI-Trypsin Binding Free Energy

**Status: GOLD STANDARD REFERENCE**

**(a) Feature Description**

The Bovine Pancreatic Trypsin Inhibitor (BPTI)–trypsin complex is the most extensively characterized protease-inhibitor system in structural biology. BPTI is a 58-residue Kunitz-type inhibitor (not Kazal-type) but shares the canonical serine protease inhibition mechanism with SPINK7: the reactive site loop inserts into the protease active site in a substrate-like manner, and the resulting complex is kinetically trapped. The BPTI-trypsin Kd of 6 × 10⁻¹⁴ M (ΔG ≈ −18 kcal/mol) is the gold-standard benchmark for protease-inhibitor free energy calculations. While significantly tighter than the SPINK7-KLK5 interaction (ΔG ≈ −9.4 kcal/mol), it provides a secondary validation point on the affinity scale and tests whether the pipeline can reproduce ultra-tight binding.

**(b) Literature Source(s)**

- **Vincent & Lazdunski (1972)** — *Biochemistry*, 11(16), 2967–2977. **Table I, pp. 2970–2971.** Kd determination from equilibrium inhibition; pH dependence (His46 pK 7.05).
- **Castro & Anderson (1996)** — *Biochemistry*, 35(35), 11435–11446. **Table 1, Table 2.** Confirms Ki = 5 × 10⁻¹⁴ M with full kinetic/thermodynamic characterization.
- **Source Reports:** `vincent_lazdunski_1972_report.md` (S-16), `castro_anderson_1996_report.md` (S-17)

**(c) Experimental Value(s)**

| Parameter | S-16 Value | S-17 Value | Consensus |
|-----------|-----------|-----------|-----------|
| Kd / Ki | 6 × 10⁻¹⁴ M | 5 × 10⁻¹⁴ M | ~6 × 10⁻¹⁴ M |
| ΔG_bind (298 K) | −18.0 kcal/mol | −18.1 kcal/mol | −18.0 kcal/mol |
| Conditions (S-16) | 50 mM Tris, 50 mM CaCl₂, 0.1 M NaCl, pH 8.0, 25°C | — | — |
| Conditions (S-17) | — | 50 mM Tris-HCl, 20 mM CaCl₂, 0.005% Triton X-100, pH 8.2, 22°C | — |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.5 kcal/mol | Two independent labs agree within 0.1 kcal/mol; estimated from typical assay variability |
| σ_comp | 1.5 kcal/mol | Ultra-tight binding requires extensive sampling |
| σ_method | 3.0 kcal/mol | Protein-protein US at −18 kcal/mol is extremely challenging (S-25, S-26) |
| **σ_combined** | **3.4 kcal/mol** | √(0.5² + 1.5² + 3.0²) |
| **95% CI** | **−24.7 to −11.3 kcal/mol** | −18.0 ± (1.96 × 3.4) |

---

### F-05: PSTI-Chymotrypsinogen Binding Free Energy

**Status: KAZAL ANALOG — CLOSEST STRUCTURAL REFERENCE**

**(a) Feature Description**

Pancreatic Secretory Trypsin Inhibitor (PSTI, also known as SPINK1) is a Kazal-type inhibitor belonging to the same protein family as SPINK7. The crystal structure of two recombinant PSTI variants (PSTI3 and PSTI4) in complex with bovine chymotrypsinogen A at 2.3 Å resolution provides the most directly relevant structural and thermodynamic reference for the SPINK7-KLK5 system, as both the inhibitor (Kazal-type) and protease (serine protease) are from the same mechanistic class. The Ki values (1.6–2.4 × 10⁻¹¹ M) establish the expected affinity range for a well-optimized Kazal-serine protease interaction.

**(b) Literature Source(s)**

- **Hecht et al. (1991)** — *J. Mol. Biol.*, 220, 711–722. **Table 2, Section 3.** Ki values from Szardenings 1989 thesis; crystal structure contacts.
- **Source Report:** `hecht_et_al_1991_report.md` (S-19)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Ki (PSTI3-chymotrypsin) | 1.6 × 10⁻¹¹ M |
| Ki (PSTI4-chymotrypsin) | 2.4 × 10⁻¹¹ M |
| ΔG_bind (298 K, PSTI3) | −14.7 kcal/mol |
| ΔG_bind (298 K, PSTI4) | −14.5 kcal/mol |
| Structure resolution | 2.3 Å |
| Binding loop rmsd vs. other Kazal | 0.25–0.35 Å |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.3 kcal/mol | Two variants agree within 0.2 kcal/mol |
| σ_comp | 1.0 kcal/mol | Co-crystal starting structure reduces sampling uncertainty |
| σ_method | 2.5 kcal/mol | Tight Kazal-protease binding, US/WHAM |
| **σ_combined** | **2.7 kcal/mol** | √(0.3² + 1.0² + 2.5²) |
| **95% CI** | **−20.0 to −9.4 kcal/mol** | −14.7 ± (1.96 × 2.7) |

---

### F-06: SPINK1-Trypsin Binding Free Energy

**Status: KAZAL FAMILY BENCHMARK**

**(a) Feature Description**

SPINK1 (identical protein to PSTI, but measured against trypsin rather than chymotrypsinogen) provides a Kazal-type inhibitor–trypsin-like protease benchmark. The Ki ≈ 7.6 nM corresponds to ΔG ≈ −11.1 kcal/mol, intermediate between the SPINK7-KLK5 value (−9.4 kcal/mol) and the PSTI-chymotrypsinogen value (−14.7 kcal/mol). This positions SPINK1-trypsin as a calibration point on the Kazal affinity ladder. Importantly, the clinically significant N34S mutation does NOT alter Ki, providing a negative control for FEP validation (see F-32).

**(b) Literature Source(s)**

- **Kuwata et al. (2002)** — *J. Gastroenterol.*, 37, 928–934. **Table 2, p. 930.** Green & Work method for Ki.
- **Source Report:** `kuwata_et_al_2002_report.md` (S-20)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Ki (recombinant wt SPINK1) | 7.6 × 10⁻⁹ M |
| Ki (natural PSTI) | 8.7 × 10⁻⁹ M |
| ΔG_bind (310 K) | −11.1 kcal/mol (recombinant) |
| Conditions | 0.1 M Tris-HCl pH 8.0, 0.02 M CaCl₂, 0.01% Triton X-100, 37°C |
| Method | Green & Work method |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.5 kcal/mol | Single method; recombinant vs natural differ by ~0.1 kcal/mol |
| σ_comp | 1.0 kcal/mol | Standard US/WHAM uncertainty |
| σ_method | 2.0 kcal/mol | Protein-protein US/WHAM (§25.4) |
| **σ_combined** | **2.3 kcal/mol** | √(0.5² + 1.0² + 2.0²) |
| **95% CI** | **−15.6 to −6.6 kcal/mol** | −11.1 ± (1.96 × 2.3) |

---

### F-07: P1 Residue Energetic Contribution to Binding

**(a) Feature Description**

In canonical serine protease-inhibitor complexes, the P1 residue (the amino acid immediately N-terminal to the scissile bond) is the dominant contributor to binding affinity. Structural and mutagenesis studies show that P1 accounts for approximately 70% of the total association energy and ~50% of the interface contact surface area. For trypsin-like proteases (including KLK5), P1 is typically Arg or Lys, which forms a salt bridge with Asp189 at the base of the S1 specificity pocket. This feature tests the pipeline's ability to decompose per-residue energetic contributions at the interface and identifies the P1 hot spot.

**(b) Literature Source(s)**

- **Krowarsch et al. (2003)** — *Cell. Mol. Life Sci.*, 60, 2427–2444. **Section 5.1, Table 3.** Compilation of P1 energetic contributions across multiple inhibitor families.
- **Castro & Anderson (1996)** — *Biochemistry*, 35(35), 11435–11446. **Table 2.** BPTI K15A (P1 → Ala) causes ΔΔG ≈ 10 kcal/mol — quantitative confirmation of P1 dominance.
- **Source Reports:** `krowarsch_et_al_2003_report.md` (S-14, Feature 4), `castro_anderson_1996_report.md` (S-17)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| P1 contribution to ΔG_bind | ~70% of total |
| P1 contribution to contact area | ~50% of interface BSA |
| BPTI K15A ΔΔG | ≈ 10 kcal/mol (confirmatory) |
| P1 residue for KLK5 inhibitors | Arg (exclusive selectivity, S-07) |

**(d) 95% CI for Pipeline Prediction**

This is a semi-quantitative structural/energetic feature. The pipeline should predict:
- The P1 residue contributes the largest per-residue ΔG_bind.
- Energy decomposition: P1 accounts for 50–90% of total binding energy.
- PASS: P1 is identified as the dominant hot-spot residue with >40% of total binding energy.
- FAIL: P1 is not the top contributor, or contributes <30%.

---

### F-08: Interfacial Hydrogen Bond Energy

**(a) Feature Description**

Individual hydrogen bonds at protease-inhibitor interfaces contribute approximately 1.5 kcal/mol each to the binding free energy. This value, averaged from the Laskowski additivity framework, provides a per-interaction energy benchmark that can be tested by computational H-bond energy decomposition or by comparing total H-bond count × 1.5 kcal/mol with the total ΔG_bind.

**(b) Literature Source(s)**

- **Krowarsch et al. (2003)** — *Cell. Mol. Life Sci.*, 60, 2427–2444. **Section 5.2, pp. 2436–2437.** Additivity framework for protease-inhibitor binding.
- **Source Report:** `krowarsch_et_al_2003_report.md` (S-14, Feature 5)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Mean H-bond energy at interface | ~1.5 kcal/mol per H-bond |
| Typical interface H-bond count | 8–15 (S-14, S-24) |
| Expected total H-bond contribution | 12–22.5 kcal/mol |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.5 kcal/mol per bond | Variation across different interface H-bonds |
| σ_comp | 0.3 kcal/mol per bond | From H-bond occupancy fluctuations in MD |
| σ_method | 0.5 kcal/mol per bond | Force field H-bond energy accuracy |
| **95% CI** | **0.7 to 2.3 kcal/mol per H-bond** | 1.5 ± (1.96 × 0.76) |

---

## 4. Category II — Kinetic Features

### F-09: BPTI-Trypsin Association and Dissociation Rate Constants

**(a) Feature Description**

The BPTI-trypsin association kinetics (kon = 1.1 × 10⁶ M⁻¹s⁻¹, koff = 6.6 × 10⁻⁸ s⁻¹) define the kinetic landscape of the gold-standard protease-inhibitor system. The association rate is 3–4 orders of magnitude below the diffusion limit (~10⁹–10¹⁰ M⁻¹s⁻¹), indicating significant orientational and conformational barriers. The extremely slow koff (t₁/₂ ≈ 120 days) reflects the essentially irreversible nature of the complex on biological timescales. These rate constants can, in principle, be estimated from Markov State Model (MSM) analysis of long unbiased trajectories or inferred from the PMF barrier heights along the binding/unbinding coordinate.

**(b) Literature Source(s)**

- **Vincent & Lazdunski (1972)** — *Biochemistry*, 11(16), 2967–2977. **Table I, pp. 2970–2971.** Temperature-dependent kon and koff measurements.
- **Castro & Anderson (1996)** — *Biochemistry*, 35(35), 11435–11446. **Table 1.** Confirms kon = (9.9 ± 2.5) × 10⁵ M⁻¹s⁻¹, koff = 5 × 10⁻⁸ s⁻¹.
- **Source Reports:** `vincent_lazdunski_1972_report.md` (S-16, Feature 2), `castro_anderson_1996_report.md` (S-17, Feature 1)

**(c) Experimental Value(s)**

| Parameter | S-16 Value | S-17 Value | Consensus |
|-----------|-----------|-----------|-----------|
| kon (M⁻¹s⁻¹) | 1.1 × 10⁶ | (9.9 ± 2.5) × 10⁵ | ~1 × 10⁶ |
| koff (s⁻¹) | 6.6 × 10⁻⁸ | 5 × 10⁻⁸ | ~6 × 10⁻⁸ |
| Conditions (S-16) | pH 8.0, 25°C | pH 8.2, 22°C | — |

**(d) 95% CI for Pipeline Prediction**

Direct rate constant prediction from MD is extremely challenging. PASS criteria:
- kon: within 1–2 orders of magnitude of experimental (10⁴–10⁸ M⁻¹s⁻¹)
- koff: within 2–3 orders of magnitude (10⁻¹¹–10⁻⁵ s⁻¹)
- Alternatively, PMF barrier height provides indirect kinetic validation (see F-10)

---

### F-10: BPTI-Trypsin Association Activation Energy

**(a) Feature Description**

The activation energy for BPTI-trypsin association (Ea = 10.5 kcal/mol) reflects the energetic barrier on the association PMF — the height of the transition state separating the unbound and bound states along the binding coordinate. This barrier arises from desolvation of the interface, orientational alignment of the RSL with the active site, and conformational adjustment. The pipeline can validate this by computing the PMF along the unbinding reaction coordinate (ξ) and measuring the barrier height between the bound minimum and the transition state.

**(b) Literature Source(s)**

- **Vincent & Lazdunski (1972)** — *Biochemistry*, 11(16), 2967–2977. **Fig. 3, p. 2973.** Arrhenius plot of kon vs. 1/T; Ea derived from slope.
- **Source Report:** `vincent_lazdunski_1972_report.md` (S-16, Feature 4)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Ea (association) | 10.5 kcal/mol |
| Temperature range | 5–40°C |
| Method | Arrhenius plot: ln(kon) vs. 1/T |
| pH dependence | kon maximal at pH 8–10; controlled by His46 pK 7.05 |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 1.0 kcal/mol | Estimated from Arrhenius fit uncertainty |
| σ_comp | 2.0 kcal/mol | PMF barrier height depends critically on umbrella window spacing and sampling |
| σ_method | 2.0 kcal/mol | Transition state identification is challenging for protein-protein systems |
| **σ_combined** | **3.0 kcal/mol** | √(1.0² + 2.0² + 2.0²) |
| **95% CI** | **4.6 to 16.4 kcal/mol** | 10.5 ± (1.96 × 3.0) |

---

### F-11: KLK5 Substrate Catalytic Parameters

**(a) Feature Description**

Michaelis-Menten kinetic parameters (Km, kcat, kcat/Km) for KLK5 against six fluorogenic substrates characterize the enzyme's catalytic efficiency and subsite preferences. While Km and kcat are not directly predictable by equilibrium MD simulations, they provide qualitative validation: the pipeline should predict stable substrate-like binding in the KLK5 active site with geometries consistent with catalytic competence (proper alignment of the scissile bond with the catalytic triad). The substrate rank order by kcat/Km validates the computed binding energetics hierarchy.

**(b) Literature Source(s)**

- **Michael et al. (2005)** — *J. Biol. Chem.*, 280(15), 14628–14635. **Table II, p. 14632.** Full Michaelis-Menten kinetics for 7 substrates.
- **Source Report:** `michael_et_al_2005_report.md` (S-06, Feature 1)

**(c) Experimental Value(s)**

| Substrate | Km (mM) | kcat (min⁻¹) | kcat/Km (mM⁻¹min⁻¹) |
|-----------|---------|--------------|----------------------|
| Boc-VPR-AMC | 0.20 ± 0.01 | 196.8 | 946.5 |
| Boc-FSR-AMC | 0.19 ± 0.01 | 169.5 | 877.4 |
| Boc-QAR-AMC | 0.61 ± 0.03 | 107.0 | 175.3 |
| Boc-LKR-AMC | 1.01 ± 0.10 | 49.4 | 48.9 |
| Tos-GPR-AMC | 1.69 ± 0.32 | 20.8 | 12.2 |
| Boc-VLK-AMC | 0.64 ± 0.17 | 4.6 | 7.1 |

KLK5 optimal pH = 8.0 (S-06).

**(d) 95% CI for Pipeline Prediction**

Not directly benchmarkable by equilibrium MD. Qualitative validation:
- **PASS:** Pipeline predicts substrate binding in correct orientation (P1 Arg in S1 pocket); catalytic geometry maintained (His57–Asp102–Ser195 triad distances < 3.5 Å).
- **MARGINAL:** Correct general binding mode but catalytic triad geometry distorted.
- **FAIL:** Substrate does not bind in the active site, or P1 residue does not occupy S1 pocket.

---

### F-12: LEKTI-KLK5 pH-Dependent Binding Kinetics

**(a) Feature Description**

The LEKTI-KLK5 interaction exhibits dramatic pH dependence: association rates (ka) and dissociation rates (kd) vary 10–100× between pH 7.5 (tight binding, physiological) and pH 4.5 (weak binding, stratum corneum surface). This pH switch provides the molecular mechanism for desquamation regulation — LEKTI releases KLK5 at the acidic skin surface, enabling controlled proteolysis. The pipeline can test pH effects by simulating at different protonation states of titratable residues (His, Glu, Asp).

**(b) Literature Source(s)**

- **Deraison et al. (2007)** — *Mol. Biol. Cell*, 18, 3607–3619. **Fig. 2, Fig. 3, Table 1.** SPR sensorgrams at pH 4.5, 5.5, 6.5, and 7.5; full ka/kd for each fragment at each pH.
- **Source Report:** `deraison_et_al_2007_report.md` (S-09, Feature 3)

**(c) Experimental Value(s)**

| LEKTI Fragment | pH | ka (M⁻¹s⁻¹) | kd (s⁻¹) | KD (M) |
|----------------|-----|-------------|----------|--------|
| D8–D11 | 7.5 | High | Very low | 1.1 × 10⁻¹² |
| D8–D11 | 4.5 | Reduced ~10–100× | Increased ~10–100× | ~10⁻⁸–10⁻⁹ |

Detailed pH-resolved rate constants for all four fragments (D5, D6, D8-D11, D9-D15) are tabulated in the source report.

**(d) 95% CI for Pipeline Prediction**

Semi-quantitative: the pipeline should predict that ΔG_bind is weaker at pH 4.5 vs. pH 7.5 by at least 3–5 kcal/mol for the relevant LEKTI fragments. PASS: ΔΔG(pH4.5 − pH7.5) > 2 kcal/mol (weaker binding at acidic pH). FAIL: No pH effect or reversed direction.

---

### F-13: Barnase-Barstar Association and Dissociation Kinetics

**(a) Feature Description**

The barnase-barstar PPI is the fastest-associating protein pair characterized at the time (kon = 3.7 × 10⁸ M⁻¹s⁻¹, approaching the diffusion limit). This reflects strong long-range electrostatic steering between the highly complementary charged surfaces. As a PPI kinetics reference, it provides a calibration point for MSM-derived rate constants and validates whether the pipeline can distinguish between electrostatically steered (barnase-barstar, kon ~ 10⁸) and conformationally gated (BPTI-trypsin, kon ~ 10⁶) association mechanisms.

**(b) Literature Source(s)**

- **Schreiber & Fersht (1995)** — *J. Mol. Biol.*, 248, 478–486. **Table 1, p. 480.** Wild-type kinetics and DMC methodology.
- **Source Report:** `schreiber_fersht_1995_report.md` (S-22, Features 1 and 2)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Ka | 10¹⁴ M⁻¹ (Kd = 0.01 pM) |
| ΔG_bind (298 K) | −19.0 kcal/mol |
| kon | 3.7 × 10⁸ M⁻¹s⁻¹ |
| koff | 3.7 × 10⁻⁶ s⁻¹ |
| Conditions | 50 mM Tris-HCl pH 8.0, 25°C |

**(d) 95% CI for Pipeline Prediction**

Primarily a reference system. If tested:
- kon within 1 order of magnitude (10⁷–10⁹ M⁻¹s⁻¹): PASS
- ΔG_bind within ±5 kcal/mol of −19.0: PASS (given the extreme tightness)

---

## 5. Category III — Structural Features

### F-14: Canonical Binding Loop Backbone Geometry

**(a) Feature Description**

All canonical serine protease inhibitors maintain a highly conserved backbone conformation at the P3–P1 residues of the reactive site loop. The backbone dihedral angles (φ, ψ) at these positions are constrained to narrow ranges that ensure the scissile bond is presented to the protease active site in the correct orientation for catalytic engagement. During MD simulation, the binding loop should maintain these backbone angles in the equilibrated complex. Deviation from the canonical ranges indicates either force field inadequacy or insufficient equilibration.

**(b) Literature Source(s)**

- **Bode & Huber (1992)** — *Eur. J. Biochem.*, 204, 433–451. **Table 1, Section 4.** Compiled from multiple crystal structures.
- **Krowarsch et al. (2003)** — *Cell. Mol. Life Sci.*, 60, 2427–2444. **Section 3, Fig. 2.** Confirms loop geometry conservation.
- **Hecht et al. (1991)** — *J. Mol. Biol.*, 220, 711–722. **Section 3.3.** PSTI binding loop rmsd 0.25–0.35 Å.
- **Source Reports:** `bode_huber_1992_report.md` (S-13, Feature 1), `krowarsch_et_al_2003_report.md` (S-14), `hecht_et_al_1991_report.md` (S-19)

**(c) Experimental Value(s)**

| Position | φ Range (°) | ψ Range (°) |
|----------|-------------|-------------|
| P3 | −140 to −120 | 140 to 170 |
| P2 | −100 to −60 | 139 to 180 |
| P1 | −120 to −95 | 9 to 50 |

**(d) 95% CI for Pipeline Prediction**

Allow ±15° beyond the crystallographic ranges for MD thermal fluctuations:
- **PASS:** Mean φ/ψ values during production trajectory fall within the expanded ranges (±15° from crystal bounds).
- **MARGINAL:** Mean values deviate by 15–30° from crystal bounds.
- **FAIL:** Mean values deviate by >30° or the binding loop unfolds.

---

### F-15: Catalytic Contact Distance (P1 Carbonyl to Serine Oγ)

**(a) Feature Description**

In all canonical protease-inhibitor complexes, the carbonyl carbon of the P1 residue sits at a sub-van der Waals distance from the Oγ of the catalytic serine (Ser195, chymotrypsin numbering). This close contact (~2.7 Å) reflects the tetrahedral intermediate geometry of the trapped acyl-enzyme intermediate. Maintaining this distance in MD simulation is essential for validating that the complex is structurally sound and the inhibitory mechanism is correctly represented.

**(b) Literature Source(s)**

- **Laskowski & Kato (1980)** — *Annu. Rev. Biochem.*, 49, 593–626. **Section on Standard Mechanism.** Reports catalytic contact at 2.6 Å.
- **Bode & Huber (1992)** — *Eur. J. Biochem.*, 204, 433–451. **Section 4, p. 440.** Reports 2.7 Å.
- **Source Reports:** `laskowski_kato_1980_report.md` (S-12, Feature 2), `bode_huber_1992_report.md` (S-13, Feature 2)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| P1 C=O carbon to Ser195 Oγ distance | 2.6–2.7 Å |
| Covalent bond C–O | 1.43 Å (for reference) |
| van der Waals contact C···O | ~3.2 Å |
| Significance | Sub-vdW distance indicates incipient covalent bond (tetrahedral intermediate) |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate |
|-----------|----------|
| Crystal value | 2.65 ± 0.05 Å (consensus of S-12 and S-13) |
| MD thermal fluctuation | ±0.3 Å |
| Force field accuracy | ±0.2 Å (non-covalent potential cannot fully capture covalent intermediate) |
| **95% CI** | **2.0 to 3.3 Å** |

**Note:** Standard non-reactive force fields (AMBER ff14SB) represent the complex as a non-covalent encounter complex, not the true tetrahedral intermediate. The distance may be slightly larger (3.0–3.5 Å) in MD compared to the crystal structure, which captures the covalent intermediate. A distance of 3.0–3.5 Å in MD should be classified as MARGINAL, not FAIL, due to this known force field limitation.

---

### F-16: Protease-Inhibitor Buried Surface Area

**(a) Feature Description**

The buried surface area (BSA) at protein-protein interfaces is a key structural descriptor of binding complementarity. For protease-inhibitor complexes, BSA averages 1530 ± 170 Å² (total, both sides combined), which is within the "standard-size" PPI interface range (1600 ± 400 Å²). Specific protease-inhibitor complexes with deposited crystal structures provide individual BSA values that serve as direct validation targets for MD-computed SASA differences. The pipeline computes BSA as SASA(protease) + SASA(inhibitor) − SASA(complex).

**(b) Literature Source(s)**

- **Lo Conte et al. (1999)** — *J. Mol. Biol.*, 285, 2177–2198. **Table 3, Table 5.** BSA for 75 PPI complexes including 6 protease-inhibitor entries.
- **Krowarsch et al. (2003)** — *Cell. Mol. Life Sci.*, 60, 2427–2444. **Section 4.** Reports 1400 Å² for canonical protease-inhibitor.
- **Bode & Huber (1992)** — *Eur. J. Biochem.*, 204, 433–451. **Section 5.** Per-chain contact surface 600–900 Å².
- **Source Reports:** `lo_conte_et_al_1999_report.md` (S-24, Feature 1), `krowarsch_et_al_2003_report.md` (S-14, Feature 2), `bode_huber_1992_report.md` (S-13, Feature 3)

**(c) Experimental Value(s)**

| System | BSA (Å²) | PDB |
|--------|---------|-----|
| **Protease-inhibitor mean** | **1530 ± 170** | — |
| Trypsinogen-PSTI (1tgs) | 1730 | 1TGS |
| Kallikrein A-PTI (2kai) | 1440 | 2KAI |
| Trypsin-PTI (2ptc) | 1570 | 2PTC |
| Chymotrypsin-OMTKY3 (1cho) | 1400 | 1CHO |
| Elastase-OMTKY3 (1ppf) | 1500 | 1PPF |
| α-Chymotrypsin-eglin c (1acb) | 1460 | 1ACB |

**(d) 95% CI for Pipeline Prediction**

For the SPINK7-KLK5 complex (no co-crystal available; homology-based model):

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 170 Å² | Standard deviation across 6 protease-inhibitor complexes |
| σ_comp | 100 Å² | SASA fluctuations during equilibrium MD |
| σ_method | 150 Å² | Homology model uncertainty; no co-crystal |
| **σ_combined** | **250 Å²** | √(170² + 100² + 150²) |
| **95% CI** | **1040 to 2020 Å²** | 1530 ± (1.96 × 250) |

---

### F-17: Interface Hydrogen Bond Count

**(a) Feature Description**

Protease-inhibitor interfaces are stabilized by networks of hydrogen bonds between backbone and side-chain atoms. The mean count of ~10 H-bonds for protease-inhibitor complexes (from analysis of 75 PPI crystal structures) provides a quantitative target for the pipeline's H-bond analysis module. Both backbone-backbone (RSL to protease main chain) and side chain-side chain H-bonds contribute. The H-bond count is correlated with, but not linearly proportional to, the binding affinity.

**(b) Literature Source(s)**

- **Lo Conte et al. (1999)** — *J. Mol. Biol.*, 285, 2177–2198. **Table 5.** Mean H-bond count per interface class.
- **Krowarsch et al. (2003)** — *Cell. Mol. Life Sci.*, 60, 2427–2444. **Section 4.** Reports 8–15 H-bonds.
- **Source Reports:** `lo_conte_et_al_1999_report.md` (S-24, Feature 2), `krowarsch_et_al_2003_report.md` (S-14, Feature 3)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Mean H-bond count (protease-inhibitor) | ~10 |
| Range across complexes | 8–15 |
| H-bond criteria | Standard geometric criteria (D-A distance < 3.5 Å, D-H-A angle > 120°) |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate |
|-----------|----------|
| σ_exp | 3 (from range 8–15) |
| σ_comp | 2 (from H-bond fluctuations in MD trajectory — time-averaged) |
| **95% CI** | **3 to 17 H-bonds** | 10 ± (1.96 × 3.6) |

---

### F-18: Interface Water Molecule Count

**(a) Feature Description**

Approximately 15 water molecules are found at the protease-inhibitor interface in crystal structures, mediating hydrogen bonds between the protein surfaces and contributing to binding specificity. Explicit-solvent MD simulations should reproduce this interfacial hydration pattern. Waters at the interface can be either buried (fully enclosed by protein atoms) or partially exposed. The pipeline can count interface waters by identifying water molecules within a cutoff distance of both binding partners.

**(b) Literature Source(s)**

- **Lo Conte et al. (1999)** — *J. Mol. Biol.*, 285, 2177–2198. **Table 5, Section 3.** Interface water statistics for protease-inhibitor complexes.
- **Source Report:** `lo_conte_et_al_1999_report.md` (S-24, Feature 3)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Mean interface waters (protease-inhibitor) | ~15 |
| Definition | Water within H-bond distance of both chains |

**(d) 95% CI for Pipeline Prediction**

Interface water count fluctuates significantly in MD due to water exchange dynamics.
- **PASS:** Time-averaged interface water count 8–25.
- **MARGINAL:** 4–8 or 25–35.
- **FAIL:** <4 (interface dehydration) or >35 (interface disruption).

---

### F-19: Interface Non-Polar Fraction and Packing Quality

**(a) Feature Description**

Protease-inhibitor interfaces are 61% non-polar (by atom composition), similar to protein interiors, reflecting the hydrophobic driving force for complex formation. The packing quality, measured by the ratio of observed to expected van der Waals volume (V/Vo = 1.00 ± 0.01), indicates that protease-inhibitor interfaces achieve perfect close-packing comparable to the protein interior — there are no significant cavities or voids at the interface. This remarkable complementarity is a hallmark of evolved protease-inhibitor interactions and distinguishes them from crystal packing contacts.

**(b) Literature Source(s)**

- **Lo Conte et al. (1999)** — *J. Mol. Biol.*, 285, 2177–2198. **Table 6, Section 3.4.** Atom composition and packing analysis.
- **Source Report:** `lo_conte_et_al_1999_report.md` (S-24, Features 4 and 5)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Non-polar atom fraction at interface | 61% |
| Packing quality (V/Vo) | 1.00 ± 0.01 |
| Comparison: protein interior V/Vo | 1.00 |
| Comparison: crystal packing V/Vo | 1.02–1.05 |

**(d) 95% CI for Pipeline Prediction**

- Non-polar fraction: **PASS** if 50–75%; **FAIL** if <40% or >80%.
- Packing quality V/Vo: **PASS** if 0.95–1.05; **FAIL** if <0.90 or >1.10.

---

### F-20: Binding Loop Structural Conservation Across Kazal Family

**(a) Feature Description**

The reactive site loop of Kazal-type inhibitors adopts a remarkably conserved backbone conformation across all characterized family members, with pairwise Cα rmsd values of 0.25–0.35 Å for the P3–P3' hexapeptide segment. This structural conservation despite only ~30–40% sequence identity reflects the stringent geometric requirements for productive engagement with the protease active site. During MD, the SPINK7 binding loop should maintain this canonical conformation; significant deviation would suggest force field inadequacy or incorrect initial model geometry.

**(b) Literature Source(s)**

- **Hecht et al. (1991)** — *J. Mol. Biol.*, 220, 711–722. **Section 3.3, p. 718.** PSTI binding loop comparison to ovomucoid, Ascaris inhibitor, and domain 3 of Japanese quail ovomucoid.
- **Bode & Huber (1992)** — *Eur. J. Biochem.*, 204, 433–451. **Section 4.** General conservation statement.
- **Source Reports:** `hecht_et_al_1991_report.md` (S-19, Feature 3), `bode_huber_1992_report.md` (S-13)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Binding loop (P3–P3') Cα rmsd across Kazal family | 0.25–0.35 Å |
| Compared structures | PSTI, OMTKY3, OMSVP3, Ascaris inhibitor |

**(d) 95% CI for Pipeline Prediction**

SPINK7 RSL rmsd to cognate Kazal structures during MD production:
- **PASS:** Mean RSL Cα rmsd < 1.0 Å from starting conformation.
- **MARGINAL:** 1.0–2.0 Å.
- **FAIL:** >2.0 Å (loop distortion or unfolding).

---

### F-21: KLK5 Subsite Specificity Profile

**(a) Feature Description**

Positional scanning combinatorial library (PS-SCL) profiling of KLK5 reveals the most stringent P1 Arg selectivity among all characterized kallikreins. The extended subsite preferences (P1: exclusively Arg; P2: Ser slightly preferred; P3: Met, Phe, Tyr accepted; P4: unusual Gly preference) define the molecular recognition landscape that any SPINK7-like inhibitor must satisfy. The pipeline can validate this by computing per-subsite interaction energies between the inhibitor RSL and the KLK5 binding cleft.

**(b) Literature Source(s)**

- **Debela et al. (2006)** — *J. Biol. Chem.*, 281(35), 25678–25688. **Table I, Table II.** PS-SCL profiles for 7 kallikreins; KLK5 active-site residue identities.
- **Source Report:** `debela_et_al_2006_report.md` (S-07, Features 1 and 2)

**(c) Experimental Value(s)**

| Subsite | KLK5 Preference | Selectivity Strength |
|---------|-----------------|---------------------|
| P1 | Arg (exclusive) | Strongest among all 7 KLKs profiled |
| P2 | Ser > Asn, Thr, Ala; Arg and Tyr rejected | Moderate |
| P3 | Met, Phe, Tyr; Asp rejected | Moderate |
| P4 | Gly (unusual), then Tyr, Val, Pro | Moderate |

The S1 pocket of KLK5 contains Asp189 (chymotrypsin numbering), providing the electrostatic anchor for P1-Arg.

**(d) 95% CI for Pipeline Prediction**

Qualitative structural validation:
- **PASS:** P1 Arg of SPINK7 occupies the S1 pocket with stable salt bridge to Asp189; P2–P4 residues occupy the respective subsites.
- **FAIL:** P1 is not in the S1 pocket, or salt bridge to Asp189 is absent.

---

### F-22: SPINK7-KLK5 Complex Binding Geometry

**(a) Feature Description**

ClusPro docking of SPINK7 (PDB 2LEO) to KLK5 (PDB 2PSX) predicts the canonical binding mode in which the SPINK7 inhibitory loop is positioned ~2.5 Å from the catalytic serine residue in the KLK5 active site. Two additional docking algorithms confirmed the same binding mode. This docking-derived geometry serves as the starting point for all MD simulations and the initial complex structure must be validated against this prediction during equilibration.

**(b) Literature Source(s)**

- **Azouz et al. (2020)** — *Sci. Transl. Med.*, 12(545), eaaz7773. **Fig. 1C, fig. S1D.** ClusPro docking: rotational screening (~10⁹ conformations) → 1000 top scored → clustered at Cα rmsd < 9 Å → largest cluster selected. Three docking algorithms converge.
- **Source Report:** `azouz_et_al_2020_report.md` (S-03, Feature 4)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| RSL-to-catalytic-Ser distance | ~2.5 Å |
| PDB inputs | SPINK7: 2LEO (chain A); KLK5: 2PSX (chain A) |
| Method | ClusPro docking (confirmed by 2 additional algorithms) |
| Binding mode | Canonical: inhibitory loop inserts into active-site cleft |

**Note:** This is a computational docking prediction, not an experimental crystal structure. However, it is the only available structural model of the SPINK7-KLK5 complex and is supported by three independent docking algorithms converging on the same binding mode.

**(d) 95% CI for Pipeline Prediction**

After equilibration, the RSL-to-Ser distance should remain consistent with canonical protease-inhibitor geometry:
- **PASS:** RSL P1 carbonyl carbon within 2.0–4.0 Å of catalytic Ser Oγ (consistent with F-15).
- **MARGINAL:** 4.0–6.0 Å (partial engagement).
- **FAIL:** >6.0 Å (complex dissociation or non-canonical binding mode).

---

### F-23: BPTI High-Resolution Reference Structure

**(a) Feature Description**

The BPTI crystal structure at 0.94 Å X-ray resolution (joint with 1.8 Å neutron refinement) provides the ultimate structural validation target for MD simulations. At this resolution, hydrogen atom positions, methyl rotor orientations, and individual water molecule B-factors are determined. With 63 ordered water molecules (4 internal to the protein), this structure defines the gold-standard starting point for BPTI simulations. The neutron data also reveal the protonation states of titratable residues and detailed hydrogen bonding geometries, serving as a reference for force field validation.

**(b) Literature Source(s)**

- **Wlodawer et al. (1984)** — *J. Mol. Biol.*, 180, 301–329. **Table 3, Table 5, Table 7.** Complete refinement statistics, water positions, B-factors.
- **Source Report:** `wlodawer_et_al_1984_report.md` (S-18, Features 1 and 2)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| X-ray resolution | 0.94 Å (R = 0.200) |
| Neutron resolution | 1.8 Å (R = 0.197) |
| Space group | P2₁2₁2₁ |
| Unit cell | a = 74.1, b = 23.4, c = 28.9 Å |
| Ordered water molecules | 63 (4 internal) |
| PDB codes | 4PTI (form II), 5PTI (form I) |

**(d) 95% CI for Pipeline Prediction**

For BPTI structural stability during MD:
- **PASS:** Backbone Cα rmsd < 1.5 Å from crystal structure during production trajectory.
- **MARGINAL:** 1.5–2.5 Å.
- **FAIL:** >2.5 Å (significant structural deviation).

---

## 6. Category IV — Dynamic Features

### F-24: BPTI Backbone Amide H/D Exchange Protection Pattern

**(a) Feature Description**

Neutron diffraction of BPTI crystals soaked in D₂O for 3 months reveals that 11 backbone amide hydrogens resist exchange with deuterium, indicating that these amides are deeply buried within the protein core or involved in very stable hydrogen bonds. The protection pattern provides a site-specific map of backbone rigidity that can be compared to MD-derived N-H order parameters (S²) or solvent accessibility profiles. Protected amides should correspond to residues with high S² values (>0.85) and low backbone RMSF in MD simulations.

**(b) Literature Source(s)**

- **Wlodawer et al. (1984)** — *J. Mol. Biol.*, 180, 301–329. **Table 5, p. 316.** Identification of 11 exchanged and 11 protected amides from neutron density maps.
- **Source Report:** `wlodawer_et_al_1984_report.md` (S-18, Feature 3)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Protected amide hydrogens | 11 positions (after 3 months D₂O) |
| Exchange conditions | Crystal soaking in D₂O, 3 months |
| Implication | These amides are in the most rigid, solvent-inaccessible regions |

**(d) 95% CI for Pipeline Prediction**

- **PASS:** ≥8 of the 11 protected amides show S² > 0.85 or backbone RMSF < 0.5 Å.
- **MARGINAL:** 5–7 of 11 agree.
- **FAIL:** <5 of 11 agree with experimental protection pattern.

---

### F-25: BPTI Crystal Form Conformational Variability

**(a) Feature Description**

BPTI has been crystallized in two forms: orthorhombic form I (PDB 5PTI) and orthorhombic form II (PDB 4PTI). Comparison of these two independently determined structures reveals a main-chain rmsd of 0.40 Å, providing a measure of the intrinsic conformational variability of the protein. This variability represents the real structural "noise floor" — MD-computed conformational fluctuations should be at least this large (since crystal packing constrains conformational space) and ideally not much larger for a well-folded globular protein.

**(b) Literature Source(s)**

- **Wlodawer et al. (1984)** — *J. Mol. Biol.*, 180, 301–329. **Section 4, p. 322.** Comparison of form I and form II structures.
- **Source Report:** `wlodawer_et_al_1984_report.md` (S-18, Feature 4)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Form I/II main-chain rmsd | 0.40 Å |
| Form I/II all-atom rmsd | ~0.55–0.60 Å (estimated from main chain + side chain) |

**(d) 95% CI for Pipeline Prediction**

MD conformational ensemble width (mean pairwise Cα rmsd):
- **PASS:** 0.3–1.5 Å (captures crystal variability without excessive fluctuation).
- **MARGINAL:** 1.5–2.5 Å.
- **FAIL:** <0.2 Å (frozen simulation) or >2.5 Å (unfolding).

---

## 7. Category V — Biophysical Features

### F-26: BSA-ΔG Empirical Correlation for PPIs

**(a) Feature Description**

An empirical relationship exists between the buried surface area (BSA) at protein-protein interfaces and the binding free energy: ΔG_bond ≈ 0.025 × BSA (kcal/mol per Å²). This correlation, though weak (r = 0.54 for 70 rigid-body complexes), provides a rapid sanity check for MD results: if the pipeline computes a BSA of ~1530 Å² for SPINK7-KLK5, the expected ΔG from this relationship would be ~38 kcal/mol, which is far larger than the experimentally observed −9.4 kcal/mol. This discrepancy (common for all PPIs) reflects the fact that desolvation entropy penalties offset much of the favorable contact energy, and only a fraction of the buried surface contributes net binding energy. Nevertheless, the BSA-ΔG relationship serves as a coarse consistency check.

**(b) Literature Source(s)**

- **Kastritis & Bonvin (2013)** — *J. R. Soc. Interface*, 10, 20120835. **Section 3, Fig. 4.** BSA-ΔG regression for rigid-body complexes.
- **Source Report:** `kastritis_bonvin_2013_report.md` (S-23, Feature 1)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Empirical relationship | ΔG_bond ≈ 0.025 × BSA (kcal/mol per Å²) |
| Correlation coefficient | r = 0.54 (70 rigid-body complexes) |
| Experimental Kd error | 20–50% (→ 0.1–0.25 kcal/mol ΔG error) |

**(d) 95% CI for Pipeline Prediction**

This is a validation relationship, not a direct prediction target. The pipeline should:
- Compute BSA and ΔG independently.
- Confirm that the BSA/ΔG ratio is consistent with other protease-inhibitor systems (not necessarily the 0.025 regression line, which has large scatter).

---

### F-27: Protease-Inhibitor Interface Close-Packing Quality

**(a) Feature Description**

The packing quality at protease-inhibitor interfaces, measured as the ratio of observed to expected van der Waals volume (V/Vo), equals 1.00 ± 0.01 — identical to the protein interior. This perfect packing reflects the exquisite shape complementarity evolved at these interfaces and distinguishes biological recognition interfaces from crystal packing contacts (V/Vo = 1.02–1.05). The pipeline can compute Voronoi-based packing metrics from equilibrated MD snapshots.

**(b) Literature Source(s)**

- **Lo Conte et al. (1999)** — *J. Mol. Biol.*, 285, 2177–2198. **Table 6.** Voronoi packing analysis.
- **Source Report:** `lo_conte_et_al_1999_report.md` (S-24, Feature 5)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| V/Vo (protease-inhibitor) | 1.00 ± 0.01 |
| V/Vo (protein interior) | 1.00 |
| V/Vo (crystal packing) | 1.02–1.05 |

**(d) 95% CI for Pipeline Prediction**

- **PASS:** V/Vo = 0.95–1.05.
- **MARGINAL:** V/Vo = 1.05–1.10 (slightly worse packing than crystal).
- **FAIL:** V/Vo > 1.10 (cavity formation at interface) or V/Vo < 0.90 (steric clashes).

---

### F-28: Inhibitor Scaffold Energetic Contribution

**(a) Feature Description**

The non-contact scaffold of canonical protease inhibitors (disulfide bonds, hydrophobic core, secondary structure elements distant from the binding loop) contributes approximately −33 kJ/mol (−7.9 kcal/mol) to the binding free energy. This scaffold energy reflects the conformational pre-organization of the binding loop: by maintaining the RSL in its catalytically competent conformation, the scaffold eliminates the conformational entropy penalty that a flexible peptide would pay upon binding. The pipeline can estimate this contribution by comparing the binding free energy of the intact inhibitor with that of the isolated binding loop peptide.

**(b) Literature Source(s)**

- **Bode & Huber (1992)** — *Eur. J. Biochem.*, 204, 433–451. **Section 5, p. 448.** Scaffold energy estimate from comparison of intact inhibitors vs. reactive-site loop peptides.
- **Source Report:** `bode_huber_1992_report.md` (S-13, Feature 4)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Scaffold energy contribution | −33 kJ/mol (−7.9 kcal/mol) |
| Thermal stability (BPTI) | Tm = 95°C |
| Thermal stability (OMTKY3) | Tm = 85°C |
| Evidence | ΔΔG(intact inhibitor − RSL peptide) from binding studies |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 2.0 kcal/mol | Estimated from limited data on loop peptide binding |
| σ_comp | 2.0 kcal/mol | Requires two separate free energy calculations |
| σ_method | 3.0 kcal/mol | Each calculation carries ±2 kcal/mol systematic error |
| **σ_combined** | **4.1 kcal/mol** | √(2.0² + 2.0² + 3.0²) |
| **95% CI** | **−15.9 to 0.1 kcal/mol** | −7.9 ± (1.96 × 4.1) |

---

### F-29: SH3-p41 Binding Free Energy (Computational Method Validation)

**(a) Feature Description**

The SH3 domain–p41 peptide complex (PDB 1BBZ) serves as the definitive computational methods validation target for binding free energy calculations. Three independent approaches — PMF via REMD-US (ΔG = −7.8 kcal/mol), alchemical double-decoupling (ΔG = −7.7 kcal/mol), and experimental measurement (ΔG = −7.99 kcal/mol) — converge within ~0.2 kcal/mol, establishing that rigorous free energy methods can achieve sub-kcal/mol accuracy for well-characterized protein-peptide systems. This convergence validates the PMF + REMD-US approach as the recommended method for our pipeline, while also demonstrating that MM/PBSA (ΔG = −2.6 kcal/mol) fails catastrophically for this system.

**(b) Literature Source(s)**

- **Gumbart et al. (2013)** — *J. Chem. Theory Comput.*, 9(1), 794–802. **Table 1, Figs. 2–5.** Full thermodynamic decomposition and convergence analysis.
- **Deng & Roux (2009)** — *J. Phys. Chem. B*, 113, 2234–2246. **Section 3.** DDM and PMF methodology framework.
- **Source Reports:** `gumbart_et_al_2013_report.md` (S-28, Features 1–5), `deng_roux_2009_report.md` (S-25)

**(c) Experimental Value(s)**

| Method | ΔG_bind (kcal/mol) | Uncertainty |
|--------|-------------------|-------------|
| **Experimental** | **−7.99** | from Pisabarro & Serrano 1996 |
| PMF (REMD-US) | −7.8 | ±0.9 |
| Alchemical (DDM) | −7.7 | ±0.6 |
| MM/PBSA | −2.6 | N/A (systematic failure) |

**REMD-US Protocol (S-28):** ~115 ns total; 6 restraints (3 conformational, 3 orientational); CHARMM22 with CMAP; separation along COM-COM axis.

**(d) 95% CI for Pipeline Prediction**

If the pipeline is applied to SH3-p41 as a positive control:

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.3 kcal/mol | Well-characterized Kd |
| σ_comp | 0.9 kcal/mol | From S-28 bootstrap uncertainty |
| σ_method | 1.0 kcal/mol | Force field (AMBER vs CHARMM) and protocol differences |
| **σ_combined** | **1.4 kcal/mol** | √(0.3² + 0.9² + 1.0²) |
| **95% CI** | **−10.7 to −5.3 kcal/mol** | −7.99 ± (1.96 × 1.4) |

---

## 8. Category VI — Mutational Features

### F-30: BPTI Reactive Region Alanine Scan ΔΔG Panel

**(a) Feature Description**

The most comprehensive alanine scanning mutagenesis study of any protease-inhibitor system: 15 point mutations spanning the BPTI reactive region (T11A through G36A) with full kinetic (kon, koff, Ki) and thermodynamic (ΔG, ΔH, TΔS by van't Hoff analysis) characterization against both trypsin and chymotrypsin. The data reveal that (1) the P1 residue K15 is overwhelmingly dominant (K15A: ΔΔG ≈ 10 kcal/mol, Ki increases from 5 × 10⁻¹⁴ to 1.4 × 10⁻⁶ M), (2) most mutants affect koff more than kon (binding affinity is determined by complex stability, not association speed), and (3) enthalpy-entropy compensation operates across mutants. This dataset is the ideal validation target for alchemical FEP calculations of ΔΔG upon point mutation.

**(b) Literature Source(s)**

- **Castro & Anderson (1996)** — *Biochemistry*, 35(35), 11435–11446. **Tables 1, 2, and 3.** Complete dataset for all 15 mutants.
- **Source Report:** `castro_anderson_1996_report.md` (S-17, Features 2–8)

**(c) Experimental Value(s)**

Key mutants (trypsin, pH 8.2, 22°C):

| Mutant | Ki (M) | ΔΔG (kcal/mol) | Primary Effect |
|--------|--------|----------------|----------------|
| Wild-type | 5 × 10⁻¹⁴ | 0 (reference) | — |
| K15A (P1) | 1.4 × 10⁻⁶ | ≈ +10 | kon ↓ 200×; koff ↑ |
| C14A | ~10⁻⁹ | ≈ +7 | Disulfide loss |
| R17A (P2') | ~10⁻¹⁰ | ≈ +5 | Interface contact |
| I18A, I19A | ~10⁻¹¹ | ≈ +3–4 | Hydrophobic packing |
| C38A | ~10⁻⁹ | ≈ +7 | Disulfide loss |
| T11A, G12A, P13A | ~10⁻¹³ | ≈ +1–2 | Peripheral |
| G36A | ~10⁻¹⁴ | ≈ 0 | No effect |

The dataset includes full temperature-dependent ΔH and TΔS for each mutant, enabling Van't Hoff validation.

**(d) 95% CI for Pipeline Prediction**

For FEP-computed ΔΔG values:

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.3 kcal/mol | Precise Ki measurements with replicates |
| σ_comp | 0.5 kcal/mol | FEP statistical uncertainty from lambda windows |
| σ_method | 1.0 kcal/mol | Alchemical FEP for point mutations (§25.4: ±0.5–1.5) |
| **σ_combined** | **1.2 kcal/mol** | √(0.3² + 0.5² + 1.0²) |
| **95% CI per mutant** | **ΔΔG_expt ± 2.3 kcal/mol** | ΔΔG ± (1.96 × 1.2) |

For the panel as a whole: Spearman rank correlation between computed and experimental ΔΔG should be ρ > 0.7 (PASS) or ρ > 0.5 (MARGINAL).

---

### F-31: BPTI Disulfide Bridge Ablation Effect on Binding

**(a) Feature Description**

BPTI contains three disulfide bonds: Cys5-Cys55, Cys14-Cys38, and Cys30-Cys51. Selective reduction of the Cys14-Cys38 disulfide (which directly stabilizes the reactive site loop) increases the Kd by approximately 10⁵-fold (from ~6 × 10⁻¹⁴ to ~10⁻⁹ M), corresponding to ΔΔG ≈ 7 kcal/mol. This dramatic effect demonstrates that the disulfide framework is essential for maintaining the binding loop in its catalytically competent conformation. The pipeline can test this by simulating with the Cys14-Cys38 bond absent and comparing ΔG_bind to the intact inhibitor.

**(b) Literature Source(s)**

- **Vincent & Lazdunski (1972)** — *Biochemistry*, 11(16), 2967–2977. **Table II, p. 2975.** Binding of reduced BPTI variants to trypsin.
- **Source Report:** `vincent_lazdunski_1972_report.md` (S-16, Feature 5)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Kd (intact BPTI-trypsin) | 6 × 10⁻¹⁴ M |
| Kd (C14-C38 reduced BPTI-trypsin) | ~10⁻⁹ M |
| ΔΔG (reduction effect) | ≈ +7 kcal/mol |
| Selectivity of effect | C14-C38 is the active-site proximal disulfide |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 1.0 kcal/mol | Approximate Kd for reduced variant |
| σ_comp | 2.0 kcal/mol | Large structural change requires extensive re-sampling |
| σ_method | 3.0 kcal/mol | Modified topology + re-equilibration introduces uncertainty |
| **σ_combined** | **3.7 kcal/mol** | √(1.0² + 2.0² + 3.0²) |
| **95% CI** | **−0.3 to +14.3 kcal/mol** | +7.0 ± (1.96 × 3.7) |

The pipeline must predict that C14-C38 reduction weakens binding (positive ΔΔG). If the computed ΔΔG is negative (tighter binding upon reduction), this constitutes a FAIL.

---

### F-32: SPINK1 N34S Mutation Functional Impact (Negative Control)

**(a) Feature Description**

The N34S polymorphism in SPINK1 is clinically significant — it is associated with chronic pancreatitis (23% of pediatric cases vs. 0.94% controls). However, functional studies demonstrate that N34S does NOT alter trypsin inhibitor activity under any tested condition (pH 5–9, ±CaCl₂, ±EDTA, dissociation kinetics). This makes N34S an ideal negative control for FEP validation: the pipeline should predict ΔΔG ≈ 0 for the N34S mutation. Any significant computed ΔΔG would indicate a force field artifact or sampling inadequacy.

**(b) Literature Source(s)**

- **Kuwata et al. (2002)** — *J. Gastroenterol.*, 37, 928–934. **Table 2, Figs. 2–4.** Comprehensive pH/salt/kinetic panel showing no functional difference.
- **Source Report:** `kuwata_et_al_2002_report.md` (S-20, Features 2 and 3)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Ki (wild-type SPINK1) | 7.6 nM |
| Ki (N34S SPINK1) | 7.6 nM (no change) |
| ΔΔG (N34S − WT) | 0.0 kcal/mol |
| pH independence | Confirmed at pH 5, 6, 7, 8, 9 |
| Calcium independence | Confirmed ±CaCl₂, ±EDTA |

**(d) 95% CI for Pipeline Prediction**

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.2 kcal/mol | Measured under multiple conditions, consistently zero |
| σ_comp | 0.5 kcal/mol | FEP statistical uncertainty |
| σ_method | 0.5 kcal/mol | FEP for conservative mutation (N→S, similar size) |
| **σ_combined** | **0.7 kcal/mol** | √(0.2² + 0.5² + 0.5²) |
| **95% CI** | **−1.4 to +1.4 kcal/mol** | 0.0 ± (1.96 × 0.7) |

**PASS:** |ΔΔG| < 1.4 kcal/mol. **FAIL:** |ΔΔG| > 2.8 kcal/mol.

---

### F-33: Barnase-Barstar Double Mutant Cycle Coupling Energies

**(a) Feature Description**

Double mutant cycle (DMC) analysis of the barnase-barstar interface provides ~45 coupling energy (ΔΔGint) measurements for pairs of interfacial residues. These coupling energies quantify the cooperative contribution of specific residue-residue interactions to binding: ΔΔGint = ΔΔG(AB) − ΔΔG(A) − ΔΔG(B) + ΔΔG(WT). The strongest coupling energies (up to 7.0 kcal/mol for Asp39b*–Arg87bn) identify electrostatic hot spots. Coupling energy decreases with inter-residue distance, vanishing beyond ~7 Å. This dataset provides the most stringent test of whether FEP calculations can capture cooperative effects at PPI interfaces.

**(b) Literature Source(s)**

- **Schreiber & Fersht (1995)** — *J. Mol. Biol.*, 248, 478–486. **Tables 2 and 3, Fig. 1.** Complete DMC coupling energy dataset.
- **Source Report:** `schreiber_fersht_1995_report.md` (S-22, Features 3–5)

**(c) Experimental Value(s)**

| Parameter | Value |
|-----------|-------|
| Number of DMC pairs analyzed | ~45 |
| Maximum coupling energy | 7.0 kcal/mol (Asp39b*–Arg87bn) |
| Coupling energy range | 0.3–7.0 kcal/mol |
| Distance dependence | Coupling vanishes beyond ~7 Å inter-Cα distance |
| Wild-type Ka | 10¹⁴ M⁻¹ |

**(d) 95% CI for Pipeline Prediction**

For each DMC coupling energy:

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.3 kcal/mol | Estimated from replicate measurements |
| σ_comp | 0.5 kcal/mol | FEP for each mutant |
| σ_method | 1.5 kcal/mol | DMC requires 4 independent FEP calculations, errors compound |
| **σ_combined** | **1.6 kcal/mol** | √(0.3² + 0.5² + 1.5²) |
| **95% CI per pair** | **ΔΔGint_expt ± 3.2 kcal/mol** | ΔΔGint ± (1.96 × 1.6) |

Panel-level: Spearman rank correlation between computed and experimental coupling energies should be ρ > 0.6 (PASS) or ρ > 0.4 (MARGINAL).

---

## 9. Summary Statistics

### 9.1 Feature Counts by Category

| Category | Count | Feature IDs |
|----------|-------|-------------|
| Thermodynamic | 8 | F-01 through F-08 |
| Kinetic | 5 | F-09 through F-13 |
| Structural | 10 | F-14 through F-23 |
| Dynamic | 2 | F-24, F-25 |
| Biophysical | 4 | F-26 through F-29 |
| Mutational | 4 | F-30 through F-33 |
| **Total** | **33** | **F-01 through F-33** |

### 9.2 Feature Counts by System

| System | Count | Feature IDs |
|--------|-------|-------------|
| SPINK7-KLK5 | 3 | F-01, F-02, F-22 |
| LEKTI-KLK5 | 2 | F-03, F-12 |
| BPTI-trypsin | 8 | F-04, F-09, F-10, F-23, F-24, F-25, F-30, F-31 |
| PSTI/SPINK1-chymotrypsinogen | 1 | F-05 |
| SPINK1-trypsin | 2 | F-06, F-32 |
| KLK5 (enzyme alone) | 2 | F-11, F-21 |
| Canonical protease-inhibitor (general) | 8 | F-07, F-08, F-14, F-15, F-16, F-17, F-18, F-19 |
| Kazal family (general) | 1 | F-20 |
| Barnase-barstar | 2 | F-13, F-33 |
| SH3-p41 | 1 | F-29 |
| PPI dataset (general) | 1 | F-26 |

### 9.3 Feature Counts by Source

| Source ID | Source | Features Contributed | Feature IDs |
|-----------|--------|---------------------|-------------|
| S-01 | Rothenberg 2009 | 0 | — |
| S-02 | Azouz et al. 2018 | 0 | — |
| S-03 | Azouz et al. 2020 | 3 | F-01, F-02, F-22 |
| S-04 | Azouz & Rothenberg 2019 | 0 | — |
| S-05 | Brattsand & Egelrud 1999 | 0 | — |
| S-06 | Michael et al. 2005 | 1 | F-11 |
| S-07 | Debela et al. 2006 | 1 | F-21 |
| S-08 | Borgono et al. 2007 | 1 | F-03 |
| S-09 | Deraison et al. 2007 | 2 | F-03, F-12 |
| S-10 | Chavanas et al. 2000 | 0 | — |
| S-11 | Egelrud et al. 2005 | 1 | F-03 |
| S-12 | Laskowski & Kato 1980 | 1 | F-15 |
| S-13 | Bode & Huber 1992 | 4 | F-14, F-15, F-16, F-28 |
| S-14 | Krowarsch et al. 2003 | 4 | F-07, F-08, F-16, F-17 |
| S-15 | Rawlings et al. 2004 | 0 | — |
| S-16 | Vincent & Lazdunski 1972 | 4 | F-04, F-09, F-10, F-31 |
| S-17 | Castro & Anderson 1996 | 3 | F-04, F-09, F-30 |
| S-18 | Wlodawer et al. 1984 | 3 | F-23, F-24, F-25 |
| S-19 | Hecht et al. 1991 | 2 | F-05, F-20 |
| S-20 | Kuwata et al. 2002 | 2 | F-06, F-32 |
| S-21 | Witt et al. 2000 | 0 | — |
| S-22 | Schreiber & Fersht 1995 | 2 | F-13, F-33 |
| S-23 | Kastritis & Bonvin 2013 | 1 | F-26 |
| S-24 | Lo Conte et al. 1999 | 4 | F-16, F-17, F-18, F-19 |
| S-25 | Deng & Roux 2009 | 1 | F-29 |
| S-26 | Mobley & Gilson 2017 | 0 | — |
| S-27 | Park et al. 2003 | 0 | — |
| S-28 | Gumbart et al. 2013 | 1 | F-29 |

**Sources contributing benchmarkable features:** 20 of 28 (71%)  
**Sources providing context/methodology only:** 8 of 28 (S-01, S-02, S-04, S-05, S-10, S-15, S-21, S-26, S-27)

### 9.4 Pipeline Method Coverage

| Pipeline Method | Features Addressable | Feature IDs |
|-----------------|---------------------|-------------|
| SMD + Jarzynski | 6 | F-01, F-02, F-04, F-05, F-06, F-29 |
| Umbrella Sampling + WHAM | 7 | F-01, F-02, F-03, F-04, F-05, F-06, F-29 |
| MBAR | 7 | F-01, F-02, F-03, F-04, F-05, F-06, F-29 |
| Alchemical FEP | 5 | F-30, F-31, F-32, F-33, F-07 |
| Production MD + structural analysis | 12 | F-14, F-15, F-16, F-17, F-18, F-19, F-20, F-21, F-22, F-23, F-24, F-25 |
| MSM/TICA | 2 | F-09, F-13 |
| Energy decomposition | 3 | F-07, F-08, F-28 |
| Multi-protonation-state simulations | 1 | F-12 |
| Voronoi analysis | 2 | F-19, F-27 |

---

## 10. Coverage Assessment

### 10.1 Category Coverage

| Category | §13 Requirement | Achieved | Assessment |
|----------|----------------|----------|------------|
| Thermodynamic | ΔG, ΔH, ΔS, Kd, Ki | 8 features covering Ki, Kd, ΔG, ΔΔG, energy decomposition | **Excellent** |
| Kinetic | kon, koff, kcat, Km | 5 features covering kon, koff, Ea, Km, kcat, pH-dependent kinetics | **Good** |
| Structural | BSA, contacts, H-bonds, RMSD, Rg | 10 features covering BSA, H-bonds, waters, packing, loop geometry, docking, high-resolution reference | **Excellent** |
| Dynamic | B-factors, exchange timescales, S² | 2 features covering H/D exchange and conformational variability | **Adequate** |
| Biophysical | pH effects, salt effects, Tm | 4 features covering BSA-ΔG correlation, packing, scaffold energy, method validation | **Good** |
| Mutational | ΔΔG, P1 mutation, disulfide | 4 features covering alanine scan, disulfide ablation, negative control, DMC coupling | **Excellent** |

### 10.2 Feature Gaps

| Gap | Description | Mitigation |
|-----|-------------|------------|
| SPINK7-KLK5 co-crystal structure | No deposited crystal structure of the exact complex | Use ClusPro docking from S-03 (PDB 2LEO + 2PSX) as starting model; validate against canonical protease-inhibitor geometry (F-14, F-15) |
| SPINK7-KLK5 calorimetric data | No ITC ΔH/ΔS data for the exact system | Use Ki from Morrison equation (S-03); cross-validate with LEKTI-KLK5 data (S-08, S-09) |
| SPINK7-KLK5 kon/koff | No SPR kinetics for the exact system | Benchmark kon/koff methods against BPTI-trypsin (F-09) and LEKTI-KLK5 (F-12) first |
| Salt concentration dependence | No ionic strength titration data | Simulate at 0.15 M NaCl (physiological); qualitative comparison only |
| NMR relaxation/order parameters | No NMR data for any SPINK-KLK system | Use BPTI H/D exchange (F-24) as proxy |

### 10.3 Confidence Assessment

| Confidence Tier | Features | Justification |
|----------------|----------|---------------|
| **High** — well-characterized experimental benchmark with multiple confirming sources | F-01, F-03, F-04, F-14, F-15, F-16, F-17, F-30 | Multiple labs, precise measurements, well-understood systems |
| **Medium** — single reliable experimental measurement or consensus from reviews | F-02, F-05, F-06, F-07, F-08, F-09, F-10, F-11, F-18, F-19, F-20, F-21, F-22, F-23, F-24, F-25, F-29, F-31, F-32 | Single-source data but from reputable labs with reproducible methods |
| **Low** — estimate, correlation, or methodological validation | F-12, F-13, F-26, F-27, F-28, F-33 | Indirect measurements, weak correlations, or reference systems far from target |

---

## 11. Key PDB Codes and Structural References

| PDB | Structure | Resolution | Relevance | Features |
|-----|-----------|------------|-----------|----------|
| **2LEO** | SPINK7 (NMR ensemble) | NMR | Primary inhibitor model | F-01, F-02, F-22 |
| **2PSX** | KLK5 (X-ray) | Crystal | Primary protease model | F-01, F-02, F-11, F-21, F-22 |
| 4PTI | BPTI form II (X-ray + neutron) | 0.94 Å | Gold-standard inhibitor structure | F-04, F-23, F-24, F-25 |
| 5PTI | BPTI form I (X-ray) | ~1.5 Å | Second BPTI form | F-25 |
| 2PTC | BPTI-trypsin complex | 1.9 Å | Gold-standard protease-inhibitor complex | F-04, F-09, F-16, F-30, F-31 |
| 1TGS | Trypsinogen-PSTI/SPINK1 complex | 1.8 Å | Kazal-protease reference complex | F-05, F-16, F-20 |
| 2KAI | Kallikrein A-PTI complex | 2.1 Å | Kallikrein-inhibitor reference | F-16 |
| 1CHO | Chymotrypsin-OMTKY3 complex | 1.8 Å | Kazal Kazal-protease reference | F-05, F-16 |
| 2OVO | OMSVP3 (free Kazal inhibitor) | 1.5 Å | Kazal family fold reference | F-14, F-20 |
| 1BBZ | SH3-p41 peptide complex | NMR | Computational methods benchmark | F-29 |

---

## 12. References

1. Azouz NP et al. (2020) Functional role of kallikrein 5 and proteinase-activated receptor 2 in eosinophilic esophagitis. *Sci. Transl. Med.* 12(545), eaaz7773.
2. Azouz NP et al. (2018) The antiprotease SPINK7 serves as an inhibitory checkpoint for esophageal epithelial inflammatory responses. *Sci. Transl. Med.* 10(444), eaap9736.
3. Borgono CA et al. (2007) A potential role for multiple tissue kallikrein serine proteases in epidermal desquamation. *J. Biol. Chem.* 282(6), 3640–3652.
4. Deraison C et al. (2007) LEKTI fragments specifically inhibit KLK5, KLK7, and KLK14 and control desquamation through a pH-dependent interaction. *Mol. Biol. Cell* 18, 3607–3619.
5. Egelrud T et al. (2005) hK5 and hK7, two serine proteinases abundant in human skin, are inhibited by LEKTI domain 6. *Br. J. Dermatol.* 153, 1200–1203.
6. Michael IP et al. (2005) Biochemical and enzymatic characterization of human kallikrein 5 (hK5). *J. Biol. Chem.* 280(15), 14628–14635.
7. Debela M et al. (2006) Specificity profiling of seven human tissue kallikreins reveals individual subsite preferences. *J. Biol. Chem.* 281(35), 25678–25688.
8. Laskowski M Jr, Kato I (1980) Protein inhibitors of proteinases. *Annu. Rev. Biochem.* 49, 593–626.
9. Bode W, Huber R (1992) Natural protein proteinase inhibitors and their interaction with proteinases. *Eur. J. Biochem.* 204, 433–451.
10. Krowarsch D et al. (2003) Canonical protein inhibitors of serine proteases. *Cell. Mol. Life Sci.* 60, 2427–2444.
11. Vincent JP, Lazdunski M (1972) Trypsin-pancreatic trypsin inhibitor association. *Biochemistry* 11(16), 2967–2977.
12. Castro MJM, Anderson S (1996) Alanine point-mutations in the reactive region of bovine pancreatic trypsin inhibitor. *Biochemistry* 35(35), 11435–11446.
13. Wlodawer A et al. (1984) Structure of bovine pancreatic trypsin inhibitor: results of joint neutron and X-ray refinement. *J. Mol. Biol.* 180, 301–329.
14. Hecht HJ et al. (1991) Three-dimensional structure of the complexes between bovine chymotrypsinogen A and two recombinant variants of human pancreatic secretory trypsin inhibitor. *J. Mol. Biol.* 220, 711–722.
15. Kuwata K et al. (2002) Functional analysis of recombinant pancreatic secretory trypsin inhibitor protein with amino-acid substitution. *J. Gastroenterol.* 37, 928–934.
16. Schreiber G, Fersht AR (1995) Energetics of protein-protein interactions: analysis of the barnase-barstar interface by single mutations and double mutant cycles. *J. Mol. Biol.* 248, 478–486.
17. Kastritis PL, Bonvin AMJJ (2013) On the binding affinity of macromolecular interactions: daring to ask why proteins interact. *J. R. Soc. Interface* 10, 20120835.
18. Lo Conte L, Chothia C, Janin J (1999) The atomic structure of protein-protein recognition sites. *J. Mol. Biol.* 285, 2177–2198.
19. Deng Y, Roux B (2009) Computations of standard binding free energies with molecular dynamics simulations. *J. Phys. Chem. B* 113, 2234–2246.
20. Mobley DL, Gilson MK (2017) Predicting binding free energies: frontiers and benchmarks. *Annu. Rev. Biophys.* 46, 531–558.
21. Park S et al. (2003) Free energy calculation from steered molecular dynamics simulations using Jarzynski's equality. *J. Chem. Phys.* 119(6), 3559–3566.
22. Gumbart JC, Roux B, Chipot C (2013) Standard binding free energies from computer simulations: what is the best strategy? *J. Chem. Theory Comput.* 9(1), 794–802.
23. Rawlings ND, Tolle DP, Barrett AJ (2004) Evolutionary families of peptidase inhibitors. *Biochem. J.* 378, 705–716.
24. Rothenberg ME (2009) Biology and treatment of eosinophilic esophagitis. *Gastroenterology* 137(4), 1238–1249.
25. Azouz NP, Rothenberg ME (2019) Mechanisms of gastrointestinal allergic disorders. *J. Clin. Invest.* 129(4), 1523–1534.
26. Brattsand M, Egelrud T (1999) Purification, molecular cloning, and expression of a human stratum corneum trypsin-like serine protease. *J. Biol. Chem.* 274(42), 30033–30040.
27. Chavanas S et al. (2000) Mutations in SPINK5, encoding a serine protease inhibitor, cause Netherton syndrome. *Nat. Genet.* 25, 141–142.
28. Witt H et al. (2000) Mutations in the gene encoding the serine protease inhibitor, Kazal type 1, are associated with chronic pancreatitis. *Nat. Genet.* 25, 213–216.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
