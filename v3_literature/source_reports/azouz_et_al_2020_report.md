# Source Report: S-03 — Azouz et al. 2020

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-03  
**PDF:** `v3_literature/sources/S-03.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** Functional role of kallikrein 5 and proteinase-activated receptor 2 in eosinophilic esophagitis

**Authors:** Nurit P. Azouz, Andrea M. Klingler, Purnima Pathre, John A. Besse, Netali Ben Baruch-Morgenstern, Adina Y. Ballaban, Garrett A. Osswald, Michael Brusilovsky, Jeff E. Habel, Julie M. Caldwell, Mario A. Ynga-Durand, Pablo J. Abonia, Yueh-Chiang Hu, Ting Wen, Marc E. Rothenberg

**Journal:** *Science Translational Medicine*, 12(545), eaaz7773, 2020

**DOI:** 10.1126/scitranslmed.aaz7773

**Affiliation:** Division of Allergy and Immunology, Cincinnati Children's Hospital Medical Center, Department of Pediatrics, University of Cincinnati College of Medicine

---

## 2. Summary

This is the critical follow-up to S-02 (Azouz et al. 2018) that directly identifies KLK5 as a functional target of SPINK7, making it the single most important source for our MD pipeline benchmarking. The paper reports the inhibition constant ($K_i$) for SPINK7 against KLK5, structural docking predictions for the SPINK7-KLK5 complex, and extensive in vivo validation using Klk5−/− mice and PAR2 antagonism in a murine EoE model. It also demonstrates that A1AT delivery to the esophagus attenuates experimental EoE.

---

## 3. Identified Features

### Feature 1: SPINK7-KLK5 Inhibition Constant (**BENCHMARK TARGET**)

| Parameter | Value |
|-----------|-------|
| **Measurement** | Inhibition constant ($K_i$) of SPINK7 against KLK5 |
| **Value** | $K_i = 132$ nM |
| **System** | Purified recombinant SPINK7 + recombinant KLK5 + fluorogenic substrate |
| **Method** | Dose-response curve fitted to Morrison equation (tight-binding inhibitor model) |
| **Replicates** | Three independent experiments performed in duplicates |
| **Feature Type** | Thermodynamic (inhibition constant) |
| **Pipeline Benchmarkable** | **YES** — This is the primary benchmark for the SPINK7-KLK5 MD pipeline. The Ki can be converted to $\Delta G_{\text{bind}}$ via $\Delta G = RT \ln(K_i)$, yielding $\Delta G \approx -9.4$ kcal/mol at 298 K |
| **Proposed 95% CI** | Literature tight-binding Ki assays typically have 2-3× uncertainty; estimated range: 50–350 nM ($\Delta G$: −10.0 to −8.8 kcal/mol) |
| **Relevant Figure/Table** | Fig. 1B |

### Feature 2: SPINK7-KLK12 Inhibition Constant

| Parameter | Value |
|-----------|-------|
| **Measurement** | Inhibition constant ($K_i$) of SPINK7 against KLK12 |
| **Value** | $K_i = 1500$ nM |
| **System** | Purified recombinant SPINK7 + recombinant KLK12 + fluorogenic substrate |
| **Method** | Morrison equation |
| **Feature Type** | Thermodynamic (inhibition constant) |
| **Pipeline Benchmarkable** | Yes (secondary target) — selectivity ratio: SPINK7 is 11.5-fold more potent toward KLK5 than KLK12 |
| **Relevant Figure/Table** | Fig. 1B, fig. S1B |

### Feature 3: SPINK7 Non-inhibition of KLK7, KLK11, KLK13

| Parameter | Value |
|-----------|-------|
| **Measurement** | SPINK7 inhibitory activity against KLK7, KLK11, KLK13 |
| **Value** | No detectable inhibition |
| **System** | Purified recombinant proteins |
| **Feature Type** | Selectivity profile (negative controls) |
| **Pipeline Benchmarkable** | Yes — MD simulations should predict weak/no binding for these targets |
| **Relevant Figure/Table** | fig. S1, B and C |

### Feature 4: ClusPro Docking — SPINK7-KLK5 Complex Structure

| Parameter | Value |
|-----------|-------|
| **Measurement** | Computational docking of SPINK7-KLK5 complex |
| **Value** | Inhibitory loop of SPINK7 positioned ~2.5 Å from serine residue in KLK5 catalytic triad |
| **PDB Inputs** | SPINK7: PDB 2LEO (chain A); KLK5: PDB 2PSX (chain A) |
| **Method** | ClusPro (rotational screening → ~10⁹ conformations → 1000 top scored → clustered at Cα RMSD < 9 Å → largest cluster). Two additional docking algorithms confirmed same binding mode |
| **Feature Type** | Structural (computational docking) |
| **Pipeline Benchmarkable** | **YES** — provides starting complex geometry and PDB identifiers for MD simulations |
| **Relevant Figure/Table** | Fig. 1C, fig. S1D |

### Feature 5: Trypsin-like Activity Increase upon SPINK7 Depletion

| Parameter | Value |
|-----------|-------|
| **Measurement** | Trypsin-like activity fold-change in SPINK7-depleted cells |
| **Value** | 2-fold increase (P ≤ 0.004) in EPC2 cells; (P = 0.0007) in primary cells; (P = 0.0001) in CRISPR-Cas9 KO |
| **System** | EPC2 cells + primary esophageal cells, ALI differentiation |
| **Feature Type** | Enzymatic activity (cellular) |
| **Pipeline Benchmarkable** | Indirect — confirms functional inhibition but cellular readout |
| **Relevant Figure/Table** | Fig. 1A |

### Feature 6: KLK5 Overexpression — TEER Reduction

| Parameter | Value |
|-----------|-------|
| **Measurement** | TEER reduction upon KLK5 overexpression |
| **Value** | P = 0.0014 (significant decrease) |
| **System** | EPC2 cells overexpressing KLK5 |
| **Feature Type** | Biological effect |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 1H |

### Feature 7: KLK5 Overexpression — DSG1 Decrease

| Parameter | Value |
|-----------|-------|
| **Measurement** | DSG1 protein decrease upon KLK5 overexpression |
| **Value** | 33% decrease (P = 0.01) |
| **System** | KLK5-overexpressing EPC2 cells |
| **Feature Type** | Biological effect |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 1E |

### Feature 8: Klk5−/− In Vivo — Esophageal Eosinophil Reduction

| Parameter | Value |
|-----------|-------|
| **Measurement** | Eosinophil reduction in Klk5−/− mice after OVA challenge |
| **Value** | 3-fold decrease vs. Klk5+/+ OVA-treated (P = 0.0003) |
| **System** | Murine EoE model (OVA sensitization + intranasal challenge) |
| **Feature Type** | In vivo phenotype |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 2C |

### Feature 9: Klk5−/− In Vivo — Esophageal Protease Activity

| Parameter | Value |
|-----------|-------|
| **Measurement** | Trypsin-like activity increase in esophagus after OVA |
| **Value** | 50% increase in Klk5+/+ OVA vs. saline (P = 0.025); no increase in Klk5−/− |
| **System** | Murine esophageal tissue |
| **Feature Type** | In vivo enzymatic activity |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 2B |

### Feature 10: SPINK5-KLK5 Inhibition (Confirmatory)

| Parameter | Value |
|-----------|-------|
| **Measurement** | SPINK5 inhibits KLK5 activity in vitro |
| **Value** | Confirmed (SPINK5 at 500 nM abolished KLK5-mediated MUC4 degradation) |
| **System** | Purified SPINK5 + KLK5 + esophageal biopsy lysates |
| **Feature Type** | Enzymatic inhibition |
| **Pipeline Benchmarkable** | Partially — confirms SPINK5-KLK5 inhibitory interaction but Ki not reported in this paper |
| **Relevant Figure/Table** | Fig. 5F, fig. S4 |

### Feature 11: MUC4 Degradation Sensitivity in EoE

| Parameter | Value |
|-----------|-------|
| **Measurement** | MUC4 degradation by KLK5 in EoE vs. control biopsies |
| **Value** | 2-fold increase in degradation (20% vs. 40%; P = 0.03) in EoE biopsies with 10 nM KLK5 |
| **System** | Human esophageal biopsy lysates + recombinant KLK5 |
| **Feature Type** | Substrate processing |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 5D |

### Feature 12: A1AT-KLK5 Inhibition

| Parameter | Value |
|-----------|-------|
| **Measurement** | A1AT blocks KLK5 activity in vitro |
| **Value** | Dose-dependent inhibition confirmed (fig. S3A) |
| **System** | Purified A1AT + KLK5 |
| **Feature Type** | Enzymatic inhibition |
| **Pipeline Benchmarkable** | No (A1AT is a serpin, not a kazal-type inhibitor; different mechanism) |
| **Relevant Figure/Table** | fig. S3A |

### Feature 13: A1AT In Vivo — Esophageal Eosinophil Reduction

| Parameter | Value |
|-----------|-------|
| **Measurement** | Eosinophil reduction after A1AT delivery |
| **Value** | 3.6-fold decrease (P = 0.0005) |
| **System** | Murine EoE model + intraperitoneal A1AT |
| **Feature Type** | In vivo phenotype |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 4F |

### Feature 14: PAR2 Antagonist — Eosinophil Reduction

| Parameter | Value |
|-----------|-------|
| **Measurement** | Eosinophil reduction after ENMD-1068 (PAR2 antagonist) |
| **Value** | 2.4-fold decrease (P = 0.009) |
| **System** | Murine EoE model + ENMD-1068 |
| **Feature Type** | In vivo phenotype |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 7E |

### Feature 15: F2RL1 (PAR2) Overexpression in EoE

| Parameter | Value |
|-----------|-------|
| **Measurement** | F2RL1 mRNA fold-change in EoE vs. control |
| **Value** | 1.85-fold increase (P = 0.0035) |
| **System** | Human esophageal biopsies (bulk RNA-seq) |
| **Feature Type** | Gene expression |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 6C |

---

## 4. Benchmarkability Assessment

**This is the most important source for pipeline benchmarking.** It provides:

1. **Primary benchmark:** $K_i = 132$ nM for SPINK7-KLK5 (Feature 1) → $\Delta G_{\text{bind}} \approx -9.4$ kcal/mol
2. **Selectivity benchmark:** $K_i = 1500$ nM for SPINK7-KLK12 → $\Delta G_{\text{bind}} \approx -7.9$ kcal/mol; SPINK7 does not inhibit KLK7, KLK11, KLK13
3. **Structural starting point:** PDB files 2LEO (SPINK7) and 2PSX (KLK5); ClusPro docking shows inhibitory loop–catalytic triad interaction at ~2.5 Å
4. **SPINK5 cross-validation:** SPINK5 inhibits KLK5 (demonstrated at 500 nM), providing a secondary system for MD validation

**Directly benchmarkable features: 3** (Features 1, 2, and 4)

---

## 5. Overall Assessment

**Usefulness for Benchmarking: HIGH**

This is the single most important paper for our pipeline. The SPINK7-KLK5 $K_i = 132$ nM is the primary experimental benchmark against which all MD free-energy calculations will be validated. The PDB identifiers (2LEO, 2PSX) provide the structural inputs for simulation. The selectivity data (KLK12 Ki, non-inhibition of KLK7/11/13) provide additional validation targets. The docking prediction (~2.5 Å from catalytic triad) provides a structural constraint for the predicted binding mode.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
