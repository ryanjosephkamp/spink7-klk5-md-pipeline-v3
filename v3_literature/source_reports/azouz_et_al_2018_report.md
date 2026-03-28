# Source Report: S-02 — Azouz et al. 2018

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-02  
**PDF:** `v3_literature/sources/S-02.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** The antiprotease SPINK7 serves as an inhibitory checkpoint for esophageal epithelial inflammatory responses

**Authors:** Netali Paz Azouz, Mark A. Ynga-Durand, Julie M. Caldwell, Ayele Jain, Manisha C. Barski, Marc Rochman, Huan He, Leah C. Kottyan, Lisa J. Martin, Daniela M. Fischesser, Lydia M. Rosen, Jeffrey T. Knife, Jocelyn M. Biagini Myers, Vincent A. Mukkada, Philip E. Putnam, Gurjit K. Khurana Hershey, Michael K. Mingler, Melissa Evans, Marc E. Rothenberg

**Journal:** *Science Translational Medicine*, 10(444), eaap9736, 2018

**DOI:** 10.1126/scitranslmed.aap9736

**Affiliation:** Division of Allergy and Immunology, Cincinnati Children's Hospital Medical Center, Department of Pediatrics, University of Cincinnati College of Medicine

---

## 2. Summary

This is the landmark study identifying SPINK7 as a critical antiprotease checkpoint in esophageal epithelial inflammatory responses and the central motivation for our entire MD pipeline. The paper demonstrates that SPINK7 is markedly down-regulated in EoE and that its loss is sufficient to impair the epithelial barrier, increase protease activity (uPA, MMPs), release proinflammatory cytokines, and activate eosinophils. Genetic epistasis between TSLP and PLAU further supports the biological relevance of the pathway.

The study is primarily cell biology and translational immunology — it establishes the biological significance of SPINK7s antiprotease role but reports functional cellular readouts (TEER, cytokine levels, protease activity fold-changes) rather than purified-system biophysical measurements (Ki, Kd, kcat). The quantitative data it provides are biological effect sizes in complex cellular or tissue systems, not molecular-level thermodynamic/kinetic parameters directly comparable to MD simulation outputs.

---

## 3. Identified Features

### Feature 1: SPINK7 mRNA Down-regulation in EoE

| Parameter | Value |
|-----------|-------|
| **Measurement** | SPINK7 mRNA fold-change in active EoE vs. control |
| **Value** | 16-fold decrease (P = 3 × 10⁻⁸) |
| **System** | Human esophageal biopsies (n ≥ 10 per group) |
| **Method** | Microarray + qPCR validation |
| **Feature Type** | Biological context (gene expression) |
| **Pipeline Benchmarkable** | No — gene expression, not a biophysical parameter |
| **Relevant Figure/Table** | Fig. 1A |

### Feature 2: SPINK5 mRNA Down-regulation in EoE

| Parameter | Value |
|-----------|-------|
| **Measurement** | SPINK5 mRNA fold-change in active EoE vs. control |
| **Value** | 1.9-fold decrease (P = 0.0005) |
| **System** | Human esophageal biopsies |
| **Method** | Microarray |
| **Feature Type** | Biological context (gene expression) |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 1A |

### Feature 3: SPINK7-SPINK5 Kazal Domain Sequence Identity

| Parameter | Value |
|-----------|-------|
| **Measurement** | Sequence identity between SPINK7 kazal domain and SPINK5 domain 15 |
| **Value** | 35% identity |
| **System** | Protein sequence alignment |
| **Method** | Sequence comparison |
| **Feature Type** | Structural context |
| **Pipeline Benchmarkable** | No — sequence identity, not a simulation target; but useful for homology modeling rationale |
| **Relevant Figure/Table** | Fig. 1B |

### Feature 4: TEER Reduction upon SPINK7 Silencing

| Parameter | Value |
|-----------|-------|
| **Measurement** | Transepithelial electrical resistance (TEER) fold-change |
| **Value** | ~2.3-fold decrease (P < 0.0001) upon SPINK7 shRNA silencing in differentiated EPC2 cells |
| **System** | EPC2 esophageal epithelial cells, ALI differentiation |
| **Method** | Resistance measurement across cell monolayer |
| **Feature Type** | Biological effect (barrier function) |
| **Pipeline Benchmarkable** | No — cellular phenotype, not molecular measurement |
| **Relevant Figure/Table** | Fig. 2B |

### Feature 5: TEER Reduction upon SPINK5 Silencing

| Parameter | Value |
|-----------|-------|
| **Measurement** | TEER fold-change upon SPINK5 silencing |
| **Value** | ~2.5-fold decrease (P < 0.0001) |
| **System** | EPC2 esophageal epithelial cells, ALI differentiation |
| **Method** | Resistance measurement |
| **Feature Type** | Biological effect (barrier function) |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 2B |

### Feature 6: TEER Reduction in CRISPR/Cas9 SPINK7 KO

| Parameter | Value |
|-----------|-------|
| **Measurement** | TEER fold-change in SPINK7 CRISPR/Cas9 knockout |
| **Value** | 2-fold decrease (P < 0.0001) |
| **System** | EPC2 SPINK7-KO cells |
| **Method** | CRISPR/Cas9 gene deletion + TEER measurement |
| **Feature Type** | Biological effect (barrier function) |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 5A |

### Feature 7: uPA Activity Increase upon SPINK7 Silencing

| Parameter | Value |
|-----------|-------|
| **Measurement** | uPA activity fold-change in SPINK7-silenced cells |
| **Value** | 1.8-fold increase (P = 0.013) |
| **System** | EPC2 differentiated cells, SPINK7 shRNA |
| **Method** | uPA activity assay (supernatants) |
| **Feature Type** | Enzymatic activity (indirect) |
| **Pipeline Benchmarkable** | Indirect — demonstrates SPINK7 inhibits uPA but does not report Ki/Kd |
| **Relevant Figure/Table** | Fig. 3A |

### Feature 8: uPA Activity Increase in EoE Biopsies

| Parameter | Value |
|-----------|-------|
| **Measurement** | uPA activity in EoE biopsies vs. control |
| **Value** | 10-fold increase (P = 0.043) |
| **System** | Human esophageal biopsies (EoE vs. control) |
| **Method** | uPA activity assay on tissue homogenates |
| **Feature Type** | Enzymatic activity (tissue-level) |
| **Pipeline Benchmarkable** | No — tissue-level, not purified system |
| **Relevant Figure/Table** | Fig. 3B |

### Feature 9: uPA Inhibition by SPINK7 KO Cell Supernatants

| Parameter | Value |
|-----------|-------|
| **Measurement** | Purified uPA activity inhibition by cell supernatants |
| **Value** | Control supernatants inhibited uPA by ~70%; SPINK7-KO supernatants showed increased uPA activity |
| **System** | Purified uPA + conditioned media from SPINK7-KO vs. control EPC2 cells |
| **Method** | uPA activity assay with exogenous purified uPA |
| **Feature Type** | Enzymatic inhibition (semi-quantitative) |
| **Pipeline Benchmarkable** | Partially — demonstrates SPINK7 is a uPA inhibitor, but no Ki value reported |
| **Relevant Figure/Table** | Fig. 5D |

### Feature 10: MMP Activity Increase in SPINK7-KO Cells

| Parameter | Value |
|-----------|-------|
| **Measurement** | Total MMP activity increase |
| **Value** | 23% increase (P = 0.008) |
| **System** | SPINK7-KO EPC2 cells |
| **Method** | Fluorogenic MMP substrate assay |
| **Feature Type** | Enzymatic activity |
| **Pipeline Benchmarkable** | No — MMP activity, not SPINK7-protease binding |
| **Relevant Figure/Table** | Fig. 5B |

### Feature 11: MMP9 Activity Increase in SPINK7-KO Cells

| Parameter | Value |
|-----------|-------|
| **Measurement** | MMP9 activity fold-change |
| **Value** | 2-fold increase (P = 0.0048) |
| **System** | SPINK7-KO EPC2 cells |
| **Method** | Gelatin zymography |
| **Feature Type** | Enzymatic activity |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 5C |

### Feature 12: IL-8 Release upon SPINK7 Silencing

| Parameter | Value |
|-----------|-------|
| **Measurement** | IL-8 release fold-change |
| **Value** | 12-fold increase (P = 0.03) |
| **System** | SPINK7-silenced EPC2 differentiated cells |
| **Method** | Cytokine multiplex assay |
| **Feature Type** | Biological effect (cytokine release) |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 4B |

### Feature 13: A1AT Rescue of Barrier Function

| Parameter | Value |
|-----------|-------|
| **Measurement** | Alpha-1-antitrypsin (A1AT) reversal of SPINK7 KO barrier defect |
| **Value** | A1AT administration reversed pathologic features of SPINK7 silencing (TEER recovery, reduced FITC-dextran flux) |
| **System** | SPINK7-silenced EPC2 cells + exogenous A1AT |
| **Method** | TEER + permeability assays |
| **Feature Type** | Biological rescue (therapeutic proof-of-concept) |
| **Pipeline Benchmarkable** | No — demonstrates protease inhibitor therapeutic concept but no molecular parameters |
| **Relevant Figure/Table** | Fig. 6 |

### Feature 14: Genetic Epistasis — TSLP × PLAU

| Parameter | Value |
|-----------|-------|
| **Measurement** | Odds ratio for EoE susceptibility with combined TSLP + PLAU minor alleles |
| **Value** | OR = 2.73 (P = 0.0003) |
| **System** | Human genetic study (725 EoE cases, 412 controls) |
| **Method** | Logistic regression with interaction term, custom Illumina GoldenGate SNP chip |
| **Feature Type** | Genetic epidemiology |
| **Pipeline Benchmarkable** | No |
| **Relevant Figure/Table** | Fig. 8B |

---

## 4. Benchmarkability Assessment

Despite reporting 14 quantitative features, **none are directly benchmarkable** by the SPINK7-KLK5 MD pipeline in the conventional sense of thermodynamic ($K_i$, $\Delta G_{\text{bind}}$) or kinetic ($k_{\text{on}}$, $k_{\text{off}}$, $k_{\text{cat}}$) parameters. The measurements are cellular/tissue-level readouts that demonstrate SPINK7's biological function but do not provide purified-system biophysical data.

However, this paper provides **critical indirect value**:
1. Confirms SPINK7 is a functional serine protease inhibitor of uPA (Feature 9: supernatant inhibition assay)
2. Establishes that 35% sequence identity between SPINK7 kazal and SPINK5 domain 15 justifies homology modeling approaches
3. Provides the biological rationale for why accurate MD simulation of SPINK7 inhibitory interactions matters

**Features extractable for pipeline context (not benchmarking):**
- SPINK7 is a functional uPA inhibitor (demonstrated but Ki not measured)
- The ~70% inhibition of purified uPA by control cell supernatants (Feature 9) could theoretically constrain expected Ki ranges, but the assay conditions (crude supernatant, unknown SPINK7 concentration) preclude quantitative extraction

---

## 5. Overall Assessment

**Usefulness for Benchmarking: Medium**

This source is the foundational paper for our entire project — it establishes SPINK7's antiprotease role in EoE and motivates the pipeline. However, as a cell biology and translational immunology study, it reports functional cellular readouts rather than purified-system biophysical measurements. The paper demonstrates that SPINK7 inhibits uPA and that loss of SPINK7 causes barrier dysfunction, protease activation, and inflammation, but it does not measure binding constants, enzymatic kinetics, or structural features that can be directly compared to MD simulation outputs.

**Key value:** Biological motivation, pathway validation, homology rationale (35% identity between SPINK7 and SPINK5 kazal domains). The paper points toward purified biochemical studies of SPINK7/SPINK5-protease interactions as the next step, which is what our pipeline aims to computationally model.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
