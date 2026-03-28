# Source Report: S-01 — Rothenberg 2009

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-01  
**PDF:** `v3_literature/sources/S-01.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** Biology and Treatment of Eosinophilic Esophagitis

**Authors:** Marc E. Rothenberg

**Journal:** *Gastroenterology*, 137(4), 1238–1249, 2009

**DOI:** 10.1053/j.gastro.2009.07.007

**Affiliation:** Division of Allergy and Immunology, Department of Pediatrics, Cincinnati Children's Hospital Medical Center, University of Cincinnati College of Medicine, Cincinnati, Ohio

---

## 2. Summary

This article is a comprehensive narrative review of Eosinophilic Esophagitis (EoE), authored by the leading authority in the field. The review covers the epidemiology, symptoms, histopathology, disease pathogenesis, genetics, the role of eosinophils, and therapeutic approaches for EoE. The article establishes the biological framework in which the SPINK7-KLK5 interaction operates: EoE is an antigen-driven, Th2-mediated inflammatory disease characterized by eosinophil accumulation in the esophagus, IL-13-induced epithelial gene expression changes (including eotaxin-3 overexpression and filaggrin down-regulation), and impaired esophageal barrier function.

The review does not itself report original experimental measurements of protease-inhibitor interactions, binding affinities, or enzymatic kinetics. Its value to our benchmarking project lies in (a) establishing the disease context that motivates the entire SPINK7-KLK5 pipeline, (b) summarizing the molecular pathways that connect SPINK7 loss to downstream pathology, and (c) providing references to primary experimental studies that may contain quantitative data. The key mechanistic insight relevant to our pipeline is that IL-13 transcriptionally silences SPINK7 in esophageal epithelium (described more fully in the companion paper by Azouz et al. 2018, S-02), leading to unregulated KLK5 protease activity and subsequent barrier disruption through PAR2 activation.

---

## 3. Identified Features

### Features Not Found / Not Applicable

This source is a narrative review article and does not contain original quantitative experimental values suitable for direct benchmarking of the SPINK7-KLK5 MD pipeline. Specifically:

- **No binding affinity data** ($K_d$, $K_i$, $\Delta G_{\text{bind}}$) are reported for any protease-inhibitor pair.
- **No enzymatic kinetic parameters** ($K_m$, $k_{\text{cat}}$, $k_{\text{on}}$, $k_{\text{off}}$) are reported.
- **No structural data** (buried surface area, interface contacts, RMSD) are reported.
- **No dynamic data** (B-factors, RMSF, conformational timescales) are reported.
- **No mutagenesis data** ($\Delta\Delta G$) are reported.

The article references many primary studies but presents their findings qualitatively or semi-quantitatively (e.g., eosinophil counts per high-power field for histopathological diagnosis, prevalence estimates for epidemiology). These clinical measurements are not directly benchmarkable by our molecular dynamics pipeline.

The one quantitative framework mentioned that has indirect relevance is the sibling recurrence risk ratio ($\lambda_S \approx 80$ for EoE; page 5, paragraph on genetics), which indicates a strong genetic component. However, this is a population-level statistic, not a molecular-level measurement amenable to MD simulation.

---

## 4. Overall Assessment

**Usefulness for Benchmarking: Low**

This review provides essential disease context and biological motivation for the SPINK7-KLK5 pipeline project but does not contain quantitative experimental benchmarks usable for validation of MD simulation predictions. Its primary value is contextual: it establishes why the SPINK7-KLK5 interaction matters clinically and provides the disease framework within which our pipeline's molecular predictions are interpreted. The article is appropriately included in the source list for completeness of the biological narrative, but it contributes zero directly benchmarkable features to the V3 experimental validation plan.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
