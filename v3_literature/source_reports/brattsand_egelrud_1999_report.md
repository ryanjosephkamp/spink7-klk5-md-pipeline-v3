# Source Report: S-05 — Brattsand & Egelrud 1999

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-05  
**PDF:** `v3_literature/sources/S-05.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** Purification, Molecular Cloning, and Expression of a Human Stratum Corneum Trypsin-like Serine Protease with Possible Function in Desquamation

**Authors:** Maria Brattsand, Torbjörn Egelrud

**Journal:** *The Journal of Biological Chemistry*, 274(42), 30033–30040, 1999

**DOI:** 10.1074/jbc.274.42.30033

**Affiliation:** Department of Public Health and Clinical Medicine, Dermatology and Venereology, Umeå University, Sweden

---

## 2. Summary

This is the original paper identifying and characterizing KLK5 (here named SCTE — stratum corneum tryptic enzyme). The authors purified the enzyme from human plantar stratum corneum, cloned its cDNA, and expressed recombinant protein. The study establishes the primary biochemical identity of KLK5 — its molecular mass, catalytic triad, substrate specificity, glycosylation, tissue expression, and activation mechanism — but does not report quantitative kinetic parameters ($K_m$, $k_{\text{cat}}$, $K_i$).

---

## 3. Identified Features

### Feature 1: KLK5 Molecular Mass

| Parameter | Value |
|-----------|-------|
| **Measurement** | Calculated and observed molecular mass |
| **Value** | Calculated: 25.2 kDa (from sequence); Observed: ~33 kDa (glycosylated, reduced samples) |
| **Feature Type** | Structural (basic characterization) |
| **Pipeline Benchmarkable** | No — but useful for validating MD system setup (correct mass) |

### Feature 2: KLK5 Primary Structure

| Parameter | Value |
|-----------|-------|
| **Measurement** | Protein architecture |
| **Value** | Signal peptide: 29 residues; Propeptide: 37 residues; Active enzyme: 227 residues; 12 cysteines (6 predicted S–S bonds); 4 potential N-glycosylation sites |
| **Feature Type** | Structural |
| **Pipeline Benchmarkable** | No — but essential for correct system setup |

### Feature 3: KLK5 Catalytic Triad

| Parameter | Value |
|-----------|-------|
| **Measurement** | Catalytic triad residues |
| **Value** | His-41, Asp-87, Ser-179 (in mature sequence); Asp-173 at S1 pocket → trypsin-like specificity |
| **Feature Type** | Structural/functional |
| **Pipeline Benchmarkable** | No — but confirms expected binding mode for SPINK7 inhibitory loop interaction |

### Feature 4: KLK5 Tissue Expression

| Parameter | Value |
|-----------|-------|
| **Measurement** | Tissue distribution (RT-PCR) |
| **Value** | Highest in skin; lower in brain, placenta, kidney |
| **Feature Type** | Biological context |
| **Pipeline Benchmarkable** | No |

### Feature 5: Sequence Homology

| Parameter | Value |
|-----------|-------|
| **Measurement** | Amino acid sequence identity to related proteases |
| **Value** | 55% to porcine enamel matrix protease; 50% to TLSP; 49% to neuropsin; 48% to SCCE (KLK7) |
| **Feature Type** | Structural context |
| **Pipeline Benchmarkable** | No |

---

## 4. Overall Assessment

**Usefulness for Benchmarking: Low**

This paper provides the foundational biochemical characterization of KLK5 (SCTE) — its sequence, structure, catalytic mechanism, and tissue expression — but does not report quantitative kinetic or thermodynamic parameters. Its value is structural: confirming the catalytic triad, disulfide bond pattern, and activation mechanism that our MD simulations must correctly represent. No directly benchmarkable features.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
