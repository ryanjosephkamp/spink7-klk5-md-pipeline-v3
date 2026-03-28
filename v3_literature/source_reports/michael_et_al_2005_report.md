# Source Report: S-06 — Michael et al. 2005

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-06  
**PDF:** `v3_literature/sources/S-06.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** Biochemical and Enzymatic Characterization of Human Kallikrein 5 (hK5), a Novel Serine Protease Potentially Involved in Cancer Progression

**Authors:** Iacovos P. Michael, Georgia Sotiropoulou, Georgios Pampalakis, Angeliki Magklara, Manik Ghosh, Greg Wasney, Eleftherios P. Diamandis

**Journal:** *The Journal of Biological Chemistry*, 280(15), 14628–14635, 2005

**DOI:** 10.1074/jbc.M408132200

---

## 2. Summary

First comprehensive enzymatic characterization of human KLK5. Reports Michaelis-Menten kinetics for seven fluorogenic substrates, serpin inhibition constants, optimal pH, and molecular mass characterization. Establishes KLK5 as a trypsin-like serine protease with strong Arg specificity at P1 and optimal activity at pH 8.0.

---

## 3. Identified Features

### Feature 1: KLK5 Substrate Kinetics (Benchmark Target)

| Substrate | $K_m$ (mM) | $k_{\text{cat}}$ (min⁻¹) | $k_{\text{cat}}/K_m$ (mM⁻¹min⁻¹) |
|-----------|-----------|-----------------|--------------------------------|
| Boc-Val-Pro-Arg-AMC | 0.20 ± 0.01 | 196.76 | 946.45 |
| Boc-Phe-Ser-Arg-AMC | 0.19 ± 0.01 | 169.50 | 877.37 |
| Boc-Gln-Ala-Arg-AMC | 0.61 ± 0.03 | 106.97 | 175.33 |
| Boc-Leu-Lys-Arg-AMC | 1.01 ± 0.10 | 49.38 | 48.89 |
| Tos-Gly-Pro-Arg-AMC | 1.69 ± 0.32 | 20.77 | 12.22 |
| Boc-Val-Leu-Lys-AMC | 0.64 ± 0.17 | 4.56 | 7.07 |

- **Pipeline Benchmarkable:** Partially — validates active-site substrate preferences for computational docking; not direct protease-inhibitor binding
- **Relevant Table:** Table I

### Feature 2: Serpin Inhibition of KLK5

| Inhibitor | $K_I$ (µM) | $k_{-2}$ (min⁻¹) |
|-----------|-----------|-----------------|
| α₂-Antiplasmin | 1.07 | 1.09 |
| Antithrombin III | 9.1 | 0.39 |
| α₁-Antitrypsin | No inhibition | — |
| α₁-Antichymotrypsin | No inhibition | — |

- **Pipeline Benchmarkable:** No — serpin mechanism differs from kazal-type inhibition
- **Relevant Table:** Table II

### Feature 3: Optimal pH

| Parameter | Value |
|-----------|-------|
| **Optimal pH** | 8.0 (phosphate buffer) |
| **Pipeline Benchmarkable** | Indirectly — pH for MD simulation conditions |

---

## 4. Overall Assessment

**Usefulness for Benchmarking: Medium**

Provides the first quantitative kinetic characterization of KLK5 substrate preferences and confirms trypsin-like specificity. The $K_m$ and $k_{\text{cat}}$ values are useful for validating computational substrate docking but are not direct protease-inhibitor binding benchmarks. No SPINK7/SPINK5 inhibition data.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
