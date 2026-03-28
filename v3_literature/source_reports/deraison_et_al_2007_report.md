# Source Report: S-09 — Deraison et al. 2007

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-09  
**PDF:** `v3_literature/sources/S-09.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** LEKTI Fragments Specifically Inhibit KLK5, KLK7, and KLK14 and Control Desquamation through a pH-dependent Interaction

**Authors:** Celine Deraison, Chrystelle Bonnart, Frederic Lopez, Celine Besson, Ross Robinson, Arumugam Jayakumar, Fredrik Wagberg, Maria Brattsand, Jean Pierre Hachem, Goran Leonardsson, Alain Hovnanian

**Journal:** *Molecular Biology of the Cell*, 18, 3607–3619, 2007

**DOI:** 10.1091/mbc.E07-02-0124

---

## 2. Summary

Comprehensive study of LEKTI fragment inhibition of KLK5, KLK7, and KLK14 with both functional inhibition ($K_i$) and biophysical binding ($K_D$ by SPR/BIAcore) data. Critically demonstrates pH-dependent binding kinetics — LEKTI binds tightly to KLK5 at physiological pH (7.5) but releases KLK5 at acidic pH (4.5–5.5), providing the molecular mechanism for stratum corneum desquamation regulation.

---

## 3. Identified Features

### Feature 1: LEKTI Domain–KLK5 Inhibition and Binding (**BENCHMARK TARGET**)

| LEKTI Fragment | $K_i$ (nM) | $K_D$ by SPR (M) | $\Delta G_{\text{bind}}$ at 298 K (from $K_D$) |
|----------------|----------|-------------|------|
| D5 | 32.8 | 9.3 × 10⁻¹⁰ | −12.3 kcal/mol |
| D6 | 83.3 | 3.6 × 10⁻⁹ | −11.5 kcal/mol |
| **D8–D11** | **3.7** | **1.1 × 10⁻¹²** | **−16.3 kcal/mol** |
| D9–D15 | 118.7 | 9.5 × 10⁻¹⁰ | −12.3 kcal/mol |

- **Pipeline Benchmarkable:** **YES** — orthogonal $K_i$ and $K_D$ measurements provide primary benchmarks
- **D8–D11/KLK5 $K_D$ = 1.1 pM** — essentially irreversible at neutral pH; exceptionally tight binding (note: $K_D$ by SPR may overestimate affinity for irreversible binders)
- **Proposed 95% CI:** $K_i$: ±3-fold; $K_D$: within order of magnitude for very tight binders
- **Relevant Table:** Table 2

### Feature 2: pH-Dependent Kinetics (**BENCHMARK TARGET**)

| pH | $k_a$ (M⁻¹s⁻¹) | $k_d$ (s⁻¹) |
|----|----------------|-------------|
| 7.5 | ~10⁴ | ~10⁻⁸ (irreversible) |
| 5.5 | — | ~10⁻⁵ |
| 4.5 | ~10³ | ~10⁻⁴ |

- **Pipeline Benchmarkable:** **YES** — pH-dependent binding kinetics can validate constant-pH MD or pH-dependent PMF calculations
- **Feature Type:** Kinetic / thermodynamic

### Feature 3: LEKTI Domain–KLK7 Inhibition

| Fragment | $K_i$ (nM) | $K_D$ (M) |
|----------|----------|---------|
| D5 | 77.2 | 6.2 × 10⁻⁹ |
| D6 | 296.4 | 6.7 × 10⁻⁸ |
| D8–D11 | 34.8 | 3.2 × 10⁻⁸ |

### Feature 4: KLK5 Substrate $K_m$

| Parameter | Value |
|-----------|-------|
| KLK5 $K_m$ for Ile-Pro-Arg-pNA | 0.6 mM |

### Feature 5: LEKTI Domain D1 Non-inhibition

D1 domain shows no inhibition and no binding to any KLK — negative control.

---

## 4. Overall Assessment

**Usefulness for Benchmarking: HIGH**

This is one of the two most important benchmark sources for LEKTI/SPINK5-KLK5 binding (alongside S-08). Provides both $K_i$ (functional) and $K_D$ (biophysical/SPR) measurements for the same interactions, enabling cross-validation. The pH-dependent binding kinetics ($k_a$, $k_d$ at pH 4.5–7.5) are uniquely valuable for validating pH-dependent MD simulations. The D8–D11/KLK5 picomolar binding is the tightest measured LEKTI-KLK interaction.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
