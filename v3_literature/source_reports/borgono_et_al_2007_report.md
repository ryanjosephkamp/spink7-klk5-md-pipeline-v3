# Source Report: S-08 — Borgoño et al. 2007

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-08  
**PDF:** `v3_literature/sources/S-08.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** A Potential Role for Multiple Tissue Kallikrein Serine Proteases in Epidermal Desquamation

**Authors:** Carla A. Borgoño, Iacovos P. Michael, Nahoko Komatsu, Arumugam Jayakumar, Ravi Kapadia, Gary L. Clayman, Georgia Sotiropoulou, Eleftherios P. Diamandis

**Journal:** *The Journal of Biological Chemistry*, 282(6), 3640–3652, 2007

**DOI:** 10.1074/jbc.M607567200

---

## 2. Summary

Comprehensive study of LEKTI (SPINK5) domain fragment inhibition of multiple kallikreins. Reports $K_i$ values for four recombinant LEKTI fragments (domains 1–6, 6–9', 9–12, 12–15) against KLK5, KLK6, KLK13, and KLK14. Also identifies DSG1 cleavage sites by multiple KLKs. This is a primary benchmark source for SPINK5-KLK5 binding.

---

## 3. Identified Features

### Feature 1: LEKTI Fragment–KLK5 $K_i$ Values (**BENCHMARK TARGET**)

| LEKTI Fragment | $K_i$ (nM) | Mechanism | $\Delta G_{\text{bind}}$ at 298 K |
|----------------|----------|-----------|-------------------------------|
| rLEKTI(1–6) | **2.35 ± 0.22** | Mixed | −11.8 kcal/mol |
| rLEKTI(6–9') | **4.68 ± 0.66** | Mixed | −11.4 kcal/mol |
| rLEKTI(9–12) | **2.75 ± 0.24** | Mixed | −11.7 kcal/mol |
| rLEKTI(12–15) | **21.80 ± 2.40** | Mixed | −10.5 kcal/mol |

- **Pipeline Benchmarkable:** **YES** — direct $K_i$ values convertible to $\Delta G_{\text{bind}}$ for SPINK5-KLK5 validation
- **Proposed 95% CI:** ±3-fold for enzyme inhibition assays → 0.8–7 nM for domains 1–6
- **Relevant Table:** Table 1

### Feature 2: LEKTI Fragment–KLK14 $K_i$ Values

| LEKTI Fragment | $K_i$ (nM) |
|----------------|----------|
| rLEKTI(1–6) | 0.22 ± 0.021 |
| rLEKTI(6–9') | 3.41 ± 0.80 |
| rLEKTI(9–12) | 10.26 ± 1.25 |
| rLEKTI(12–15) | NI |

- **Pipeline Benchmarkable:** Secondary target for selectivity validation

### Feature 3: LEKTI Fragment–KLK6 and KLK13 $K_i$ Values

| LEKTI Fragment | KLK6 $K_i$ (nM) | KLK13 $K_i$ (nM) |
|----------------|-----------|------------|
| rLEKTI(1–6) | 13.08 ± 2.36 | 24.13 ± 3.82 |
| rLEKTI(6–9') | 47.58 ± 1.40 | 222.12 ± 9.51 |
| rLEKTI(9–12) | 195.32 ± 11.66 | 408.63 ± 16.47 |
| rLEKTI(12–15) | NI | NI |

### Feature 4: SLPI and Elafin Non-inhibition

Neither SLPI nor elafin inhibited any KLK tested — useful selectivity control.

### Feature 5: DSG1 Cleavage Sites

| KLK | DSG1 Cleavage Sites |
|-----|---------------------|
| KLK5 | Arg146–Val147, Arg242–Asp243, Arg422–Tyr423 |
| KLK6 | Lys368–Ala369, Arg422–Tyr423 |
| KLK14 | Tyr528–Ser529 |
| KLK1 | Lys197–Ile198, Lys445–Asn446 |

---

## 4. Overall Assessment

**Usefulness for Benchmarking: HIGH**

Primary source for LEKTI (SPINK5)–KLK5 inhibition constants. The four LEKTI fragment $K_i$ values (2.35–21.80 nM against KLK5) provide direct benchmarks for MD free-energy calculations of the SPINK5-KLK5 system, which is a homologous system to our primary target SPINK7-KLK5. The cross-KLK selectivity data provides additional validation dimensions.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
