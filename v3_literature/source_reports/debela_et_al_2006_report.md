# Source Report: S-07 — Debela et al. 2006

**Document Class:** V3 Literature Review — Per-Source Analysis Report  
**Source ID:** S-07  
**PDF:** `v3_literature/sources/S-07.pdf`  
**Date Reviewed:** 2026-03-22

---

## 1. Bibliographic Information

**Title:** Specificity Profiling of Seven Human Tissue Kallikreins Reveals Individual Subsite Preferences

**Authors:** Mekdes Debela, Viktor Magdolen, Norman Schechter, Martina Valachova, Friedrich Lottspeich, Charles S. Craik, Youngchool Choe, Wolfram Bode, Peter Goettig

**Journal:** *The Journal of Biological Chemistry*, 281(35), 25678–25688, 2006

**DOI:** 10.1074/jbc.M602372200

---

## 2. Summary

Comprehensive substrate specificity profiling of seven human kallikreins (KLK3/PSA, KLK4, KLK5, KLK6, KLK7, KLK10, KLK11) using positional scanning combinatorial libraries (PS-SCL). For KLK5, the paper reports exclusive P1 Arg selectivity, S2–S4 subsite preferences, and structural analysis of the active-site architecture. Provides cross-comparison between kallikrein family members.

---

## 3. Identified Features

### Feature 1: KLK5 Extended Subsite Preferences (PS-SCL)

| Subsite | Preferences |
|---------|------------|
| **P1** | Exclusively Arg (most specific among all 7 KLKs) |
| **P2** | Ser slightly preferred > Asn, Thr, Ala; Arg and Tyr rejected |
| **P3** | Met, Phe, Tyr accepted; Asp rejected |
| **P4** | Unusual Gly preference, then Tyr, Val, Pro |

- **Pipeline Benchmarkable:** Partially — validates computational substrate/inhibitor docking into KLK5 active site
- **Feature Type:** Structural/selectivity

### Feature 2: KLK5 Active-Site Structural Residues

| Parameter | Value |
|-----------|-------|
| S1 pocket key residues | Asp-189, Ser-190, Gly-216, Gly-226 |
| **Feature Type** | Structural |
| **Pipeline Benchmarkable** | Yes — validates correct active-site representation in MD models |

### Feature 3: PDB Codes Referenced

Related kallikrein structures: 1SPJ (KLK1), 1LO6 (KLK6), 2PKA, 2KAI (porcine KLK), 1GVZ, 1TON, 1AO5

---

## 4. Overall Assessment

**Usefulness for Benchmarking: Medium**

Provides detailed substrate specificity profiles for KLK5 and comparative data across seven kallikreins. The subsite preference data can validate in silico substrate docking. Structural analysis of S1–S4 pockets is useful for active-site modeling. No direct $K_i$ or $K_d$ measurements for protease-inhibitor interactions.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
