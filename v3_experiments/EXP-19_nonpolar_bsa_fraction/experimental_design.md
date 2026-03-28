# EXP-19: Nonpolar Fraction of Buried Surface Area

**Experiment ID:** EXP-19  
**Feature ID:** F-19 (benchmarks.md)  
**Category:** Structural  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the fraction of buried surface area (BSA) at the protease-inhibitor interface that is nonpolar. For BPTI-trypsin, the nonpolar fraction is 61% (Krowarsch et al. 2003, Janin & Chothia 1990), with a 95% CI of [55%, 67%]. This is higher than the average protein surface (~57% nonpolar), reflecting the hydrophobic contribution to binding. The pipeline computes this from the per-atom SASA decomposition used in EXP-16.

---

## 2. Hypothesis

**H₁:** The nonpolar fraction of BSA in BPTI-trypsin will be 61 ± 6% (95% CI: [55%, 67%]).

**H₂:** The nonpolar fraction will be stable throughout the trajectory (CV < 5%), reflecting a persistent interface composition.

---

## 3. Background and Rationale

The hydrophobic effect is the dominant driving force for most protein-protein associations. A high nonpolar BSA fraction indicates substantial hydrophobic burial upon binding, contributing favorable ΔG through the release of ordered water. The 61% value for BPTI-trypsin is intermediate between purely hydrophobic interfaces (~70%) and charge-dominated interfaces (~50%), consistent with the mixed character of protease-inhibitor binding (hydrophobic core + electrostatic rim).

---

## 4. Experimental Protocol

### 4.1 System

BPTI-trypsin (PDB 2PTC) from EXP-04 production trajectory. Shared data with EXP-16 (BSA calculation).

### 4.2 SASA Decomposition

Atom classification:
- **Nonpolar:** Carbon atoms (C, CA, CB, CG, CD, CE, CZ, etc.)
- **Polar:** Nitrogen (N, NE, NZ, NH1, NH2, etc.), Oxygen (O, OG, OH, OD, OE, etc.), Sulfur (SG, SD)

### 4.3 Analysis Protocol

For each frame (last 50 ns at 100 ps intervals = 500 frames):

1. Compute total BSA = SASA(free protease) + SASA(free inhibitor) − SASA(complex).
2. Compute nonpolar BSA: use only nonpolar (carbon) atoms in SASA calculation.
3. Compute polar BSA: total BSA − nonpolar BSA.
4. Nonpolar fraction = BSA_nonpolar / BSA_total × 100%.

### 4.4 Production MD Parameters (from config.py)

- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P
- Probe radius: 1.4 Å
- Duration: 100 ns; analysis: last 50 ns
- Temperature: 310 K
- Save interval: 10 ps

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Nonpolar Fraction |
|---------------|-------------------|
| **PASS** | [55%, 67%] |
| **MARGINAL** | [48%, 74%] |
| **FAIL** | Outside marginal range |

---

## 5. Control Conditions

### 5.1 Crystal Structure Reference

Compute nonpolar fraction directly from PDB 2PTC (single frame). Should agree with literature value of 61%.

### 5.2 Random Surface Control

Compute nonpolar fraction of the total protein surface (not just the interface). Expected: ~57% for average protein surface, confirming the interface is enriched in nonpolar contacts.

### 5.3 Cross-System Comparison

Compute for SPINK7-KLK5 and PSTI-chymotrypsinogen. Canonical protease-inhibitor interfaces should all show 55–70% nonpolar fraction.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Nonpolar fraction | 61% | Krowarsch 2003 |
| 95% CI | [55%, 67%] | benchmarks.md |
| Nonpolar BSA | ~930 Å² | 61% of 1530 Å² |
| Polar BSA | ~600 Å² | 39% of 1530 Å² |
| Surface average | ~57% | Protein surface typical |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Atom classification error** | Wrong atoms labeled polar/nonpolar | Backbone C and amide N misassigned | Medium |
| **Interface rearrangement** | Nonpolar fraction drifts over time | Conformational change | Low |
| **Probe radius sensitivity** | Different fractions with different probes | SASA algorithm parameter | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Crystal value | Nonpolar fraction from PDB 2PTC = ~61% | Matches literature |
| Atom classification check | Manual verification of C/N/O assignments | All atoms correctly classified |
| Total BSA consistency | Total BSA matches EXP-16 | Exact agreement |
| Time stability | Fraction CV < 5% | Stable composition |
| Surface enrichment | Interface > total surface nonpolar fraction | 61% > 57% |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
