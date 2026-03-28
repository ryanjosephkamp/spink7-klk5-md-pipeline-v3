# EXP-16: Buried Surface Area (BSA) Upon Complex Formation

**Experiment ID:** EXP-16  
**Feature ID:** F-16 (benchmarks.md)  
**Category:** Structural  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the buried surface area (BSA) upon protease-inhibitor complex formation and validates it against the experimental benchmark of 1530 ± 170 Å² for BPTI-trypsin (Krowarsch et al. 2003, Janin & Chothia 1990). BSA is a fundamental structural descriptor of protein-protein interfaces, correlating with binding affinity (EXP-26). The pipeline computes BSA from SASA difference: BSA = SASA(free protease) + SASA(free inhibitor) − SASA(complex).

---

## 2. Hypothesis

**H₁:** The time-averaged BSA for BPTI-trypsin will be 1530 ± 170 Å² (95% CI from benchmarks.md), matching the crystal structure value.

**H₂:** BSA will be stable throughout 100 ns of production MD (CV < 15%), indicating a persistent, well-defined interface.

**H₃:** The BSA breakdown (polar vs. nonpolar) will show ~61% nonpolar contribution (cross-reference EXP-19).

---

## 3. Background and Rationale

BSA quantifies the extent of the protein-protein interface. For canonical serine protease-inhibitor pairs, BSA ≈ 1200–1800 Å², typical of tight, protease-substrate-like interfaces. BSA correlates with ΔG_bind (Horton & Lewis 1992), hydrophobic burial, and shape complementarity. Reproducing BSA validates the structural integrity of the simulated complex and provides a sanity check for binding free energy calculations.

---

## 4. Experimental Protocol

### 4.1 Systems

- **Primary:** BPTI-trypsin (PDB 2PTC, EXP-04 trajectory)
- **Secondary:** SPINK7-KLK5 (EXP-01), PSTI-chymotrypsinogen (EXP-05)

### 4.2 SASA Calculation Parameters

- Probe radius: 1.4 Å (standard water probe)
- Algorithm: Shrake-Rupley or LCPO
- Atom radii: AMBER ff14SB van der Waals radii

### 4.3 BSA Computation

For each frame from the last 50 ns of production MD (at 100 ps intervals = 500 frames):

1. **SASA(complex):** Total SASA of the protease-inhibitor complex.
2. **SASA(protease):** SASA of protease chain only (remove inhibitor, no re-equilibration).
3. **SASA(inhibitor):** SASA of inhibitor chain only (remove protease, no re-equilibration).
4. **BSA = SASA(protease) + SASA(inhibitor) − SASA(complex).**
5. Decompose into polar (N, O, S atoms) and nonpolar (C atoms) contributions.

### 4.4 Per-Residue BSA Decomposition

1. For each interface residue, compute per-residue BSA contribution.
2. Rank residues by BSA contribution.
3. Cross-reference: P1 should be the top contributor (EXP-07).

### 4.5 Production MD Parameters (from config.py)

- Force field: AMBER ff14SB
- Duration: 100 ns
- Temperature: 310 K
- Save interval: 10 ps

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | BSA (Å²) |
|---------------|----------|
| **PASS** | [1360, 1700] (1530 ± 170) |
| **MARGINAL** | [1100, 1950] |
| **FAIL** | Outside marginal range |

---

## 5. Control Conditions

### 5.1 Crystal Structure BSA

Compute BSA directly from PDB 2PTC (single frame, no MD). This provides the reference value and serves as a check that the SASA algorithm is implemented correctly.

### 5.2 Free Protein Control

For free BPTI and free trypsin simulations, SASA should match the values used in BSA calculation (self-consistency).

### 5.3 Random Contact Control

Create a non-interacting placement of BPTI near trypsin (touching but not docked). BSA should be < 200 Å² (minimal random contact), confirming BSA measures specific interface burial.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| BSA (BPTI-trypsin) | 1530 ± 170 Å² | Krowarsch 2003, Janin 1990 |
| BSA from crystal | ~1520 Å² | PDB 2PTC direct |
| Nonpolar fraction | ~61% | Cross-ref EXP-19 |
| BSA CV | < 15% | Stable interface |
| Top BSA residue | P1 (K15) | Cross-ref EXP-07 |
| SPINK7-KLK5 BSA | ~1200–1600 Å² | Kazal-type range |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Interface disruption** | BSA drops by >30% during MD | Equilibration failure | Critical |
| **Probe radius sensitivity** | Different from crystal BSA | SASA algorithm choice | Low |
| **Conformational sampling** | BSA fluctuates widely (CV > 20%) | Flexible loops | Medium |
| **Water penetration** | Crystal waters at interface inflate BSA | Explicit water treatment | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Crystal BSA | BSA(PDB 2PTC) = ~1520 Å² | Within 5% of literature |
| SASA algorithm | Compare Shrake-Rupley vs MDTraj | Agreement within 3% |
| Time stability | BSA time series: no monotonic drift | Slope < 50 Å²/100 ns |
| Decomposition check | Polar + nonpolar = total BSA | Exact sum |
| Per-residue sum | Σ per-residue BSA ≈ total BSA | Within 2% (overlap corrections) |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
