# EXP-15: Catalytic Contact Distance (Ser195–P1 Carbonyl)

**Experiment ID:** EXP-15  
**Feature ID:** F-15 (benchmarks.md)  
**Category:** Structural  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates that the pipeline's MD simulations maintain the critical catalytic contact distance between the Ser195 Oγ and the P1 residue scissile bond carbonyl carbon. In crystal structures of canonical serine protease-inhibitor complexes, this distance is ~2.7 Å (Radisky & Bhatt, Vincent & Bhatt 2007), consistent with the tetrahedral intermediate geometry. The time-averaged distance from production MD must reproduce this with 95% CI of [2.4, 3.0] Å.

---

## 2. Hypothesis

**H₁:** The time-averaged distance between Ser195 Oγ and the P1 carbonyl carbon will be 2.7 ± 0.3 Å (95% CI: [2.4, 3.0] Å) in the BPTI-trypsin complex.

**H₂:** This distance will be maintained with < 10% of frames showing distances > 3.5 Å (indicating occasional but not persistent loss of catalytic contact).

---

## 3. Background and Rationale

The Ser195 Oγ → P1 carbonyl carbon distance is the most fundamental structural parameter in serine protease catalysis. In the Michaelis complex (substrate-like binding), this distance is ~2.7 Å, poised for nucleophilic attack. In inhibitor complexes, this distance is maintained but the reaction does not proceed (inhibitors are poor substrates). Reproducing this distance validates the active-site geometry, force field accuracy, and the physical realism of the inhibitor binding mode.

---

## 4. Experimental Protocol

### 4.1 Systems

- **Primary:** BPTI-trypsin (PDB 2PTC) — K15 Cα carbonyl to Ser195 Oγ
- **Secondary:** PSTI-chymotrypsinogen (PDB 1TGS), SPINK7-KLK5 (EXP-01)

### 4.2 Crystal Structure Reference

Measure Ser195 Oγ → P1 C distance in:
- 2PTC: K15 C → S195 Oγ (expected: ~2.7 Å)
- 1TGS: P1 C → S195 Oγ

### 4.3 MD Trajectory Analysis

Production MD parameters (from config.py):
- Duration: 100 ns
- Temperature: 310 K
- Timestep: 2 fs
- Save interval: 10 ps (10,000 frames)

Analysis:
1. For each frame, compute d(Ser195 Oγ, P1 C).
2. Report: time series, running average, histogram, mean ± std.
3. Block averaging: divide trajectory into 10 blocks of 10 ns each, compute mean per block, report standard error.

### 4.4 Additional Catalytic Triad Distances

Monitor the complete catalytic triad geometry:
1. d(Ser195 Oγ, His57 Nε2): expected ~2.8 Å
2. d(His57 Nδ1, Asp102 Oδ): expected ~2.8 Å
3. All three distances must remain < 3.5 Å simultaneously for catalytic competence.

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | d(Ser195 Oγ, P1 C) |
|---------------|---------------------|
| **PASS** | Mean in [2.4, 3.0] Å |
| **MARGINAL** | Mean in [2.0, 3.5] Å |
| **FAIL** | Mean outside [2.0, 3.5] Å or > 20% frames > 4.0 Å |

---

## 5. Control Conditions

### 5.1 Crystal Structure Control

The measured distance in the crystal structure (2PTC) provides the ground truth. The MD average should converge to within 0.3 Å of this reference.

### 5.2 Free Protease Control

In the unbound KLK5/trypsin, the Ser195 Oγ should be solvent-accessible (no contact partner). This confirms the distance metric is specific to the bound complex.

### 5.3 Cross-System Consistency

All three systems (BPTI-trypsin, PSTI-chymotrypsinogen, SPINK7-KLK5) should show the same catalytic distance distribution, confirming it's a universal feature of canonical inhibitor complexes.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| d(Ser195 Oγ, P1 C) | 2.7 ± 0.3 Å | Radisky & Bhatt, Vincent & Bhatt |
| 95% CI | [2.4, 3.0] Å | benchmarks.md |
| d(S195 Oγ, H57 Nε2) | ~2.8 Å | Catalytic triad |
| d(H57 Nδ1, D102 Oδ) | ~2.8 Å | Catalytic triad |
| Frames with d > 3.5 Å | < 10% | Rare fluctuations only |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Catalytic contact loss** | P1 drift from active site | Loop instability (see EXP-14) | Critical |
| **His57 flip** | Imidazole ring rotates, breaking triad | Chi2 dihedral sampling | Medium |
| **Force field distance bias** | Systematic shift from crystal | AMBER ff14SB vdW parameters | Low |
| **Protonation state error** | His57 protonated on wrong nitrogen | pKa prediction | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Crystal reference | d = 2.7 Å in 2PTC | Matches literature |
| First frame | d matches minimized structure | No large perturbation |
| Time convergence | Running average stabilizes within 20 ns | Equilibrated |
| Block averaging | Standard error < 0.2 Å | Well-sampled |
| Triad integrity | All 3 distances < 3.5 Å simultaneously in >80% frames | Catalytically competent |
| Cross-system | Same distance in BPTI, PSTI, SPINK7 complexes | Universal feature |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
