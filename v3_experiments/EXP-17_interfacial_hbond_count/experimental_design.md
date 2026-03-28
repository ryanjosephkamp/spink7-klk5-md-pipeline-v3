# EXP-17: Interfacial Hydrogen Bond Count

**Experiment ID:** EXP-17  
**Feature ID:** F-17 (benchmarks.md)  
**Category:** Structural  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the number of persistent intermolecular hydrogen bonds at the protease-inhibitor interface. For BPTI-trypsin, the crystal structure shows ~10 direct H-bonds, and the literature consensus is ~10 persistent interfacial H-bonds with a 95% CI of [7, 13] (Krowarsch et al. 2003, Radisky & Bhatt). The pipeline should reproduce this count from production MD trajectory analysis, providing essential structural validation of the binding interface electrostatics.

---

## 2. Hypothesis

**H₁:** The time-averaged number of persistent interfacial H-bonds in BPTI-trypsin will be 10 ± 3 (95% CI: [7, 13]).

**H₂:** The H-bond network will be stable (no net gain or loss > 2 bonds over 100 ns), indicating a well-maintained interface.

**H₃:** The total energetic contribution of H-bonds (count × ~1.5 kcal/mol from EXP-08) will be consistent with 20–35% of ΔG_bind.

---

## 3. Background and Rationale

Interfacial hydrogen bonds are critical determinants of protease-inhibitor binding specificity and affinity. In canonical complexes, a network of ~10 H-bonds forms between the binding loop (P3-P3') and the protease active-site cleft, including contacts to the oxyanion hole, the specificity pocket, and the enzyme backbone β-sheet. The number and quality of H-bonds reflect the complementarity of the binding interface and serve as a direct structural validator.

---

## 4. Experimental Protocol

### 4.1 System

BPTI-trypsin (PDB 2PTC) from EXP-04 production trajectory.

### 4.2 H-Bond Definition

Geometric criteria:
- Donor-acceptor distance: d(D···A) ≤ 3.5 Å
- D-H···A angle: θ ≥ 135°
- Only intermolecular H-bonds (donor on one chain, acceptor on the other)

### 4.3 Analysis Protocol

1. **Crystal structure count:** Count H-bonds in PDB 2PTC using the same criteria.
2. **Per-frame count:** For each of the 10,000 frames (100 ns at 10 ps intervals):
   - Count all intermolecular H-bonds satisfying criteria.
   - Record donor-acceptor pair identity.
3. **Persistent H-bonds:** H-bonds present in >50% of frames.
4. **H-bond occupancy table:** For each unique D-A pair, report occupancy (% frames present).
5. **Time series:** Plot total H-bond count vs. time.

### 4.4 Production MD Parameters (from config.py)

- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P
- Duration: 100 ns
- Temperature: 310 K, friction 1.0 ps⁻¹
- Timestep: 2 fs
- Save interval: 10 ps

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Persistent H-Bond Count |
|---------------|------------------------|
| **PASS** | [7, 13] |
| **MARGINAL** | [5, 16] |
| **FAIL** | < 5 or > 16 |

---

## 5. Control Conditions

### 5.1 Crystal Structure Control

Count H-bonds directly in PDB 2PTC. This serves as the reference and validates the geometric criteria implementation.

### 5.2 Intra-molecular H-Bond Control

Count intra-molecular H-bonds in BPTI and trypsin separately. These should be substantially more numerous than interfacial H-bonds, confirming the metric specifically captures the interface.

### 5.3 Water-Mediated H-Bonds

Count water-mediated bridges at the interface (D–H···O_water···H–A, both contacts ≤ 3.5 Å). Expected: ~15 water-mediated H-bonds (cross-reference EXP-18). These are separate from direct H-bonds.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Persistent H-bond count | ~10 | Krowarsch 2003 |
| 95% CI | [7, 13] | benchmarks.md |
| Crystal structure count | ~10 | PDB 2PTC |
| Key H-bonds | P1 backbone ↔ oxyanion hole, P3 ↔ enzyme β-sheet | Structural consensus |
| Stability | CV < 25% over 100 ns | Stable interface |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Geometric cutoff sensitivity** | Count varies significantly with d/θ thresholds | Arbitrary criteria | Medium |
| **Transient H-bonds overcounted** | Many H-bonds with <20% occupancy inflate count | Persistence threshold | Low |
| **Interface disruption** | Large fluctuations in count (CV > 40%) | Equilibration issue | High |
| **Water competition** | Some direct H-bonds replaced by water-mediated | Dynamic interface | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Crystal structure | Direct H-bond count ~10 | Matches literature |
| Cutoff sensitivity | Test d = 3.0, 3.5, 4.0 Å → counts should differ by < 30% | Robust result |
| Time stability | Count vs. time: no monotonic trend | Equilibrated |
| Key H-bonds present | P1→oxyanion and P3→β-sheet listed | Expected contacts found |
| Occupancy distribution | Most persistent H-bonds have >70% occupancy | Bimodal distribution |
| Cross-ref EXP-08 | Count × 1.5 kcal/mol ≈ total H-bond energy | Quantitative consistency |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
