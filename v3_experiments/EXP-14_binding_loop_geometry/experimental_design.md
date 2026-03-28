# EXP-14: Canonical Binding Loop Geometry

**Experiment ID:** EXP-14  
**Feature ID:** F-14 (benchmarks.md)  
**Category:** Structural  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates that the pipeline preserves the canonical binding loop geometry of Kazal-type inhibitors during MD simulation. The binding (reactive-site) loop adopts a conserved substrate-like conformation that inserts into the protease active site. Key metrics include the backbone RMSD of the P3-P3' loop residues relative to crystal structures, the Ramachandran angles of the P1 residue, and the maintenance of the canonical disulfide-stabilized loop architecture. BPTI-trypsin (PDB 2PTC) and SPINK1/PSTI structures serve as structural references.

---

## 2. Hypothesis

**H₁:** The binding loop (P3-P3' residues) backbone RMSD relative to the starting crystal structure will remain < 1.5 Å throughout 100 ns of production MD.

**H₂:** The P1 residue backbone dihedral angles (φ, ψ) will remain in the extended β-sheet region (φ ≈ −120°, ψ ≈ +130°), consistent with substrate-like conformation.

**H₃:** The loop backbone hydrogen bond pattern (P1 C=O → oxyanion hole, P3 N-H → enzyme backbone) will persist with >80% occupancy.

---

## 3. Background and Rationale

The Laskowski standard mechanism requires the inhibitor binding loop to adopt a substrate-like conformation, presenting the scissile bond to the catalytic Ser. This geometry is remarkably conserved across all canonical serine protease inhibitors despite low sequence identity. Maintaining this geometry in MD is essential — if the loop distorts, all downstream binding free energy calculations are compromised. This experiment tests the most fundamental structural prerequisite.

---

## 4. Experimental Protocol

### 4.1 Systems

- **Primary:** BPTI-trypsin (PDB 2PTC) from EXP-04 production trajectory
- **Secondary:** SPINK7-KLK5 (EXP-01) and SPINK1-trypsin (EXP-06) trajectories

### 4.2 Production MD Parameters

From config.py:
- Duration: 100 ns
- Temperature: 310 K (Langevin, friction 1.0 ps⁻¹)
- Timestep: 2 fs
- Save interval: 10 ps (10,000 frames)
- Checkpoint: 100 ps

### 4.3 Analysis Protocol

**4.3.1 Binding Loop RMSD:**
1. Define binding loop as P3-P3' residues (e.g., BPTI: residues 13–20).
2. Align trajectory to protein backbone (global fit).
3. Compute backbone RMSD of loop residues for each frame vs. crystal structure.
4. Report: time series, mean, standard deviation, maximum deviation.

**4.3.2 P1 Backbone Dihedrals:**
1. Compute φ and ψ for P1 residue (K15 in BPTI) across all frames.
2. Plot Ramachandran distribution.
3. Check: >90% of frames in extended β region (φ ∈ [−180°, −60°], ψ ∈ [+60°, +180°]).

**4.3.3 Key H-Bond Occupancy:**
1. P1 C=O → oxyanion hole (Ser195/Gly193 backbone NH): persistent H-bond.
2. P3 backbone → enzyme backbone: antiparallel β-sheet pairing.
3. Report occupancy (% frames with D-A < 3.5 Å, angle > 135°).

**4.3.4 Disulfide Bond Integrity:**
1. Monitor all disulfide S-S distances (BPTI: C14-C38, C5-C55, C30-C51).
2. All must remain within 1.8–2.2 Å (covalent bond range).

### 4.4 Acceptance Criteria (from benchmarks.md)

| Classification | Loop RMSD | P1 φ/ψ in β |
|---------------|-----------|-------------|
| **PASS** | < 1.5 Å (mean) | > 90% frames |
| **MARGINAL** | 1.5–2.5 Å | 70–90% frames |
| **FAIL** | > 2.5 Å or loop distortion | < 70% frames |

---

## 5. Control Conditions

### 5.1 Positive Control

**Free BPTI (unbound):** The binding loop should also maintain its geometry in the free inhibitor (known from NMR/crystal), demonstrating intrinsic stability.

### 5.2 Structural Reference

**PDB superposition:** Overlay BPTI (4PTI), PSTI (1TGS), SPINK1 crystal structures — all share the canonical loop. Simulation should preserve this convergent geometry.

### 5.3 Negative Control

**High-temperature run (500 K):** At denaturing temperatures, the binding loop should lose its canonical geometry, confirming the metric is sensitive to structural disruption.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Loop RMSD (BPTI) | < 1.0 Å | Crystal 2PTC reference |
| P1 φ | −120 ± 30° | β-sheet region |
| P1 ψ | +130 ± 30° | β-sheet region |
| Oxyanion H-bond occupancy | > 85% | Canonical mechanism |
| Disulfide S-S distances | 2.03 ± 0.05 Å | Covalent bond |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Loop flexibility overestimated** | RMSD > 2 Å, large fluctuations | Force field flexibility | High |
| **P1 dihedral flip** | P1 exits β-sheet region | Rare event sampling | Medium |
| **Disulfide strain** | S-S distance fluctuates > 2.5 Å | Force field bonded terms | Low |
| **Homology model loop error** | SPINK7 loop geometry wrong from start | No crystal structure for free SPINK7 | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Starting structure | Loop RMSD = 0.0 Å vs crystal | Correct reference |
| First 1 ns | RMSD < 1.0 Å | Early stability |
| Disulfide check (all frames) | All S-S in [1.8, 2.2] Å | Bond integrity |
| Visual inspection | Superpose frame 0 / 50 ns / 100 ns | No visible distortion |
| Cross-system | All 3 inhibitors maintain canonical loop | Robustness |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
