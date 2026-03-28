# EXP-27: Interface Packing Density (V/Vo Ratio)

**Experiment ID:** EXP-27  
**Feature ID:** F-27 (benchmarks.md)  
**Category:** Biophysical  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the packing density at the protease-inhibitor interface, measured as the volume ratio V/Vo, where V is the actual volume occupied by interface atoms and Vo is the Voronoi volume. For well-packed protein-protein interfaces, V/Vo ≈ 1.00 (Krowarsch et al. 2003), indicating near-crystalline packing density comparable to the protein interior. The 95% CI is [0.95, 1.05]. This validates that the simulated interface has no pathological voids or clashes.

---

## 2. Hypothesis

**H₁:** The interface packing density V/Vo for BPTI-trypsin will be 1.00 ± 0.05 (95% CI: [0.95, 1.05]).

**H₂:** The V/Vo ratio will be comparable to the protein core packing density (~0.74–0.78 for vdW volume fraction), confirming the interface is as well-packed as the interior.

---

## 3. Background and Rationale

The packing density at protein-protein interfaces is a key indicator of binding quality. Well-evolved interfaces (like protease-inhibitor pairs) are tightly packed, with complementary shapes that exclude voids and maximize van der Waals contacts. A V/Vo ≈ 1.0 indicates optimal packing, similar to close-packed spheres in the protein interior. Values significantly <1.0 suggest voids (poor complementarity), while values >1.0 indicate steric strain (atomic overlaps in the equilibrated structure suggest force field issues).

---

## 4. Experimental Protocol

### 4.1 System

BPTI-trypsin (PDB 2PTC) from EXP-04 production trajectory.

### 4.2 Volume Calculation Method

**4.2.1 Voronoi Tessellation:**
1. Select interface atoms (within 5.0 Å of the partner chain).
2. Compute Voronoi tessellation of the interface region.
3. Vo = sum of Voronoi cell volumes for interface atoms.

**4.2.2 van der Waals Volume:**
1. V = sum of atomic vdW volumes (4/3 πr³ for each atom, using AMBER radii).
2. Correct for atomic overlap in bonded atoms.

**4.2.3 V/Vo Ratio:**
V/Vo = total vdW volume / total Voronoi volume for interface atoms.

### 4.3 Alternative Method: Cavity Detection

1. Roll a 1.4 Å probe over the interface.
2. Count unoccupied volume (cavities).
3. Packing efficiency = 1 − (cavity volume / total interface volume).
4. For well-packed interface, packing efficiency should be ~0.74 (close-packed fraction).

### 4.4 Time-Averaging

From last 50 ns of production MD (100 ps intervals = 500 frames):
1. Compute V/Vo per frame.
2. Report mean ± std.
3. Block averaging: 5 blocks of 10 ns each.

### 4.5 Production MD Parameters (from config.py)

- Force field: AMBER ff14SB
- Duration: 100 ns
- Temperature: 310 K
- Save interval: 10 ps

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | V/Vo |
|---------------|------|
| **PASS** | [0.95, 1.05] |
| **MARGINAL** | [0.90, 1.10] |
| **FAIL** | Outside marginal range |

---

## 5. Control Conditions

### 5.1 Protein Interior Control

Compute V/Vo for the BPTI core (interior atoms >5 Å from surface). Expected: ~0.74 (close-packed fraction for protein interiors). This validates the Voronoi tessellation implementation.

### 5.2 Crystal Structure Reference

Compute V/Vo directly from PDB 2PTC (no MD). Should match literature value of ~1.00.

### 5.3 Artificial Gap Control

Shift BPTI 0.5 Å away from trypsin and compute V/Vo. Expected V/Vo < 0.90 (voids at interface), confirming the metric detects poor packing.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| V/Vo (BPTI-trypsin interface) | 1.00 ± 0.05 | Krowarsch 2003 |
| V/Vo (protein interior) | ~0.74 | Close-packed reference |
| Time stability | CV < 5% | Stable packing |
| Crystal reference | ~1.00 | PDB 2PTC |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Voronoi boundary effects** | Artifacts at interface edges | Tessellation of mixed regions | Medium |
| **Water penetration** | Interface waters reduce V/Vo | Dynamic solvation | Low |
| **Force field steric errors** | V/Vo > 1.1 (overlap) or < 0.9 (voids) | vdW radius assignment | Medium |
| **Computational method choice** | Voronoi vs. cavity detection disagree | Implementation-dependent | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Interface atom selection | ~30–60 residues identified as interfacial | Reasonable count |
| Voronoi quality | No degenerate cells or infinite volumes | Algorithm converged |
| Crystal V/Vo | ~1.00 from PDB 2PTC | Matches literature |
| Interior V/Vo | ~0.74 for protein core | Known reference |
| Time stability | V/Vo CV < 5% over 50 ns | Well-sampled |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
