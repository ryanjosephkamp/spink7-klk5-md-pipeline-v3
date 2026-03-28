# EXP-18: Interfacial Water Molecule Count

**Experiment ID:** EXP-18  
**Feature ID:** F-18 (benchmarks.md)  
**Category:** Structural  
**Status:** SEMI-QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the number of ordered water molecules at the protease-inhibitor interface. Crystal structures of BPTI-trypsin show ~15 crystallographic waters within 3.5 Å of both binding partners (Krowarsch et al. 2003). The 95% pipeline CI is [10, 20] water molecules. These interfacial waters mediate hydrogen bonds, fill cavities, and contribute to binding thermodynamics. Reproducing this count from explicit-solvent MD validates the solvation shell structure and the water model's behavior at the interface.

---

## 2. Hypothesis

**H₁:** The time-averaged number of water molecules within 3.5 Å of both the protease and inhibitor (interfacial waters) will be 15 ± 5 (95% CI: [10, 20]) in BPTI-trypsin.

**H₂:** A subset (~5–8) of these waters will be ordered (low B-factor equivalent, i.e., high occupancy and low positional RMSD) and correspond to crystallographic water sites.

---

## 3. Background and Rationale

Interfacial water molecules play a dual role: they mediate hydrogen bonds between the protease and inhibitor (bridging interactions) and contribute to the desolvation penalty upon binding. Crystal structures consistently show ~15 ordered waters at canonical protease-inhibitor interfaces. In explicit-solvent MD, these waters should be identifiable as high-density, low-mobility solvent molecules that persist at the interface. Accurately reproducing the interfacial water count validates the TIP3P water model and the force field's ability to describe protein-water interactions.

---

## 4. Experimental Protocol

### 4.1 System

BPTI-trypsin (PDB 2PTC) from EXP-04 production MD trajectory.

### 4.2 Interfacial Water Definition

A water molecule is "interfacial" if:
- At least one of its atoms is within 3.5 Å of a heavy atom on the protease (chain A)
- AND at least one of its atoms is within 3.5 Å of a heavy atom on the inhibitor (chain B)

### 4.3 Analysis Protocol

1. **Crystal structure count:** Identify crystallographic waters in PDB 2PTC satisfying the dual-contact criterion.
2. **Per-frame count:** For each frame (last 50 ns, 100 ps intervals = 500 frames):
   - Count all water molecules satisfying the interfacial criterion.
   - Record oxygen positions.
3. **Time-averaged count:** Mean ± std across all analyzed frames.
4. **Water density map:** Compute 3D water density on a 0.5 Å grid at the interface, averaged over all frames. High-density peaks (>2× bulk density) correspond to ordered water sites.
5. **Residence time analysis:** For each interfacial water site, compute mean residence time (how long individual water molecules remain at that site).

### 4.4 Production MD Parameters (from config.py)

- Water model: TIP3P (`amber14/tip3p.xml`)
- Force field: AMBER ff14SB
- Duration: 100 ns
- Temperature: 310 K
- Save interval: 10 ps

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Interfacial Water Count |
|---------------|------------------------|
| **PASS** | [10, 20] |
| **MARGINAL** | [7, 25] |
| **FAIL** | < 7 or > 25 |

---

## 5. Control Conditions

### 5.1 Crystal Water Comparison

Compare the 3D water density map from MD with crystallographic water positions in PDB 2PTC. High-density peaks should overlap with crystal water sites (distance < 1.5 Å).

### 5.2 Bulk Water Control

Count waters in a shell of equivalent volume but far from the interface (e.g., 15–18.5 Å from interface COM). This gives the bulk density reference (ρ_bulk), and the interfacial count should be distinguishable from random.

### 5.3 Free Protein Control

For free BPTI in solution, count waters within 3.5 Å of the binding loop surface. These should be more numerous than interfacial waters (exposed surface), confirming the interface excludes water upon binding.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Interfacial water count | ~15 | Krowarsch 2003 |
| 95% CI | [10, 20] | benchmarks.md |
| Crystal water sites matched | >60% of crystal waters reproduced | Convergent density |
| Ordered water residence time | >100 ps for top 5–8 sites | Low mobility |
| Bulk density ratio (interface/bulk) | >1.5 for ordered sites | Structured solvation |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **TIP3P over-diffusion** | Interfacial waters exchange too quickly | TIP3P diffusion ~2× too fast | Medium |
| **Cutoff sensitivity** | Count varies with 3.0 vs 3.5 vs 4.0 Å | Distance criterion arbitrary | Low |
| **Interface opening** | Water floods interface → high count | Complex partially dissociated | High |
| **Crystal packing artifacts** | Crystal waters not present in solution | Lattice contacts stabilize waters | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Crystal waters | Count in PDB 2PTC | ~15 (matches literature) |
| Frame-to-frame | Standard deviation < 5 | Stable count |
| Density map | Peaks at crystal water positions | Visual overlap |
| Bulk comparison | Interface count > random volume count | Signal above noise |
| Residence analysis | Some waters with τ > 100 ps | Ordered sites exist |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
