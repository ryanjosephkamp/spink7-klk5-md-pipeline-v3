# EXP-20: Binding Loop Conservation Across Kazal Inhibitors

**Experiment ID:** EXP-20  
**Feature ID:** F-20 (benchmarks.md)  
**Category:** Structural  
**Status:** QUALITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates that the pipeline reproduces the structural conservation of the canonical binding loop across multiple Kazal-family inhibitors. Despite low sequence identity (20–40%), the binding loop (P3-P3') backbone geometry is remarkably conserved with RMSD < 1.0 Å across BPTI, SPINK1, PSTI, LEKTI D6, and SPINK7 (Laskowski & Qasim 2000, Krowarsch et al. 2003). The pipeline should reproduce this convergent evolution signature by showing that equilibrated structures of different inhibitors superpose with sub-angstrom loop RMSD.

---

## 2. Hypothesis

**H₁:** The pairwise backbone RMSD of the binding loop (P3-P3' residues) between any two Kazal-family inhibitors, after equilibration, will be < 1.0 Å.

**H₂:** The loop conservation will persist throughout production MD (time-averaged pairwise RMSD < 1.5 Å for all pairs).

---

## 3. Background and Rationale

The canonical serine protease inhibitor binding loop is one of the most conserved structural motifs in protein evolution. This conservation reflects the stringent geometric requirements for productive binding to the protease active site (Laskowski standard mechanism). Reproducing this conservation computationally validates that the force field and equilibration protocol preserve biologically meaningful structural features across different protein sequences.

---

## 4. Experimental Protocol

### 4.1 Systems

Compare binding loops from equilibrated structures of:
1. BPTI (PDB 2PTC / 4PTI) — P3-P3' = residues 13–20
2. SPINK1/PSTI (PDB 1TGS) — equivalent P3-P3' residues
3. SPINK7 (PDB 2LEO + homology) — equivalent P3-P3' residues
4. LEKTI D6 (homology model) — equivalent P3-P3' residues

### 4.2 Structural Alignment

1. Identify P3-P3' residues in each inhibitor by sequence-structure alignment to BPTI.
2. Superpose using only backbone atoms (N, Cα, C, O) of P3-P3'.
3. Compute pairwise RMSD for all 6 pairs (4 inhibitors → C(4,2) = 6 pairs).

### 4.3 MD-Based Comparison

From the last 50 ns of production MD for each system:
1. Extract 500 frames (100 ps intervals).
2. For each frame, align the binding loop to the BPTI reference (PDB 2PTC).
3. Compute time-averaged RMSD per system.
4. Report pairwise comparison across all four inhibitors.

### 4.4 Production MD Parameters (from config.py)

- Force field: AMBER ff14SB (`amber14-all.xml`)
- Duration: 100 ns per system
- Temperature: 310 K
- Save interval: 10 ps

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Pairwise Loop RMSD |
|---------------|-------------------|
| **PASS** | < 1.0 Å for all pairs (crystal); < 1.5 Å (MD average) |
| **MARGINAL** | < 1.5 Å for all pairs (crystal); < 2.0 Å (MD average) |
| **FAIL** | Any pair > 2.0 Å |

---

## 5. Control Conditions

### 5.1 Positive Control: Crystal-Crystal Superposition

Superpose the crystal structures directly (no MD) — this is the ground truth. BPTI vs. PSTI loop RMSD should be ~0.5 Å from crystallography.

### 5.2 Negative Control: Non-Kazal Protein Loop

Select a random surface loop of similar length from a non-inhibitor protein. Its RMSD to the BPTI binding loop should be > 2.0 Å, confirming the metric detects true conservation.

### 5.3 Disulfide Framework Control

Verify that the disulfide bonds flanking the binding loop (C14-C38 in BPTI equivalent) are intact in all systems. Disulfide integrity is prerequisite for loop geometry conservation.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| BPTI vs. PSTI loop RMSD | ~0.5 Å | Crystal structures |
| BPTI vs. SPINK7 loop RMSD | < 1.0 Å | Structural conservation |
| All-pairs max RMSD | < 1.0 Å (crystal) | Laskowski 2000 |
| MD average RMSD to reference | < 1.5 Å per system | Equilibrated dynamics |
| Loop φ/ψ conservation | All P1 in β-sheet region | Cross-ref EXP-14 |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Homology model loop bias** | SPINK7/LEKTI loops modeled from BPTI template | Circular validation | High |
| **Sequence alignment error** | Wrong P3-P3' residues identified | Manual alignment needed | Medium |
| **Loop flexibility** | One inhibitor's loop is significantly more flexible | Sequence-specific dynamics | Low |
| **Crystal packing effects** | Crystal vs. solution geometries differ | Environment difference | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Sequence alignment | P3-P3' mapped correctly in all 4 inhibitors | Consistent with UniProt annotations |
| Crystal superposition | BPTI-PSTI RMSD ≈ 0.5 Å | Matches published structures |
| Disulfide check | All conserved disulfides intact | S-S distances 2.0 ± 0.1 Å |
| Per-system RMSD from crystal | Each inhibitor's MD loop < 1.5 Å from own crystal | Structural stability |
| All-pairs table | 6 pairwise RMSDs all < 1.0 Å | Conservation confirmed |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
