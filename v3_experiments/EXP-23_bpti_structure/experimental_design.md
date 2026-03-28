# EXP-23: BPTI Crystal Structure Reproduction

**Experiment ID:** EXP-23  
**Feature ID:** F-23 (benchmarks.md)  
**Category:** Structural  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates that the pipeline's MD simulations of BPTI faithfully reproduce the crystal structure geometry. BPTI (PDB 4PTI / 5PTI) is one of the best-characterized small proteins, with ultra-high-resolution crystal structures (1.0 Å) and extensive NMR data. The backbone RMSD of the equilibrated MD structure relative to the crystal should be < 1.5 Å, and key structural features (disulfide bonds, secondary structure, hydrogen bonds) should be preserved. This is the most fundamental structural validation for the force field and simulation protocol.

---

## 2. Hypothesis

**H₁:** The time-averaged backbone RMSD of BPTI relative to the crystal structure (PDB 4PTI) will be < 1.5 Å throughout 100 ns of production MD.

**H₂:** All three disulfide bonds (C5-C55, C14-C38, C30-C51) will maintain S-S distances within 2.0 ± 0.15 Å.

**H₃:** Secondary structure elements (α-helix, β-sheet) will persist with >90% occupancy via DSSP analysis.

---

## 3. Background and Rationale

BPTI is the "hydrogen atom of protein folding" — a 58-residue protein that has been extensively studied by X-ray, NMR, and MD simulations for decades. Reproducing its crystal structure in MD validates the entire simulation workflow (force field, water model, equilibration protocol, thermostat/barostat). Any deviation > 2 Å would indicate systematic problems with the pipeline.

---

## 4. Experimental Protocol

### 4.1 System Preparation

- PDB: 4PTI (1.0 Å resolution, free BPTI)
- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P (`amber14/tip3p.xml`)
- Box padding: 1.2 nm, cubic, ionic strength 0.15 M NaCl
- pH: 7.4

### 4.2 Equilibration (from config.py)

- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT: 500 ps at 310 K, friction 1.0 ps⁻¹, restraints 1000 kJ/mol/nm² on heavy atoms
- NPT: 1000 ps at 310 K, 1.0 atm, barostat interval 25 steps
- Timestep: 2 fs

### 4.3 Production MD

- Duration: 100 ns
- Temperature: 310 K (Langevin, friction 1.0 ps⁻¹)
- Pressure: 1.0 atm (Monte Carlo barostat, interval 25)
- Timestep: 2 fs
- Save interval: 10 ps (10,000 frames)
- Checkpoint: 100 ps

### 4.4 Structural Analysis

**4.4.1 Backbone RMSD:**
1. Align each frame to crystal structure (4PTI).
2. Compute backbone (N, Cα, C, O) RMSD.
3. Report: time series, mean, std, max.

**4.4.2 Per-Residue RMSF:**
1. Compute Cα RMSF per residue from MD.
2. Compare to crystallographic B-factors: RMSF_xtal = √(3B/(8π²)).
3. Correlation coefficient should be > 0.5.

**4.4.3 Disulfide Bonds:**
1. Monitor S-S distances for C5-C55, C14-C38, C30-C51.
2. All should remain in [1.85, 2.15] Å.

**4.4.4 Secondary Structure:**
1. DSSP assignment per frame.
2. α-helix (residues 47-56) and β-sheet (residues 18-24, 29-35) should persist >90% of frames.

**4.4.5 Hydrogen Bond Network:**
1. Count intramolecular H-bonds per frame.
2. Key backbone H-bonds should match crystal structure pattern.

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Backbone RMSD |
|---------------|---------------|
| **PASS** | < 1.5 Å (mean) |
| **MARGINAL** | 1.5–2.5 Å |
| **FAIL** | > 2.5 Å |

---

## 5. Control Conditions

### 5.1 NMR Comparison

Compare RMSF profile to NMR order parameters or chemical shift data (published extensively for BPTI). Flexible regions (loop 8-12, termini) should show higher RMSF.

### 5.2 Temperature Control

Run identical simulation at 300 K (closer to crystal conditions). RMSD should be slightly lower than the 310 K run.

### 5.3 Water Structure

Verify that crystallographic water sites (well-characterized in 4PTI at 1.0 Å) are reproduced as high-density peaks in the MD water density map (cross-ref EXP-18 methodology).

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Backbone RMSD | 0.8–1.2 Å | Literature MD studies |
| Disulfide S-S distances | 2.03 ± 0.05 Å | Covalent bond |
| α-helix occupancy | >95% | Stable secondary structure |
| β-sheet occupancy | >90% | Stable secondary structure |
| B-factor correlation | r > 0.5 | RMSF-B agreement |
| RMSF (loop 8-12) | > 1.0 Å | Known flexible region |
| RMSF (core) | < 0.5 Å | Rigid core |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **RMSD drift** | Monotonic increase beyond 2 Å | Force field systematic error | Critical |
| **Disulfide distortion** | S-S distance > 2.3 Å | Bond parameter error | High |
| **Secondary structure loss** | α-helix or β-sheet breaks | Equilibration insufficient | Medium |
| **Crystal packing absence** | Some crystal contacts stabilize structure | In-solution differs | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Minimized structure | RMSD < 0.3 Å from crystal | Minimal perturbation |
| Post-NVT | RMSD < 1.0 Å | Stable with restraints |
| Post-NPT | RMSD < 1.5 Å | Stable without restraints |
| 10 ns checkpoint | RMSD equilibrated, no drift | Converged |
| Disulfide distances | All three in [1.85, 2.15] Å at all times | Bond integrity |
| DSSP | Secondary structure present throughout | No unfolding |
| Final 100 ns | Mean RMSD < 1.5 Å | PASS criterion |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
