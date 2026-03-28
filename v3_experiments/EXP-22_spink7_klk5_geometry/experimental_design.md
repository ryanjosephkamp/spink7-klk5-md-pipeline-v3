# EXP-22: SPINK7–KLK5 Interface Geometry Validation

**Experiment ID:** EXP-22  
**Feature ID:** F-22 (benchmarks.md)  
**Category:** Structural  
**Status:** QUALITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates the overall structural geometry of the SPINK7-KLK5 complex, the primary target of the pipeline. Since no experimental co-crystal structure exists for SPINK7-KLK5, the complex is modeled by docking (ClusPro or HADDOCK) followed by MD refinement. This experiment validates the docked model by checking: (1) canonical binding mode geometry, (2) consistency with known Kazal-protease interaction patterns, (3) structural stability during MD, and (4) internal consistency with other experimental benchmarks.

---

## 2. Hypothesis

**H₁:** The SPINK7-KLK5 docked complex will adopt the canonical Laskowski binding mode with the SPINK7 reactive-site loop inserted in the KLK5 active-site cleft.

**H₂:** The complex will be structurally stable during 100 ns MD (backbone RMSD < 3.0 Å from the equilibrated structure).

**H₃:** The interface geometry will satisfy all structural benchmarks: BSA (~1200–1600 Å²), H-bond count (7–13), catalytic distance (~2.7 Å), and P1 in S1.

---

## 3. Background and Rationale

SPINK7-KLK5 is the central system of this pipeline. All thermodynamic predictions (EXP-01, EXP-02) depend on the quality of the structural model. Since no experimental structure is available, careful validation against multiple structural criteria is essential. The complex should adhere to the canonical protease-inhibitor binding geometry observed in all well-characterized Kazal-protease pairs.

---

## 4. Experimental Protocol

### 4.1 Complex Construction

1. **SPINK7 model:** PDB 2LEO (NMR, 20 models, select model 1 or best DOPE score).
2. **KLK5 model:** PDB 2PSX (X-ray) or homology model.
3. **Docking:** ClusPro protein-protein docking, constrain SPINK7 binding loop toward KLK5 active site.
4. **Pose selection:** Choose pose with P1 Arg in S1 pocket, canonical orientation.

### 4.2 MD Refinement and Validation

Using EXP-01 production trajectory:

**4.2.1 Global Stability:**
- Backbone RMSD vs. equilibrated structure: time series over 100 ns.
- Radius of gyration: should remain stable.

**4.2.2 Interface-Specific Metrics:**
- BSA: compute per-frame, compare to EXP-16 range (1200–1600 Å²).
- H-bond count: compare to EXP-17 range (7–13).
- Catalytic distance: compare to EXP-15 (2.4–3.0 Å).
- P1 in S1: salt bridge d(Asp, P1 Arg) < 4.0 Å in >80% frames.

**4.2.3 Binding Mode Validation:**
- Superpose SPINK7-KLK5 onto BPTI-trypsin (2PTC): the binding loops should align.
- Compute interface RMSD after superposition: should be < 2.0 Å.

### 4.3 Simulation Parameters (from config.py)

- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P (`amber14/tip3p.xml`)
- Box padding: 1.2 nm, ionic strength 0.15 M
- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT: 500 ps, NPT: 1000 ps at 310 K
- Production: 100 ns, timestep 2 fs, save 10 ps

### 4.4 Acceptance Criteria (from benchmarks.md)

| Classification | Criterion |
|---------------|-----------|
| **PASS** | Canonical binding mode; backbone RMSD < 3.0 Å; all structural metrics in PASS range |
| **MARGINAL** | Binding mode correct; RMSD < 4.0 Å; most metrics in MARGINAL range |
| **FAIL** | Non-canonical binding mode; interface disruption; RMSD > 4.0 Å |

---

## 5. Control Conditions

### 5.1 Positive Control: BPTI-Trypsin (EXP-04)

The gold-standard complex with a crystal structure. All interface metrics from EXP-04 should fall well within PASS ranges, validating the analysis pipeline.

### 5.2 Docking Redundancy Control

Run 3 independent docking calculations (different starting orientations or different docking servers). All should yield the same binding mode. If not, the model is ambiguous.

### 5.3 SPINK7 Stability Control

Run free SPINK7 in solution (no KLK5) for 100 ns. Confirm the binding loop remains in canonical geometry (cross-ref EXP-14) — verifying that SPINK7 independently maintains its binding-competent conformation.

---

## 6. Expected Outcomes

| Metric | Expected Value | Cross-Reference |
|--------|---------------|-----------------|
| Backbone RMSD (100 ns) | < 3.0 Å | Structural stability |
| BSA | 1200–1600 Å² | EXP-16 |
| H-bond count | 7–13 | EXP-17 |
| Catalytic distance | 2.4–3.0 Å | EXP-15 |
| P1 in S1 | >80% occupancy of salt bridge | EXP-21 |
| Binding mode alignment to 2PTC | RMSD < 2.0 Å | Canonical mode |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Wrong docking mode** | P1 not in S1, non-canonical orientation | Docking scoring function | Critical |
| **SPINK7 model quality** | NMR ensemble model has loop artifacts | 2LEO resolution/ensemble | High |
| **KLK5 model quality** | Active site distorted in homology model | Template-target differences | High |
| **Interface instability** | BSA drops >30% during MD | Weak predicted interactions | Medium |
| **No co-crystal validation** | Cannot directly verify binding mode | Fundamental data gap | Acknowledged limitation |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| SPINK7 model selection | Best model from 2LEO NMR ensemble | Lowest energy / best stereochemistry |
| KLK5 model validation | Catalytic triad intact, S1 pocket formed | Triad distances < 3.5 Å |
| Docking result | P1 in S1, canonical orientation | Visual + salt bridge check |
| Minimization | No steric clashes at interface | Energy decrease > 1000 kJ/mol |
| Equilibration | Backbone RMSD plateaus < 3.0 Å | Converged |
| 10 ns checkpoint | All interface metrics in expected range | Early validation |
| 50 ns checkpoint | Stable time series, no drift | Mid-run check |
| Final 100 ns | All metrics within PASS range | Complete validation |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
