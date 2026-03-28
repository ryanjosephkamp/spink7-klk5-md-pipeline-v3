# EXP-11: KLK5 Substrate Kinetics (Km/kcat) — Qualitative Binding Validation

**Experiment ID:** EXP-11  
**Feature ID:** F-11 (benchmarks.md)  
**Category:** Kinetic  
**Status:** QUALITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates that the pipeline's model of KLK5 produces substrate binding geometry consistent with enzyme activity. KLK5 cleaves Boc-VPR-AMC with Km = 0.20 mM and kcat = 196.8 min⁻¹ (Brattsand et al. 2005, Michael et al. 2005). Since direct enzyme kinetics are outside the scope of classical MD, this experiment uses substrate docking and short MD simulations to confirm that the substrate adopts a productive binding pose in the S1-S3 subsites, with correct distances to the catalytic triad (Ser195, His57, Asp102). This is a qualitative/structural validation supporting the pipeline's KLK5 model accuracy.

---

## 2. Hypothesis

**H₁:** Docking of the Boc-VPR-AMC substrate into the KLK5 active site will produce a pose where the scissile bond (Arg P1 C=O) is within 3.5 Å of the Ser195 Oγ and His57 Nε2.

**H₂:** Short (10 ns) MD simulation will maintain the productive binding geometry, with the P1 Arg side chain remaining in the S1 pocket.

---

## 3. Background and Rationale

Km and kcat are Michaelis-Menten parameters that characterize enzyme-substrate interactions. While computing these directly from MD is not feasible (requiring quantum-mechanical treatment of bond breaking/forming), validating that the substrate binds productively in the correct orientation is essential. If KLK5's active site geometry is correct in the pipeline's force field representation, the substrate should naturally dock with the P1 Arg in S1, the P2 Pro in S2, and the P3 Val in S3, with catalytic distances consistent with the Ser protease mechanism.

---

## 4. Experimental Protocol

### 4.1 KLK5 Model

Use the prepared KLK5 structure from EXP-01 (PDB 2PSX or homology model).

### 4.2 Substrate Preparation

1. Build Boc-VPR-AMC peptide substrate.
2. Generate GAFF/AM1-BCC parameters for non-standard groups (Boc, AMC).
3. Alternatively, use the tripeptide VPR as a simplified substrate (standard amino acids only).

### 4.3 Docking

1. Dock substrate into KLK5 active site using rigid-body placement guided by known serine protease substrate orientation.
2. P1 Arg → S1 specificity pocket (above Asp189).
3. Scissile bond carbonyl → oxyanion hole (Ser195/Gly193 backbone NH).

### 4.4 Short MD Simulation

- System: KLK5-substrate complex in explicit TIP3P water
- Force field: AMBER ff14SB (`amber14-all.xml`)
- Box padding: 1.2 nm, ionic strength: 0.15 M NaCl
- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- Equilibration: NVT 500 ps → NPT 1000 ps at 310 K
- Production: 10 ns (short — goal is geometry stability, not kinetics)
- Timestep: 2 fs
- Save interval: 10 ps

### 4.5 Analysis Metrics

1. **Catalytic distances (time-averaged):**
   - Ser195 Oγ → scissile C=O: should be < 3.5 Å
   - His57 Nε2 → Ser195 Oγ: should be < 3.5 Å  
   - His57 Nδ1 → Asp102 Oδ: should be < 3.0 Å

2. **P1 Arg position:** COM of Arg guanidinium relative to Asp189 carboxylate (should remain < 4.0 Å, indicating salt bridge).

3. **Substrate RMSD:** RMSD of substrate heavy atoms relative to initial docked pose.

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | Criterion |
|---------------|-----------|
| **PASS** | Substrate maintains productive geometry; all catalytic distances < 3.5 Å; P1 in S1 throughout |
| **MARGINAL** | Transient departures (< 20% of frames), but average geometry productive |
| **FAIL** | Substrate dissociates or P1 leaves S1 pocket |

---

## 5. Control Conditions

### 5.1 Positive Control

**BPTI-trypsin (EXP-04):** The P1 K15 maintains correct geometry in the S1 pocket throughout 100 ns. If the pipeline works for BPTI-trypsin, the KLK5 active site should similarly accommodate a P1 Arg/Lys.

### 5.2 Negative Control

**D-Arg substrate:** A D-amino acid at P1 should not fit productively in the S1 pocket. If docked, it should drift away during MD.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Km | 0.20 ± 0.05 mM | Brattsand 2005 (reference only) |
| kcat | 196.8 ± 25.9 min⁻¹ | Michael 2005 (reference only) |
| Ser195-Oγ to scissile C distance | < 3.5 Å | Structural consensus |
| P1 Arg in S1 | Maintained throughout | Trypsin-like specificity |
| Substrate RMSD | < 2.0 Å from initial pose | Stable binding |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Substrate parameterization error** | AMC/Boc parameters incorrect | Non-standard residue handling | Medium |
| **Docking artifact** | Wrong initial pose → MD unstable | No exhaustive search | Medium |
| **Active site distortion** | KLK5 loop rearrangement displaces substrate | Homology model quality | Medium |
| **Timescale limitation** | 10 ns too short to observe rebinding after transient departure | Kinetic accessibility | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| KLK5 model check | Catalytic triad geometry correct in starting structure | H57-D102 < 3.0 Å, S195-H57 < 3.5 Å |
| Docking check | Visual: P1 Arg in S1, correct chain orientation | Consistent with known structures |
| Minimization | No steric clashes (energy decrease > 1000 kJ/mol) | Converged minimization |
| Equilibration | RMSD of KLK5 backbone < 2.0 Å from start | Stable protein |
| Catalytic distances stable | Time series of key distances | All < 3.5 Å for >80% of frames |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
