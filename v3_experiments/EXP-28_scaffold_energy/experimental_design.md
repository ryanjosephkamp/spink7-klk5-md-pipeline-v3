# EXP-28: Scaffold Energy Contribution to Binding

**Experiment ID:** EXP-28  
**Feature ID:** F-28 (benchmarks.md)  
**Category:** Biophysical  
**Status:** SEMI-QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the energetic contribution of the inhibitor scaffold (the non-contact residues that maintain the binding loop geometry) to the overall binding free energy. The scaffold contributes approximately −7.9 kcal/mol (Krowarsch et al. 2003), representing the thermodynamic cost of pre-organizing the binding loop in the correct conformation. This experiment uses free-energy perturbation or energy decomposition to separate scaffold from contact contributions, validating the pipeline's ability to dissect the energetic architecture of protein-protein binding.

---

## 2. Hypothesis

**H₁:** Energy decomposition will identify a scaffold contribution of −7.9 ± 3.0 kcal/mol to the total binding free energy.

**H₂:** The scaffold contribution will be primarily entropic (pre-organizing the loop reduces configurational entropy loss upon binding, which is energetically favorable from a binding perspective).

---

## 3. Background and Rationale

In Kazal-type inhibitors, the binding loop (P3-P3') makes direct contacts with the protease, but the surrounding scaffold (disulfide bonds, hydrophobic core, secondary structure) pre-organizes the loop in the correct conformation. Without the scaffold, the free loop peptide binds much more weakly (ΔΔG ≈ 7–8 kcal/mol). This "scaffold energy" quantifies the benefit of structural pre-organization and is a hallmark of evolved protein-protein interfaces.

---

## 4. Experimental Protocol

### 4.1 Approach: Free Loop Peptide Comparison

**Method A: Loop peptide vs. intact inhibitor**

1. **Full inhibitor binding:** ΔG_bind from EXP-04 (BPTI-trypsin, expected: −18.0 kcal/mol).
2. **Free loop peptide:** Extract the P3-P3' peptide (residues 13-20 of BPTI) as an isolated 8-residue peptide.
3. **Peptide binding:** Compute ΔG_bind for the isolated loop peptide to trypsin using US/WHAM.
4. **Scaffold energy:** ΔG_scaffold = ΔG_bind(full) − ΔG_bind(peptide).

### 4.2 Free Loop Peptide Simulation

- Build 8-residue peptide (BPTI residues 13-20, no disulfide bonds).
- Cap with ACE/NME.
- System preparation: AMBER ff14SB, TIP3P, box padding 1.2 nm, 0.15 M NaCl.
- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm.
- Equilibration: NVT 500 ps, NPT 1000 ps at 310 K.
- Production: 50 ns to characterize free peptide conformational ensemble.

### 4.3 Peptide-Trypsin US/WHAM

- Dock peptide into trypsin active site (same pose as in BPTI complex).
- US: ξ range 1.5–4.0 nm, spacing 0.05 nm (51 windows), k = 1000 kJ/mol/nm².
- Per-window: 10 ns (after 200 ps equilibration).
- WHAM: tolerance 10⁻⁶, max_iterations 100,000, bins 200.
- Bootstrap: n_bootstrap = 200.

### 4.4 Method B: Energy Decomposition

From EXP-04 trajectory:
1. Partition inhibitor residues into "contact" (P3-P3') and "scaffold" (all others).
2. Compute interaction energy of scaffold with protease (should be small — scaffold doesn't directly contact protease).
3. Compute conformational energy of binding loop in bound vs. free inhibitor (strain energy).
4. ΔG_scaffold ≈ ΔG_conformational(free→bound loop) — the cost avoided by pre-organization.

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Scaffold Energy |
|---------------|----------------|
| **PASS** | [−4.9, −10.9] kcal/mol (−7.9 ± 3.0) |
| **MARGINAL** | [−2.0, −14.0] kcal/mol |
| **FAIL** | Outside marginal range or wrong sign |

---

## 5. Control Conditions

### 5.1 Free Peptide Conformational Check

The isolated loop peptide should sample multiple conformations (no pre-organization). Its Ramachandran distribution should be much broader than the intact BPTI loop — confirming the scaffold reduces conformational entropy.

### 5.2 Scaffold-Protease Interaction Control

The scaffold (non-contact residues) should have near-zero direct interaction energy with the protease. If significant, the scaffold contributes directly rather than through pre-organization.

### 5.3 Disulfide Contribution

Removing C14-C38 disulfide (computationally, e.g., by mutation to Ala-Ala) should reduce ΔG_bind by a portion of the scaffold energy. Cross-reference EXP-31 (disulfide ablation ΔΔG ≈ 7 kcal/mol).

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Scaffold energy | −7.9 kcal/mol | Krowarsch 2003 |
| 95% CI | [−4.9, −10.9] | benchmarks.md |
| Loop peptide ΔG_bind | ~−10 kcal/mol | Full ΔG minus scaffold |
| Full BPTI ΔG_bind | ~−18.0 kcal/mol | EXP-04 |
| Free peptide conformational entropy | Much larger than bound loop | Pre-organization validated |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Loop peptide instability** | Peptide unfolds before binding | No scaffold support | High |
| **Docking bias** | Peptide placed in same orientation as full BPTI | Initial condition dependence | Medium |
| **Sampling insufficiency** | Free peptide doesn't explore enough conformations | 50 ns may be short | Medium |
| **Decomposition ambiguity** | Scaffold/contact boundary arbitrary | Method B sensitivity | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Peptide construction | 8-residue peptide matches BPTI 13-20 | Correct sequence |
| Free peptide RMSD | Large fluctuations (>3 Å) — confirms no pre-organization | Flexible peptide |
| Peptide-protease complex | Peptide stays in active site during US | Successful sampling |
| Method consistency | Methods A and B agree within 3 kcal/mol | Robustness |
| EXP-31 cross-check | Disulfide ablation ΔΔG ≈ scaffold energy | Consistent picture |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~3–5 GB (loop peptide–trypsin complex: ~30,000 atoms solvated; smaller than full BPTI-trypsin) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation (loop peptide extraction) | 0.2 | 0.2 | Extract binding loop from BPTI; cap termini |
| Equilibration (peptide-trypsin complex) | 0.5 | 0.3 | 1.5 ns total |
| Umbrella Sampling (51 windows × 10.2 ns) | 38 | 19 | 520 ns total; slightly smaller system than EXP-04 |
| WHAM/MBAR analysis | 0.5 | 0.5 | CPU-bound post-processing |
| Energy decomposition (loop vs. scaffold) | 1 | 0.5 | Per-component energy comparison with EXP-04 |
| **Total** | **~42** | **~21** | §10.17: Tier 5 ≈ 42 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each US window; after equilibration |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-28/` |
| Restart procedure | Mount Drive → load latest checkpoint → verify energy drift < 0.1% → resume from next window |
| Estimated sessions needed | 2–3 (A100) or 1–2 (H100) |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Positions, velocities, box vectors, RNG state | OpenMM binary `.chk` | `checkpoint_<phase>_<step>.chk` |
| Full simulation state (portable) | OpenMM XML | `state_<phase>_<step>.xml` |
| Umbrella window ξ timeseries | NumPy `.npy` | `umbrella_window_<idx>.npy` |
| PMF results (peptide-trypsin) | NumPy `.npz` | `pmf_peptide_trypsin.npz` |

**Resume verification protocol:**
1. Reload checkpoint: `simulation.loadCheckpoint('checkpoint_<phase>_<step>.chk')`.
2. Run 1000 steps; compute energy.
3. Compare to energy at checkpoint save: drift must be < 0.1%.
4. If drift exceeds threshold, discard last 100 ps and re-equilibrate for 200 ps.

### 9.5 Platform Configuration

```python
from openmm import Platform
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
```

**Environment verification:**
```bash
!nvidia-smi
python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"
```

### 9.6 Dependency Notes

**Upstream dependency on EXP-04:**

| Required EXP-04 Output | Drive Path | Description |
|-----------------------|------------|-------------|
| ΔG_bind (full BPTI-trypsin) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/` | For ΔΔG_scaffold = ΔG(BPTI) − ΔG(loop peptide) |
| Equilibrated BPTI-trypsin complex (2PTC) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/equilibrated_complex.pdb` | Source for loop peptide extraction |
| Per-residue energy decomposition (from EXP-07) | `/content/drive/MyDrive/v3_gpu_results/EXP-07/` | Optional: compare scaffold contribution methods |

**Pre-execution check:** EXP-04 must be completed. The loop peptide is extracted from the EXP-04 equilibrated BPTI structure.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
