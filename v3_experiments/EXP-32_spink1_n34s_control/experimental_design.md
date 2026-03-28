# EXP-32: SPINK1 N34S Negative Control (ΔΔG ≈ 0)

**Experiment ID:** EXP-32  
**Feature ID:** F-32 (benchmarks.md)  
**Category:** Mutational  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment serves as a critical negative control for the FEP methodology. The SPINK1 N34S mutation is a common clinical variant (associated with chronic pancreatitis susceptibility), but structural and kinetic studies show it does NOT significantly alter trypsin binding (ΔΔG ≈ 0, Kuwata et al. 2002, Kiraly et al. 2007). The pipeline should correctly predict |ΔΔG(N34S)| < 1.0 kcal/mol. This tests the pipeline's ability to distinguish functionally neutral mutations from hot-spot mutations (EXP-30), an essential capability for clinical variant interpretation.

---

## 2. Hypothesis

**H₁:** FEP-calculated ΔΔG for SPINK1 N34S will be < 1.0 kcal/mol in magnitude (functionally neutral for binding).

**H₂:** The structural metrics (loop RMSD, H-bond count, BSA) will be unchanged between WT and N34S.

---

## 3. Background and Rationale

The N34S variant in SPINK1 (exon 3, c.101A>G) is present in ~2% of the population and is the most common SPINK1 variant associated with pancreatitis risk. However, in vitro studies consistently show that N34S SPINK1 inhibits trypsin with the same Ki as WT SPINK1. The pancreatitis association is likely through expression/secretion effects, not direct binding impairment. This makes N34S an ideal negative control: the pipeline must correctly predict that this mutation is binding-neutral.

N34 is located in the SPINK1 scaffold, away from the P1 loop and the trypsin interface. Mutation to Ser (N→S: removing the amide, adding a hydroxyl) should minimally perturb the interface.

---

## 4. Experimental Protocol

### 4.1 System

SPINK1-trypsin complex from EXP-06 (or modeled from PSTI-chymotrypsinogen PDB 1TGS with appropriate modifications).

### 4.2 FEP: N34S

**Leg 1: Complex (SPINK1-trypsin)**
1. Hybrid topology: N34 → S34 (Asn → Ser).
2. Lambda windows: n_lambda = 20 (0.0 to 1.0).
3. Per-window: 2.0 ns.
4. Temperature: 310 K.
5. Soft-core: alpha = 0.5, power = 1.
6. Save interval: 1.0 ps.
7. n_equil: 5000 steps.

**Leg 2: Free inhibitor (SPINK1 alone)**
1. Same FEP protocol on free SPINK1.

**ΔΔG = ΔG_complex(N→S) − ΔG_free(N→S)**

### 4.3 Analysis

- MBAR: solver "robust", tolerance 10⁻⁷, max_iterations 10,000.
- Bootstrap: n_bootstrap = 200.
- 95% CI from bootstrap.

### 4.4 Structural Comparison

Run 20 ns production MD of N34S mutant complex:
1. Loop RMSD vs WT: should be unchanged.
2. H-bond count: should be unchanged.
3. BSA: should be unchanged.
4. P1 in S1: should be maintained.

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | |ΔΔG| (kcal/mol) |
|---------------|------------------|
| **PASS** | < 1.0 |
| **MARGINAL** | 1.0–2.0 |
| **FAIL** | > 2.0 |

---

## 5. Control Conditions

### 5.1 Positive Control: K15A (EXP-30)

K15A ΔΔG ≈ +10 kcal/mol. The pipeline must clearly distinguish K15A (massive effect) from N34S (no effect) — at least 8 kcal/mol difference.

### 5.2 Location Control

Verify that N34 is not at the binding interface:
- d(N34, nearest trypsin atom) > 8 Å in the equilibrated complex.
- N34 BSA contribution ≈ 0 Å².
- No N34-trypsin H-bonds.

### 5.3 Self-Perturbation Control

N→N (identity transformation): ΔΔG should be exactly 0. This validates the FEP setup.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| ΔΔG(N34S) | ~0 kcal/mol | Kuwata 2002, Kiraly 2007 |
| |ΔΔG| bound | < 1.0 kcal/mol | benchmarks.md |
| N34-trypsin distance | > 8 Å | Not at interface |
| Structural change | None measurable | Neutral mutation |
| K15A vs N34S difference | > 8 kcal/mol | Hot spot vs neutral |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **False positive** | |ΔΔG| > 2 kcal/mol | FEP noise for small perturbation | High |
| **Sampling noise** | Bootstrap CI spans 0 but is wide (>3 kcal/mol) | Insufficient sampling | Medium |
| **Long-range electrostatic artifact** | N→S charge change propagates | PME finite-size effects | Medium |
| **SPINK1 model error** | N34 position wrong in homology model | Structural inaccuracy | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| N34 location | Not at interface (>8 Å from trypsin) | Scaffold position confirmed |
| Perturbation size | N→S: small chemical change (amide → hydroxyl) | Minimal perturbation |
| Lambda overlap | All adjacent windows well-overlapped | MBAR reliable |
| ΔΔG magnitude | |ΔΔG| < 1.0 kcal/mol | Neutral mutation |
| Structural metrics unchanged | ΔLoop_RMSD < 0.3 Å, ΔH-bond < 1, ΔBSA < 50 Å² | No structural effect |
| Discriminates from K15A | ΔΔG(N34S) << ΔΔG(K15A) | Hot spot vs neutral |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (SPINK1-trypsin complex: ~50,000 atoms; FEP lambda windows with alchemical transformations) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation (N34S mutation) | < 0.1 | < 0.1 | CPU-bound; mutate Asn34→Ser |
| Equilibration (NVT + NPT, WT + mutant) | 0.5 | 0.3 | 1.5 ns per leg × 2 legs |
| FEP forward leg: complex (20λ × 2 ns) | 2.5 | 1.25 | 40 ns total |
| FEP forward leg: solvent (20λ × 2 ns) | 2.5 | 1.25 | 40 ns total |
| BAR/MBAR analysis | 0.3 | 0.3 | CPU-bound free energy estimation |
| **Total** | **~6** | **~3** | §10.17: Tier 5 ≈ 6 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each λ window (every 2 ns); after equilibration |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-32/` |
| Restart procedure | Mount Drive → identify last completed λ window → verify energy drift < 0.1% → resume from next λ |
| Estimated sessions needed | 1 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| FEP window state | OpenMM binary `.chk` | `checkpoint_<leg>_lambda<idx>.chk` |
| Full simulation state | OpenMM XML | `state_<leg>_lambda<idx>.xml` |
| Per-window ΔU samples | NumPy `.npy` | `du_<leg>_lambda<idx>.npy` |
| Combined FEP results | NumPy `.npz` | `fep_n34s_results.npz` |

**Resume verification protocol:**
1. Reload checkpoint: `simulation.loadCheckpoint('checkpoint_<leg>_lambda<idx>.chk')`.
2. Run 1000 steps; compute energy.
3. Compare to energy at checkpoint save: drift must be < 0.1%.
4. If drift exceeds threshold, discard window and re-equilibrate at that λ for 200 ps.

### 9.5 Platform Configuration

```python
from openmm import Platform
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
# FEP requires careful lambda scheduling:
# for lam in np.linspace(0, 1, 20):
#     system = set_lambda(system, lam)  # Scale electrostatics + LJ
```

**Environment verification:**
```bash
!nvidia-smi
python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"
```

### 9.6 Dependency Notes

**Upstream dependency on EXP-06:**

| Required EXP-06 Output | Drive Path | Description |
|-----------------------|------------|-------------|
| Equilibrated SPINK1-trypsin complex | `/content/drive/MyDrive/v3_gpu_results/EXP-06/production/equilibrated_complex.pdb` | Starting structure for N34S mutation |
| ΔG_bind (WT SPINK1-trypsin) | `/content/drive/MyDrive/v3_gpu_results/EXP-06/analysis/` | Baseline for ΔΔG comparison |

**Pre-execution check:** EXP-06 must be completed. The N34S mutation is applied to the EXP-06 equilibrated structure. Expected result: |ΔΔG| < 1.0 kcal/mol (negative control).

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
