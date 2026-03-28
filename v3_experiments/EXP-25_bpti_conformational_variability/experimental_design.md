# EXP-25: BPTI Conformational Variability (Cα RMSD Fluctuation)

**Experiment ID:** EXP-25  
**Feature ID:** F-25 (benchmarks.md)  
**Category:** Dynamic  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the conformational variability of BPTI from MD simulations, measured as the time-averaged Cα RMSD fluctuation around the mean structure. The experimental reference value is 0.40 Å (from NMR ensemble analysis, Berndt et al. 1992), with a 95% pipeline CI of [0.25, 0.55] Å. This metric directly reports on the magnitude of thermal motions and tests whether the force field produces physically realistic protein dynamics — neither too rigid (insufficient sampling) nor too flexible (poor parameterization).

---

## 2. Hypothesis

**H₁:** The time-averaged Cα RMSD fluctuation of BPTI around the MD mean structure will be 0.40 ± 0.15 Å (95% CI: [0.25, 0.55] Å).

**H₂:** The average Cα RMSF will be consistent with the Cα RMSD fluctuation, with core residues showing RMSF < 0.3 Å and loop regions showing RMSF > 0.5 Å.

---

## 3. Background and Rationale

Conformational variability quantifies the extent of structural excursions from the average conformation. For a well-folded, disulfide-stabilized protein like BPTI, the Cα RMSD fluctuation should be modest (~0.4 Å), reflecting thermal motions within the native basin rather than large-scale conformational changes. This is an orthogonal check to backbone RMSD from crystal (EXP-23) — here we measure the spread of the conformational ensemble rather than the deviation from a reference.

---

## 4. Experimental Protocol

### 4.1 System

Free BPTI (PDB 4PTI) from EXP-23 production trajectory.

### 4.2 Conformational Variability Calculation

1. From 100 ns production MD, extract all frames (10,000 frames at 10 ps intervals).
2. Compute the average structure: mean Cα coordinates across all frames (after backbone alignment).
3. For each frame, compute Cα RMSD relative to this mean structure.
4. The conformational variability = mean(RMSD_to_average) across all frames.

### 4.3 Per-Residue RMSF

1. Compute Cα RMSF (root mean square fluctuation) per residue.
2. Average RMSF across all residues → related to conformational variability.
3. Identify rigid core (RMSF < 0.3 Å) and flexible regions (RMSF > 0.8 Å).

### 4.4 Block Averaging

1. Divide trajectory into 5 blocks of 20 ns each.
2. Compute conformational variability per block.
3. Standard error across blocks → confidence interval for convergence.

### 4.5 Production MD Parameters (from config.py)

- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P
- Duration: 100 ns
- Temperature: 310 K, friction 1.0 ps⁻¹
- Timestep: 2 fs
- Save interval: 10 ps

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | Cα RMSD Fluctuation |
|---------------|--------------------|
| **PASS** | [0.25, 0.55] Å |
| **MARGINAL** | [0.15, 0.70] Å |
| **FAIL** | Outside marginal range |

---

## 5. Control Conditions

### 5.1 NMR Ensemble Comparison

Compare the MD conformational variability to the NMR ensemble variability from PDB 2LEO (SPINK7, 20 models). The NMR ensemble spread should be comparable to the MD thermal spread for well-ordered regions.

### 5.2 Temperature Scaling Control

At 310 K, variability should be slightly larger than at 300 K. The scaling should follow ~√(T₂/T₁) for harmonic motions.

### 5.3 Disulfide Impact

Compare variability of disulfide-connected regions vs. free loops. Disulfide-bridged regions should show lower variability (< 0.3 Å), confirming the constraining effect of covalent crosslinks.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Cα RMSD fluctuation | 0.40 ± 0.15 Å | Berndt 1992, benchmarks.md |
| 95% CI | [0.25, 0.55] Å | benchmarks.md |
| Core RMSF | < 0.3 Å | Disulfide-stabilized |
| Loop RMSF | 0.5–1.2 Å | Flexible regions |
| Block-to-block CV | < 15% | Well-sampled |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Too rigid** | Variability < 0.2 Å | Restrained/frozen DOF | Medium |
| **Too flexible** | Variability > 0.7 Å | Force field softness | Medium |
| **Non-convergence** | Block averages don't converge | Insufficient sampling | Medium |
| **Temperature mismatch** | NMR at 300K, MD at 310K | ~3% difference expected | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Trajectory alignment | All frames aligned to reference | No translation artifacts |
| Mean structure quality | Mean structure has no clashes | Physical average |
| Block convergence | 5 blocks agree within 0.05 Å | Well-sampled |
| RMSF profile | Core < 0.3 Å, loops > 0.5 Å | Expected pattern |
| Temperature check | MD at 310 K ~ expected for AMBER ff14SB | Literature comparison |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~1–2 GB (trajectory analysis only; no live simulation) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Load EXP-24 trajectory (100 ns, 10,000 frames) | 0.1 | 0.1 | I/O from Drive |
| Cα RMSD calculation (per-frame vs. crystal) | 0.2 | 0.2 | MDTraj analysis; CPU-bound |
| RMSF calculation + time-block analysis | 0.1 | 0.1 | Statistical convergence |
| Visualization/output | 0.1 | 0.1 | Plots + data export |
| **Total** | **~0.5** | **~0.5** | §10.17: Tier 3 ≈ 0.5 GPU-hrs |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | Not applicable (single short analysis) |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-25/` |
| Restart procedure | Re-run analysis from start (< 30 min total) |
| Estimated sessions needed | 1 (can share session with EXP-24) |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Cα RMSD timeseries | NumPy `.npy` | `ca_rmsd_timeseries.npy` |
| RMSF per residue | NumPy `.npy` | `rmsf_per_residue.npy` |
| Summary statistics | JSON | `conformational_variability.json` |

### 9.5 Platform Configuration

No GPU simulation required. Analysis uses MDTraj and NumPy (CPU). GPU only needed if trajectory loading benefits from GPU-accelerated I/O.

```python
# No OpenMM platform needed for this experiment
import mdtraj as md
traj = md.load('bpti_production_100ns.dcd', top='bpti_4pti.pdb')
```

### 9.6 Dependency Notes

**CRITICAL — Tier 3 upstream dependency on EXP-24:**

| Required EXP-24 Output | Drive Path | Description |
|-----------------------|------------|-------------|
| 100 ns production MD trajectory | `/content/drive/MyDrive/v3_gpu_results/EXP-24/production/bpti_production_100ns.dcd` | Full trajectory for RMSD/RMSF |
| Equilibrated structure (topology) | `/content/drive/MyDrive/v3_gpu_results/EXP-24/production/bpti_4pti.pdb` | Reference structure |

**Pre-execution check:** Verify EXP-24 trajectory exists on Drive before starting EXP-25. If missing, EXP-24 must be completed first.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
