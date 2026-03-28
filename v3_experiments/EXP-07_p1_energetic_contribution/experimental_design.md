# EXP-07: P1 Residue Energetic Contribution to Binding

**Experiment ID:** EXP-07  
**Feature ID:** F-07 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** SEMI-QUANTITATIVE STRUCTURAL/ENERGETIC  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment tests whether the V2 pipeline identifies the P1 residue as the dominant energetic contributor to protease-inhibitor binding. The P1 residue (Arg or Lys at the scissile bond) accounts for ~70% of the total association energy and ~50% of the interface contact area in canonical serine protease-inhibitor complexes (Krowarsch et al. 2003). Using per-residue energy decomposition on the equilibrated BPTI-trypsin (PDB 2PTC) and/or SPINK7-KLK5 complex, the pipeline should identify P1 as the single largest per-residue contributor. This is independently confirmed by the K15A mutation in BPTI (ΔΔG ≈ 10 kcal/mol, Castro & Anderson 1996).

---

## 2. Hypothesis

**H₁:** Per-residue energy decomposition of the protease-inhibitor interface will identify the P1 residue (K15 in BPTI, Arg in SPINK7) as the single largest contributor to ΔG_bind, accounting for >40% of the total binding energy.

**H₂:** The P1 residue will contribute >40% of the total interface contact surface area.

---

## 3. Background and Rationale

The Laskowski standard mechanism places the P1 residue at the core of protease-inhibitor recognition. For trypsin-like proteases (KLK5, trypsin), P1 is Arg or Lys, forming a salt bridge with Asp189 in the S1 specificity pocket. This experiment validates the pipeline's ability to decompose binding energetics at the per-residue level — essential for identifying mutational hot spots and validating the physical realism of the computed binding interface.

---

## 4. Experimental Protocol

### 4.1 System

Use the equilibrated BPTI-trypsin complex from EXP-04 (PDB 2PTC) and/or the SPINK7-KLK5 complex from EXP-01.

### 4.2 Per-Residue Energy Decomposition

1. From the last 50 ns of production MD trajectory (to ensure equilibration), extract frames at 100 ps intervals (500 frames).
2. For each frame, compute per-residue interaction energies between the inhibitor and protease:
   - Electrostatic contribution: Coulombic interaction between each inhibitor residue and all protease atoms.
   - van der Waals contribution: Lennard-Jones interaction for same pairs.
   - Total: sum of electrostatic + vdW for each residue.
3. Time-average across all frames.
4. Rank residues by magnitude of interaction energy.

### 4.3 Contact Surface Area Decomposition

1. For each frame, compute per-residue BSA (buried surface area upon binding).
2. BSA_residue = SASA_residue(free) − SASA_residue(complex).
3. Time-average and rank.

### 4.4 Simulation Parameters

All parameters from the production run of EXP-04:
- Temperature: 310 K
- Force field: AMBER ff14SB
- Timestep: 2 fs
- Analysis frames: last 50 ns at 100 ps intervals

### 4.5 Statistical Comparison

Semi-quantitative benchmark:
- **PASS:** P1 is the #1 ranked residue in energy decomposition AND contributes >40% of total binding energy.
- **MARGINAL:** P1 is in the top 3 residues, contributing 30–40%.
- **FAIL:** P1 is not the top contributor or contributes <30%.

---

## 5. Control Conditions

### 5.1 Positive Control

**BPTI K15A data (F-30, EXP-30):** K15A ΔΔG ≈ 10 kcal/mol experimentally. The energy decomposition should show K15 contributing ≥10 kcal/mol to the interface interaction energy.

### 5.2 Negative Control

Peripheral residues (e.g., G36 in BPTI, where G36A ΔΔG ≈ 0) should show near-zero interface interaction energy.

---

## 6. Expected Outcomes

| Metric | Expected Value |
|--------|---------------|
| P1 rank in energy decomposition | #1 |
| P1 fraction of total binding energy | ~70% (accept >40%) |
| P1 fraction of interface BSA | ~50% (accept >30%) |
| Dominant interaction at P1 | Salt bridge with Asp189 |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **P1 salt bridge disrupted** | K15/Arg not in S1 pocket | Equilibration failure | Critical |
| **Force field electrostatic artifacts** | Non-P1 residue dominates due to partial charge artifacts | Fixed-charge limitation | Medium |
| **Decomposition scheme dependence** | Results vary with cutoff or decomposition method | Methodological sensitivity | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After trajectory extraction | 500 frames from production MD | Correct frame count |
| After energy decomposition | Total interaction energy negative | Attractive interaction |
| Visual check | P1 residue in S1 pocket in all frames | Salt bridge maintained |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~2–4 GB (energy decomposition on pre-computed trajectory; no live simulation of full system) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Load EXP-04 trajectory (last 50 ns, 500 frames) | 0.1 | 0.1 | I/O from Drive |
| Per-residue energy decomposition | 1.5 | 1.0 | Re-evaluation of nonbonded energies per-residue |
| Statistical analysis (bootstrap, ranking) | 0.5 | 0.5 | CPU-bound |
| **Total** | **~2–3** | **~1.5–2** | §10.17: Tier 2 ≈ 2–3 GPU-hrs |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each trajectory batch (every 100 frames) |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-07/` |
| Restart procedure | Mount Drive → load partial decomposition results → resume from next unprocessed frame |
| Estimated sessions needed | 1 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Per-residue energy arrays (partial) | NumPy `.npy` | `energy_decomp_frames_<start>_<end>.npy` |
| Final decomposition results | NumPy `.npz` | `p1_energy_decomposition.npz` |

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

**CRITICAL — Tier 2 upstream dependency on EXP-04:**

| Required EXP-04 Output | Drive Path | Description |
|-----------------------|------------|-------------|
| Production MD trajectory (last 50 ns) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/trajectory_last50ns.dcd` | 500 frames at 100 ps intervals |
| Equilibrated complex structure | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/equilibrated_complex.pdb` | Topology file |
| System XML | `/content/drive/MyDrive/v3_gpu_results/EXP-04/system/system.xml` | Force field parameters for energy re-evaluation |

**Pre-execution check:** Verify all EXP-04 outputs exist on Drive before starting EXP-07. If missing, EXP-04 must be completed first.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
