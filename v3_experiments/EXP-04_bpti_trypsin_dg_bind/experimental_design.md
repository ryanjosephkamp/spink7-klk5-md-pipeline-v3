# EXP-04: BPTI-Trypsin Binding Free Energy

**Experiment ID:** EXP-04  
**Feature ID:** F-04 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** GOLD STANDARD REFERENCE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment determines the binding free energy of the BPTI-trypsin complex, the gold-standard protease-inhibitor system (Kd = 6 × 10⁻¹⁴ M, ΔG ≈ −18 kcal/mol), using US/WHAM and SMD/Jarzynski methods. Unlike EXP-01, this system benefits from a high-resolution co-crystal structure (PDB 2PTC, 1.9 Å), removing model quality as a confounding variable. This is the definitive test of the pipeline's free energy methodology for protease-inhibitor systems — success here validates the method; failure indicates fundamental methodological limitations independent of model quality.

---

## 2. Hypothesis

**H₁:** The V2 pipeline's US/WHAM estimate of ΔG_bind for BPTI-trypsin will fall within the 95% CI [−24.7, −11.3] kcal/mol.

**H₂:** The wide CI reflects the extreme difficulty of computationally reproducing ultra-tight binding (−18 kcal/mol). The pipeline may achieve MARGINAL classification due to the fundamental challenge of sampling the deep free energy basin required for femtomolar affinity.

---

## 3. Background and Rationale

### 3.1 Scientific Context

BPTI (58-residue Kunitz-type inhibitor) and trypsin form the most extensively characterized protease-inhibitor complex. Two independent labs confirm Kd ≈ 6 × 10⁻¹⁴ M (Vincent & Lazdunski 1972, Castro & Anderson 1996). The co-crystal structure (PDB 2PTC at 1.9 Å) provides an ideal starting point. While BPTI is Kunitz-type (not Kazal), it shares the canonical serine protease inhibition mechanism with SPINK7.

### 3.2 What This Reveals

As the gold standard with no model uncertainty, this experiment isolates method accuracy from model quality. A PASS validates the computational methodology; a FAIL attributable to sampling or force field is a fundamental pipeline limitation.

---

## 4. Experimental Protocol

### 4.1 System Preparation

| Parameter | Value |
|-----------|-------|
| Complex structure | PDB 2PTC (BPTI-trypsin co-crystal, 1.9 Å) |
| Force field | AMBER ff14SB (`amber14-all.xml`) |
| Water model | TIP3P (`amber14/tip3p.xml`) |
| pH | 7.4 |
| Box padding | 1.2 nm |
| Ionic strength | 0.15 M NaCl |
| Box shape | Cubic |

### 4.2 Minimization and Equilibration

| Parameter | Value | Source |
|-----------|-------|--------|
| Minimization | 10,000 steps, tolerance 10.0 kJ/mol/nm | `config.py: MinimizationConfig` |
| NVT | 500 ps, 310 K, friction 1.0 ps⁻¹, timestep 2 fs | `config.py: EquilibrationConfig` |
| NPT | 1000 ps, 310 K, 1.0 atm, barostat 25 steps | `config.py: EquilibrationConfig` |
| Restraint | 1000 kJ/mol/nm², save interval 10 ps | `config.py: EquilibrationConfig` |

### 4.3 Production MD

| Parameter | Value | Source |
|-----------|-------|--------|
| Duration | 100 ns | `config.py: ProductionConfig.duration_ns` |
| Temperature | 310 K | `config.py: ProductionConfig.temperature_k` |
| Timestep | 0.002 ps | `config.py: ProductionConfig.timestep_ps` |
| Save interval | 10 ps | `config.py: ProductionConfig.save_interval_ps` |

### 4.4 SMD

| Parameter | Value | Source |
|-----------|-------|--------|
| Spring constant | 1000 kJ/mol/nm² | `config.py: SMDConfig.spring_constant_kj_mol_nm2` |
| Pulling velocity | 0.001 nm/ps | `config.py: SMDConfig.pulling_velocity_nm_per_ps` |
| Pull distance | 3.0 nm | `config.py: SMDConfig.pull_distance_nm` |
| Replicates | 50 | `config.py: SMDConfig.n_replicates` |
| Save interval | 1.0 ps | `config.py: SMDConfig.save_interval_ps` |
| Total time | 150 ns | |

### 4.5 Umbrella Sampling

| Parameter | Value | Source |
|-----------|-------|--------|
| ξ range | 1.5–4.0 nm (51 windows at 0.05 nm spacing) | `config.py: UmbrellaConfig` |
| Spring constant | 1000 kJ/mol/nm² | `config.py: UmbrellaConfig.spring_constant_kj_mol_nm2` |
| Per-window duration | 10.0 ns | `config.py: UmbrellaConfig.per_window_duration_ns` |
| Pre-positioning velocity | 0.01 nm/ps | `config.py: UmbrellaConfig.pre_position_velocity_nm_per_ps` |
| Equilibration per window | 200 ps | `config.py: UmbrellaConfig.equilibration_duration_ps` |
| Save interval | 1.0 ps | `config.py: UmbrellaConfig.save_interval_ps` |
| Total time | 510 ns | |

### 4.6 WHAM and MBAR Analysis

- WHAM: tolerance 10⁻⁶, max iterations 100,000, bootstraps 200, bins 200
- MBAR: solver "robust", tolerance 10⁻⁷, max iterations 10,000, bootstraps 200, bins 200

### 4.7 Statistical Comparison

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.5 kcal/mol | Two labs agree within 0.1 kcal/mol |
| σ_comp | from bootstrap | Pipeline internal uncertainty |
| σ_method | 3.0 kcal/mol | Ultra-tight binding is extremely challenging (§25.4) |
| σ_combined | 3.4 kcal/mol | √(0.5² + 1.5² + 3.0²) — assuming σ_comp ≈ 1.5 |
| 95% CI | [−24.7, −11.3] | −18.0 ± (1.96 × 3.4) |

---

## 5. Control Conditions

### 5.1 Positive Control

**EXP-29 (SH3-p41):** Method validation target where three independent computational approaches agree within 0.2 kcal/mol of experiment. If the pipeline can reproduce SH3-p41 but not BPTI-trypsin, the limitation is specific to tight protease-inhibitor binding.

### 5.2 Negative Control

1. All invariants IV-1 through IV-10 satisfied.
2. BPTI structural stability: backbone RMSD < 1.5 Å from 2PTC crystal structure (F-23).
3. K15 (P1) remains in S1 pocket throughout equilibration.
4. Disulfide bonds (Cys5-Cys55, Cys14-Cys38, Cys30-Cys51) remain intact (IV-6).

---

## 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind (US/WHAM) | −18.0 kcal/mol | [−24.7, −11.3] |
| PMF minimum depth | Deep, well-defined bound state | >10 kcal/mol below bulk |

### Classification

- **PASS:** ΔG within [−24.7, −11.3] kcal/mol
- **MARGINAL:** ΔG within [−31.4, −4.6] kcal/mol (2× CI)
- **FAIL:** ΔG outside [−31.4, −4.6] kcal/mol

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation Implied | Severity |
|-------------|--------------|-------------------|----------|
| **ΔG insufficiently negative** | ΔG > −11 kcal/mol | Method cannot capture ultra-tight binding | High |
| **Incomplete dissociation in SMD** | Plateau not reached in 3 nm pull | Pull distance insufficient for tight complex | Medium — increase pull distance |
| **WHAM convergence issues** | Deep bound-state well causes histogram undersampling | Need more windows at bound state | Medium |
| **Pulling pathway artifacts** | Non-optimal unbinding pathway | COM pulling may not find minimum-work path | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After structure loading | 2PTC correctly parsed; both chains present | BPTI + trypsin resolved |
| After minimization | IV-1 satisfied | Energy decreased |
| After equilibration | All invariants pass; K15 in S1 pocket | Stable complex |
| After production MD | RMSD < 2 Å; complex intact | Structural stability |
| After SMD | Work distributions unimodal (IV-10) | Clean pulling |
| After US | Histogram overlap ≥ 10% (IV-8) | Sufficient sampling |
| After WHAM | Convergence < 10⁻⁶ (IV-9) | WHAM converged |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (BPTI-trypsin complex: ~35,000 atoms solvated in TIP3P with 1.2 nm padding) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation | < 0.1 | < 0.1 | CPU-bound; negligible |
| Equilibration (NVT + NPT) | 0.5 | 0.3 | 1.5 ns total |
| Production MD (100 ns) | 3 | 1.5 | ~800 ns/day (A100), ~1500 ns/day (H100) for 35k atoms |
| SMD (50 replicates × 3 ns) | 5 | 2.5 | 150 ns total; parallel replicates |
| Umbrella Sampling (51 windows × 10.2 ns) | 40 | 20 | 520 ns total; sequential windows |
| WHAM/MBAR analysis | 0.5 | 0.5 | CPU-bound post-processing |
| **Total** | **~52** | **~26** | §10.17: Tier 1 ≈ 50 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each US window; after every 10 SMD replicates; after each equilibration phase |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-04/` |
| Restart procedure | Mount Drive → load latest checkpoint → verify energy drift < 0.1% → resume from next incomplete phase |
| Estimated sessions needed | 3–4 (A100) or 2 (H100) |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Positions, velocities, box vectors, RNG state | OpenMM binary `.chk` | `checkpoint_<phase>_<step>.chk` |
| Full simulation state (portable) | OpenMM XML serialization | `state_<phase>_<step>.xml` |
| Umbrella window ξ timeseries | NumPy `.npy` | `umbrella_window_<idx>.npy` |
| SMD work trajectories | NumPy `.npy` | `smd_work_replicate_<idx>.npy` |

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
# Pass to Simulation constructor:
# simulation = Simulation(topology, system, integrator, platform, properties)
```

**Environment verification:**
```bash
!nvidia-smi                          # Confirm GPU allocated
python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"  # List platforms
```

### 9.6 Dependency Notes

EXP-04 has no upstream dependencies. However, EXP-04 outputs are required by four downstream experiments:

| Downstream Experiment | Required EXP-04 Output | Drive Path |
|-----------------------|-----------------------|------------|
| EXP-07 (P1 energetic contribution) | Last 50 ns production MD trajectory, equilibrated structure | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/` |
| EXP-08 (Interfacial H-bond energy) | Last 50 ns production MD trajectory | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/` |
| EXP-09 (Kinetics) | PMF from US/WHAM (ξ, W(ξ), bootstrap samples) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/` |
| EXP-10 (Activation energy) | PMF from US/WHAM (ξ, W(ξ), bootstrap samples) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/` |
| EXP-28 (Scaffold energy) | ΔG_bind result, equilibrated complex structure | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/` |

All outputs must be persisted to Google Drive immediately after generation.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
