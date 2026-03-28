# EXP-06: SPINK1-Trypsin Binding Free Energy

**Experiment ID:** EXP-06  
**Feature ID:** F-06 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** KAZAL FAMILY BENCHMARK  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment computes the binding free energy of SPINK1 (identical to PSTI) against trypsin (Ki ≈ 7.6 nM, ΔG ≈ −11.1 kcal/mol), serving as a calibration point on the Kazal affinity ladder intermediate between SPINK7-KLK5 (−9.4 kcal/mol) and PSTI-chymotrypsinogen (−14.7 kcal/mol). The clinically significant N34S mutation does NOT alter Ki, providing a coupled negative control (see EXP-32). The trypsinogen-PSTI co-crystal structure (PDB 1TGS, 1.8 Å) provides the starting geometry, with trypsin (vs. trypsinogen) achieved by appropriate protonation and activation peptide handling.

---

## 2. Hypothesis

**H₁:** The V2 pipeline's US/WHAM estimate of ΔG_bind for SPINK1-trypsin will fall within the 95% CI [−15.6, −6.6] kcal/mol.

**H₂:** The computed ΔG will be correctly ranked between SPINK7-KLK5 (−9.4) and PSTI-chymotrypsinogen (−14.7).

---

## 3. Background and Rationale

SPINK1 is measured at 37°C (310 K) rather than 25°C (298 K), matching the pipeline's default temperature. Kuwata et al. (2002) used the Green & Work method at pH 8.0. The recombinant (Ki = 7.6 nM) and natural (Ki = 8.7 nM) forms agree within 0.1 kcal/mol, confirming measurement reliability. This intermediate-affinity system tests whether the pipeline can resolve different positions on the affinity scale.

---

## 4. Experimental Protocol

### 4.1 System Preparation

| Parameter | Value |
|-----------|-------|
| Complex structure | PDB 1TGS (trypsinogen-PSTI); trypsinogen → trypsin by activation peptide removal |
| Force field | AMBER ff14SB (`amber14-all.xml`) |
| Water model | TIP3P (`amber14/tip3p.xml`) |
| pH | 7.4 |
| Box padding | 1.2 nm |
| Ionic strength | 0.15 M NaCl |
| Box shape | Cubic |

### 4.2 Minimization, Equilibration, Production

- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT: 500 ps, 310 K, friction 1.0 ps⁻¹, timestep 2 fs
- NPT: 1000 ps, 310 K, 1.0 atm, barostat 25 steps
- Restraint: 1000 kJ/mol/nm²
- Production: 100 ns, save 10 ps

### 4.3 Umbrella Sampling

| Parameter | Value | Source |
|-----------|-------|--------|
| ξ range | 1.5–4.0 nm (51 windows) | `config.py: UmbrellaConfig` |
| Window spacing | 0.05 nm | `config.py: UmbrellaConfig.window_spacing_nm` |
| Spring constant | 1000 kJ/mol/nm² | `config.py: UmbrellaConfig.spring_constant_kj_mol_nm2` |
| Per-window duration | 10.0 ns | `config.py: UmbrellaConfig.per_window_duration_ns` |
| Pre-positioning velocity | 0.01 nm/ps | `config.py: UmbrellaConfig.pre_position_velocity_nm_per_ps` |
| Equilibration per window | 200 ps | `config.py: UmbrellaConfig.equilibration_duration_ps` |
| Save interval | 1.0 ps | `config.py: UmbrellaConfig.save_interval_ps` |

### 4.4 WHAM and MBAR

- WHAM: tolerance 10⁻⁶, max iterations 100,000, bootstraps 200, bins 200
- MBAR: solver "robust", tolerance 10⁻⁷, max iterations 10,000, bootstraps 200, bins 200

### 4.5 Statistical Comparison

| Component | Estimate |
|-----------|----------|
| σ_exp | 0.5 kcal/mol |
| σ_method | 2.0 kcal/mol |
| σ_combined | 2.3 kcal/mol |
| 95% CI | [−15.6, −6.6] |

---

## 5. Control Conditions

### 5.1 Positive Control

**EXP-05 (PSTI-chymotrypsinogen):** Same inhibitor, different protease. Comparison reveals protease-specificity effects.

### 5.2 Negative Control

**EXP-32 (SPINK1 N34S):** Same system with a mutation that experimentally has zero effect on binding. Must predict ΔΔG ≈ 0.

---

## 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind | −11.1 kcal/mol | [−15.6, −6.6] |
| T_simulation | 310 K (matches experimental 37°C) | — |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Trypsinogen vs. trypsin** | Zymogen activation loop affects binding | Structure preparation issue | Medium |
| **Temperature mismatch** | 310 K simulation vs. 310 K experiment (matched) | None expected | Low |
| **Ranking error** | Placed incorrectly on ladder | Pipeline resolution insufficient | High |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After structure prep | Trypsinogen correctly converted to trypsin form | Active site accessible |
| After equilibration | All invariants pass | Stable complex |
| After US/WHAM | Smooth PMF, converged | IV-8, IV-9 satisfied |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (SPINK1-trypsin complex: ~50,000 atoms solvated in TIP3P with 1.2 nm padding) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation | < 0.1 | < 0.1 | CPU-bound; 1TGS trypsinogen→trypsin conversion |
| Equilibration (NVT + NPT) | 0.5 | 0.3 | 1.5 ns total |
| Production MD (100 ns) | 3 | 1.5 | ~800 ns/day (A100) for ~50k atoms |
| Umbrella Sampling (51 windows × 10.2 ns) | 45 | 22 | 520 ns total; sequential windows |
| WHAM/MBAR analysis | 0.5 | 0.5 | CPU-bound post-processing |
| **Total** | **~52** | **~26** | §10.17: Tier 1 ≈ 50 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each US window; after each equilibration phase |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-06/` |
| Restart procedure | Mount Drive → load latest checkpoint → verify energy drift < 0.1% → resume from next incomplete phase |
| Estimated sessions needed | 3–4 (A100) or 2 (H100) |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Positions, velocities, box vectors, RNG state | OpenMM binary `.chk` | `checkpoint_<phase>_<step>.chk` |
| Full simulation state (portable) | OpenMM XML serialization | `state_<phase>_<step>.xml` |
| Umbrella window ξ timeseries | NumPy `.npy` | `umbrella_window_<idx>.npy` |

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

EXP-06 has no upstream dependencies. EXP-06 results contribute to:

| Downstream Experiment | Usage |
|-----------------------|-------|
| EXP-32 (SPINK1 N34S control) | Uses same SPINK1-trypsin system as baseline for ΔΔG |
| EXP-26 (BSA–ΔG correlation) | ΔG data point on Kazal affinity ladder |

All outputs must be persisted to Google Drive immediately after generation.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
