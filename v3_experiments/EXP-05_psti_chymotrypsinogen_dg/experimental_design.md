# EXP-05: PSTI-Chymotrypsinogen Binding Free Energy

**Experiment ID:** EXP-05  
**Feature ID:** F-05 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** KAZAL ANALOG — CLOSEST STRUCTURAL REFERENCE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment computes the binding free energy of the PSTI (Pancreatic Secretory Trypsin Inhibitor / SPINK1)-chymotrypsinogen complex, the closest available Kazal-type structural analog to SPINK7-KLK5. The co-crystal structure at 2.3 Å (Hecht et al. 1991, PDB 1TGS for trypsinogen-PSTI) eliminates model quality as a variable, and the Ki = 1.6 × 10⁻¹¹ M (ΔG ≈ −14.7 kcal/mol) provides a benchmark intermediate between SPINK7-KLK5 (−9.4) and BPTI-trypsin (−18.0) on the affinity scale.

---

## 2. Hypothesis

**H₁:** The V2 pipeline's US/WHAM estimate of ΔG_bind for PSTI-chymotrypsinogen will fall within the 95% CI [−20.0, −9.4] kcal/mol.

**H₂:** The computed ΔG will place PSTI-chymotrypsinogen between SPINK7-KLK5 (−9.4) and BPTI-trypsin (−18.0) on the affinity ladder, validating the pipeline's ability to reproduce the quantitative ranking of protease-inhibitor affinities.

---

## 3. Background and Rationale

PSTI (SPINK1) is a Kazal-type inhibitor in the same family as SPINK7. The crystal structure of PSTI3 and PSTI4 variants in complex with bovine chymotrypsinogen A (Hecht et al. 1991) provides the most directly relevant structural and thermodynamic reference. The binding loop rmsd of 0.25–0.35 Å across the Kazal family confirms structural conservation. This experiment calibrates the pipeline's performance on a well-characterized Kazal-protease system with experimental crystal structure data.

---

## 4. Experimental Protocol

### 4.1 System Preparation

| Parameter | Value |
|-----------|-------|
| Complex structure | PDB 1TGS (trypsinogen-PSTI complex, 1.8 Å) |
| Force field | AMBER ff14SB (`amber14-all.xml`) |
| Water model | TIP3P (`amber14/tip3p.xml`) |
| pH | 7.4 |
| Box padding | 1.2 nm |
| Ionic strength | 0.15 M NaCl |
| Box shape | Cubic |

### 4.2 Minimization, Equilibration, Production

Parameters identical to EXP-01 (§4.1.3, §4.1.4, §4.2):
- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT: 500 ps, 310 K; NPT: 1000 ps, 310 K, 1.0 atm
- Production: 100 ns, timestep 2 fs, save 10 ps

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

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.3 kcal/mol | Two PSTI variants agree within 0.2 kcal/mol |
| σ_comp | from bootstrap | Pipeline uncertainty |
| σ_method | 2.5 kcal/mol | Tight Kazal-protease, US/WHAM |
| σ_combined | 2.7 kcal/mol | √(0.3² + 1.0² + 2.5²) |
| 95% CI | [−20.0, −9.4] | −14.7 ± (1.96 × 2.7) |

---

## 5. Control Conditions

### 5.1 Positive Control

**EXP-04 (BPTI-trypsin):** Also uses co-crystal starting structure. Pair comparison reveals whether Kazal vs. Kunitz fold affects pipeline accuracy.

### 5.2 Negative Control

1. All invariants satisfied (IV-1 through IV-9).
2. Kazal fold integrity maintained (binding loop RMSD < 1.0 Å from crystal).
3. ΔG should be negative.

---

## 6. Expected Outcomes

| Metric | Expected Value | 95% CI |
|--------|---------------|--------|
| ΔG_bind | −14.7 kcal/mol | [−20.0, −9.4] |
| Affinity ranking | BPTI-trypsin > PSTI-chymotrypsinogen > SPINK7-KLK5 | Correct order |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Chymotrypsinogen activation** | Zymogen form less stable | Simulation of zymogen vs. active protease | Medium |
| **Force field bias** | Systematic over/underestimation | AMBER ff14SB limitation | High |
| **Ranking inversion** | Placed incorrectly on affinity scale | Pipeline lacks resolution | High |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| After structure loading | 1TGS correctly parsed | Both chains present |
| After equilibration | All invariants pass | Stable complex |
| After US/WHAM | PMF smooth with clear minimum | Converged |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (PSTI-chymotrypsinogen complex: ~50,000 atoms solvated in TIP3P with 1.2 nm padding) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation | < 0.1 | < 0.1 | CPU-bound; negligible |
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
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-05/` |
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

EXP-05 has no upstream dependencies. EXP-05 results contribute to:

| Downstream Experiment | Usage |
|-----------------------|-------|
| EXP-26 (BSA–ΔG correlation) | ΔG data point on affinity ladder |

All outputs must be persisted to Google Drive immediately after generation.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
