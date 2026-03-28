# EXP-29: SH3–p41 Method Validation (ΔG_bind)

**Experiment ID:** EXP-29  
**Feature ID:** F-29 (benchmarks.md)  
**Category:** Biophysical / Method Validation  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates the pipeline's free-energy methods on the Fyn SH3–p41 peptide system, a well-characterized benchmark with ΔG_bind = −7.99 kcal/mol (Isvoran et al. 2018, PDB 4EIK). Unlike the protease-inhibitor systems, SH3-p41 is a smaller, simpler protein-peptide interaction that has been extensively used for PMF/US/SMD method validation. Successfully reproducing this binding free energy validates the pipeline's methodology independently of the protease-inhibitor structural biology, providing an out-of-family positive control.

---

## 2. Hypothesis

**H₁:** The US/WHAM-derived ΔG_bind for Fyn SH3–p41 will fall within the 95% CI of [−12.7, −3.3] kcal/mol.

**H₂:** The SMD/Jarzynski-derived ΔG_bind will independently fall within the same confidence interval.

**H₃:** Both methods will agree within their respective bootstrap uncertainties.

---

## 3. Background and Rationale

Fyn SH3–p41 is a standard benchmark for PMF-based binding free energy calculations. The interaction is moderate-affinity (ΔG = −7.99 kcal/mol), involves a short proline-rich peptide binding to the SH3 domain via a well-defined PPII helix binding mode, and has a high-resolution co-crystal structure (1.5 Å). This system serves as an orthogonal validation: if the pipeline's methods work for both SH3-p41 and BPTI-trypsin, the approach is robust across different binding geometries and affinities.

---

## 4. Experimental Protocol

### 4.1 System Preparation

- PDB: 4EIK (Fyn SH3 domain in complex with p41 peptide)
- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P (`amber14/tip3p.xml`)
- pH: 7.4
- Box padding: 1.2 nm, cubic
- Ionic strength: 0.15 M NaCl

### 4.2 Equilibration (from config.py)

- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT: 500 ps at 310 K, friction 1.0 ps⁻¹, restraints 1000 kJ/mol/nm²
- NPT: 1000 ps at 310 K, 1.0 atm, barostat interval 25 steps
- Timestep: 2 fs

### 4.3 SMD Protocol

- Spring constant k: 1000 kJ/mol/nm²
- Pulling velocity v: 0.001 nm/ps
- Pull distance: 3.0 nm
- Number of replicates: 50
- Save interval: 1.0 ps
- ΔG estimation: Jarzynski equality with cumulant expansion (2nd order):
  $$\Delta G \approx \langle W \rangle - \frac{\beta}{2} \sigma_W^2$$

### 4.4 Umbrella Sampling Protocol

- Reaction coordinate ξ: COM distance SH3 ↔ p41 peptide
- ξ range: 1.5–4.0 nm, spacing 0.05 nm (51 windows)
- Spring constant: 1000 kJ/mol/nm²
- Per-window: 10 ns (after 200 ps equilibration)
- Pre-positioning velocity: 0.01 nm/ps
- Save interval: 1.0 ps
- WHAM: tolerance 10⁻⁶, max_iterations 100,000, bins 200
- MBAR: solver "robust", tolerance 10⁻⁷, max_iterations 10,000
- Bootstrap: n_bootstrap = 200

### 4.5 ΔG Extraction

- PMF well depth + volume correction (C° = 1/1660 Å⁻³)
- Bootstrap 95% CI from WHAM/MBAR

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | ΔG_bind (kcal/mol) |
|---------------|--------------------|
| **PASS** | [−12.7, −3.3] |
| **MARGINAL** | [−16.0, −1.0] |
| **FAIL** | Outside marginal range or wrong sign |

---

## 5. Control Conditions

### 5.1 Positive Control: BPTI-Trypsin (EXP-04)

Pipeline should correctly rank-order: |ΔG(BPTI-trypsin)| > |ΔG(SH3-p41)|.

### 5.2 Method Agreement Control

US/WHAM and SMD/Jarzynski should agree within 3 kcal/mol. Large disagreements indicate method-specific artifacts.

### 5.3 Convergence Control

First-half vs. second-half ΔG from US should agree within 2 kcal/mol (convergence check).

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| ΔG_bind (experimental) | −7.99 kcal/mol | Isvoran 2018 |
| 95% CI | [−12.7, −3.3] | benchmarks.md |
| US/WHAM ΔG | Within CI | Pipeline target |
| SMD/Jarzynski ΔG | Within CI | Pipeline target |
| SMD-US agreement | Within 3 kcal/mol | Method consistency |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Peptide flexibility** | p41 unfolds during pulling | Short peptide instability | Medium |
| **End-state sampling** | PMF doesn't converge at large ξ | Peptide conformational change | Medium |
| **Force field bias** | PPII helix preference over/under-estimated | ff14SB peptide parameters | Medium |
| **Jarzynski dissipation** | SMD work variance too large | Insufficient replicates | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| PDB quality | 4EIK loads, all atoms present | No missing residues |
| Equilibration | Complex RMSD < 2.0 Å | Stable |
| SMD work distribution | σW / ⟨W⟩ < 1.0 | Reasonable dissipation |
| US histogram overlap | Adjacent windows overlap >10% | Adequate sampling |
| WHAM convergence | Free energy converged within tolerance | Converged PMF |
| WHAM vs MBAR | Agreement within 1.5 kcal/mol | Internal consistency |
| Method rank | |ΔG(BPTI)| > |ΔG(SH3-p41)| | Correct ordering |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~3–5 GB (SH3-p41 complex: ~25,000 atoms solvated in TIP3P with 1.2 nm padding; smaller than protease-inhibitor systems) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation | < 0.1 | < 0.1 | CPU-bound; negligible |
| Equilibration (NVT + NPT) | 0.3 | 0.2 | 1.5 ns total; smaller system |
| Production MD (100 ns) | 2 | 1 | ~1200 ns/day (A100) for ~25k atoms |
| SMD (50 replicates × 3 ns) | 4 | 2 | 150 ns total |
| Umbrella Sampling (51 windows × 10.2 ns) | 40 | 20 | 520 ns total; sequential windows |
| WHAM/MBAR analysis | 0.5 | 0.5 | CPU-bound post-processing |
| **Total** | **~52** | **~26** | §10.17: Tier 1 ≈ 50 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each US window; after every 10 SMD replicates; after each equilibration phase |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-29/` |
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

EXP-29 has no upstream dependencies. EXP-29 results contribute to:

| Downstream Experiment | Usage |
|-----------------------|-------|
| EXP-26 (BSA–ΔG correlation) | ΔG data point (out-of-family validation) |

All outputs must be persisted to Google Drive immediately after generation.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
