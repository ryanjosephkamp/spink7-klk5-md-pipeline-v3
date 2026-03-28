# EXP-13: Barnase–Barstar Association Kinetics (kon)

**Experiment ID:** EXP-13  
**Feature ID:** F-13 (benchmarks.md)  
**Category:** Kinetic  
**Status:** SEMI-QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates the pipeline's ability to estimate ultra-fast, electrostatically-steered association kinetics using the barnase-barstar complex (PDB 1BRS). The experimental kon = 3.7 × 10⁸ M⁻¹s⁻¹ (Schreiber & Fersht 1996) is among the fastest known protein-protein associations, driven by long-range electrostatic complementarity. The pipeline should reproduce this rate within 1–2 orders of magnitude using PMF-based diffusion theory and/or Brownian dynamics estimates, and correctly identify the electrostatic steering mechanism.

---

## 2. Hypothesis

**H₁:** The PMF-derived or BD-estimated kon for barnase-barstar will fall within 10⁷–10¹⁰ M⁻¹s⁻¹ (within ~1 order of 3.7 × 10⁸).

**H₂:** Removing electrostatic interactions (charge-neutralized control) will reduce the estimated kon by ≥2 orders of magnitude, confirming electrostatic steering.

**H₃:** The electrostatic interaction energy at long range (ξ > 3 nm) will be attractive and > 2 kT, consistent with steered diffusion.

---

## 3. Background and Rationale

Barnase-barstar is the canonical example of electrostatic steering in protein-protein association. The large complementary charge patches (barnase is +2 near the active site, barstar is −3 near its binding interface) create a long-range funnel that accelerates association beyond the geometric collision limit. This experiment serves as a critical methodological control: if the pipeline correctly predicts fast, electrostatically-steered kinetics for barnase-barstar, it validates the PMF-based kinetic methodology used for the slower BPTI-trypsin system (EXP-09).

---

## 4. Experimental Protocol

### 4.1 System Preparation

- PDB: 1BRS (barnase-barstar complex, 2.0 Å resolution)
- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P (`amber14/tip3p.xml`)
- pH: 7.4 (SystemConfig default)
- Box padding: 1.2 nm, cubic box
- Ionic strength: 0.15 M NaCl

### 4.2 Equilibration

- Minimization: 10,000 steps, tolerance 10.0 kJ/mol/nm
- NVT: 500 ps at 310 K, friction 1.0 ps⁻¹, restraints 1000 kJ/mol/nm²
- NPT: 1000 ps at 310 K, 1.0 atm, barostat interval 25 steps
- Timestep: 2 fs

### 4.3 Umbrella Sampling for PMF

- Reaction coordinate ξ: COM distance, barnase ↔ barstar
- ξ range: 1.5–4.0 nm, spacing 0.05 nm (51 windows)
- Spring constant: 1000 kJ/mol/nm²
- Per-window: 10 ns (after 200 ps equilibration)
- Save interval: 1.0 ps
- WHAM: tolerance 10⁻⁶, max_iterations 100,000, bins 200
- Bootstrap: n_bootstrap = 200

### 4.4 Rate Constant Estimation

**Method A: Diffusion-influenced rate (Szabo-Shoup-Northrup):**

$$k_{on} = 4\pi D_{rel} R^* \exp\left(-\frac{W(R^*)}{k_BT}\right)$$

Where:
- D_rel = mutual diffusion coefficient (~15 Å²/ns for proteins this size)
- R* = reactive distance (contact radius)
- W(R*) = PMF at contact (if attractive, enhances rate)

**Method B: Full PMF integration (Smoluchowski with potential):**

$$\frac{1}{k_{on}} = \int_{R^*}^{\infty} \frac{\exp(W(r)/k_BT)}{4\pi r^2 D(r)} dr$$

### 4.5 Electrostatic Steering Control

Repeat the PMF calculation with all protein charges set to zero (neutral proteins):
- Same ξ range, windows, and sampling
- Expected: PMF flat at long range → no electrostatic steering
- Expected kon(neutral) << kon(charged), reduction ≥ 2 orders of magnitude

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | Pipeline kon |
|---------------|-------------|
| **PASS** | 10⁷–10¹⁰ M⁻¹s⁻¹ (within ~1 order) |
| **MARGINAL** | 10⁵–10⁶ or 10¹⁰–10¹¹ (within 2–3 orders) |
| **FAIL** | <10⁵ or >10¹¹ M⁻¹s⁻¹ |

---

## 5. Control Conditions

### 5.1 Positive Control

**EXP-09 (BPTI-trypsin):** kon ~10⁶ M⁻¹s⁻¹ (2 orders slower than barnase-barstar). The pipeline should correctly order: kon(barnase-barstar) >> kon(BPTI-trypsin).

### 5.2 Electrostatic Steering Control

**Charge-neutralized barnase-barstar:** As described in §4.5, removing charges should reduce kon by ≥2 orders, confirming the electrostatic steering mechanism.

### 5.3 Ionic Strength Control

At high ionic strength (1.0 M NaCl), electrostatic screening should reduce kon. Expected: kon(0.15 M) > kon(1.0 M), with reduction consistent with Debye-Hückel screening.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| kon (barnase-barstar) | 3.7 × 10⁸ M⁻¹s⁻¹ | Schreiber & Fersht 1996 |
| Pipeline target | 10⁷–10¹⁰ M⁻¹s⁻¹ | benchmarks.md |
| kon(neutral) / kon(charged) | < 0.01 | Electrostatic steering |
| Long-range PMF | Attractive (> 2 kT at 3 nm) | Electrostatic funnel |
| Electrostatic enhancement | ≥100× over neutral collision rate | Gabdoulline & Wade 1997 |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Ionic screening overestimated** | kon too low at 0.15 M | PME artifacts at boundary | Medium |
| **D_rel uncertainty** | Diffusion coefficient varies by 2× | Hydrodynamic interactions | Medium |
| **1D PMF misses angular steering** | No angular dependence captured | Scalar ξ insufficient | Medium |
| **Charge neutralization too drastic** | Protein unfolds without charges | Non-physical control | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Structure quality | 1BRS loads without missing atoms at interface | All heavy atoms present |
| Equilibration | RMSD < 2.0 Å, no interface disruption | Stable complex |
| PMF shape | Attractive well + plateau at large ξ | Physically reasonable |
| Long-range PMF | Attractive at ξ > 3 nm (electrostatic funnel) | W(ξ) < 0 at 3 nm |
| Rank order | kon(barnase-barstar) > kon(BPTI-trypsin) | Correct ranking |
| Electrostatic test | Neutral protein PMF flat at long range | |W(3 nm)| < 0.5 kT |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (barnase-barstar complex: ~50,000 atoms solvated in TIP3P with 1.2 nm padding) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation | < 0.1 | < 0.1 | CPU-bound; negligible |
| Equilibration (NVT + NPT) | 0.5 | 0.3 | 1.5 ns total |
| Umbrella Sampling — charged (51 windows × 10.2 ns) | 18 | 9 | 520 ns total |
| Umbrella Sampling — neutral control (51 windows × 10.2 ns) | 18 | 9 | 520 ns charge-neutralized repeat |
| Rate constant estimation (Szabo-Shoup + Smoluchowski) | 0.5 | 0.5 | CPU-bound analysis |
| Ionic strength control (1.0 M NaCl) | 5 | 2.5 | Partial US re-run at high salt |
| **Total** | **~42** | **~21** | §10.17: Tier 1 ≈ 42 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each US window; after each equilibration phase |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-13/` |
| Restart procedure | Mount Drive → load latest checkpoint → verify energy drift < 0.1% → resume from next incomplete phase |
| Estimated sessions needed | 2–3 (A100) or 1–2 (H100) |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Positions, velocities, box vectors, RNG state | OpenMM binary `.chk` | `checkpoint_<phase>_<step>.chk` |
| Full simulation state (portable) | OpenMM XML serialization | `state_<phase>_<step>.xml` |
| Umbrella window ξ timeseries | NumPy `.npy` | `umbrella_window_<condition>_<idx>.npy` |
| PMF results (charged, neutral, high-salt) | NumPy `.npz` | `pmf_<condition>.npz` |

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

EXP-13 has no upstream dependencies. EXP-13 results contribute to:

| Downstream Experiment | Usage |
|-----------------------|-------|
| EXP-33 (Barnase-barstar DMC) | Uses same 1BRS system; shares equilibrated structure |
| EXP-26 (BSA–ΔG correlation) | Kinetic data point for pipeline validation |

All outputs must be persisted to Google Drive immediately after generation.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
