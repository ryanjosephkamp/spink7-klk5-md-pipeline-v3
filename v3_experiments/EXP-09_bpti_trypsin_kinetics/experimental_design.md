# EXP-09: BPTI–Trypsin Association/Dissociation Kinetics (kon/koff)

**Experiment ID:** EXP-09  
**Feature ID:** F-09 (benchmarks.md)  
**Category:** Kinetic  
**Status:** SEMI-QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment tests whether the V2 pipeline can produce kinetic estimates for BPTI-trypsin association and dissociation that fall within the expected order-of-magnitude ranges. The experimental reference values are kon ≈ 10⁶ M⁻¹s⁻¹ and koff ≈ 6 × 10⁻⁸ s⁻¹ (Kd = 6 × 10⁻¹⁴ M). Since direct kinetic simulation of events with timescales >10⁷ s is computationally intractable, the approach uses PMF-based transition-state theory and/or MSM-derived rate estimates to infer kinetics from the free-energy landscape.

---

## 2. Hypothesis

**H₁:** The PMF-derived kon estimate will fall within 1–2 orders of magnitude of the experimental value (10⁴–10⁸ M⁻¹s⁻¹).

**H₂:** The PMF-derived koff estimate will fall within 1–2 orders of magnitude of the experimental value (10⁻¹⁰–10⁻⁶ s⁻¹).

**H₃:** The ratio koff/kon will yield a Kd consistent with the thermodynamic ΔG from EXP-04 (within 2 kcal/mol via ΔG = −RT ln Kd).

---

## 3. Background and Rationale

BPTI-trypsin is one of the tightest known protein-protein interactions (Kd ≈ 6 × 10⁻¹⁴ M). The fast kon (~10⁶ M⁻¹s⁻¹) indicates partially diffusion-controlled association, while the extremely slow koff (~10⁻⁸ s⁻¹, t₁/₂ ≈ 200 days) reflects the deep free-energy minimum. PMF-derived rate constants use Kramers' theory or transition-state theory (TST) applied to the 1D PMF along the unbinding coordinate ξ, providing mechanistic links between the barrier height (EXP-10) and kinetics.

---

## 4. Experimental Protocol

### 4.1 PMF from Umbrella Sampling

Use the PMF computed in EXP-04 (US/WHAM along COM distance ξ).
- ξ range: 1.5–4.0 nm, 51 windows, spacing 0.05 nm
- Per-window sampling: 10 ns
- Spring constant: 1000 kJ/mol/nm²
- WHAM tolerance: 10⁻⁶, max iterations: 100,000, bins: 200
- Bootstrap uncertainty: n_bootstrap = 200

### 4.2 Rate Constant Estimation

**Method A: Transition-state theory (TST)**

$$k_{off} = \frac{D(\xi^{\ddagger})}{2\pi k_BT} \cdot \sqrt{\frac{|W''(\xi_{min})| \cdot |W''(\xi^{\ddagger})|}{(k_BT)^2}} \cdot \exp\left(-\frac{\Delta W^{\ddagger}}{k_BT}\right)$$

Where:
- ξ^‡ = transition state (PMF maximum along ξ)
- ΔW^‡ = barrier height (EXP-10: ~10.5 kcal/mol)
- D(ξ) = position-dependent diffusion coefficient from autocorrelation of ξ

**Method B: Kramers' theory**

$$k_{off} = \frac{\omega_{min} \omega_{max}}{2\pi \gamma} \cdot \exp\left(-\frac{\Delta W^{\ddagger}}{k_BT}\right)$$

With friction coefficient γ from Langevin dynamics.

**kon from detailed balance:**

$$k_{on} = \frac{k_{off}}{K_d} = k_{off} \cdot C^{\circ} \cdot \exp\left(\frac{\Delta G_{bind}}{k_BT}\right)$$

### 4.3 MSM-Based Approach (if applicable)

From production trajectories:
- tICA lag time: 10 ps, 4 components
- k-means clusters: 200
- Bayesian MSM samples: 100
- Lag times: (1, 2, 5, 10, 20, 50, 100, 200, 500) ps
- Extract slowest implied timescales; if unbinding transition is captured, extract kon/koff directly.

### 4.4 Acceptance Criteria (from benchmarks.md)

| Classification | Criterion |
|---------------|-----------|
| **PASS** | Both kon and koff within 2 orders of magnitude of experiment |
| **MARGINAL** | One rate within 2 orders, other within 3 |
| **FAIL** | Either rate off by >3 orders of magnitude |

---

## 5. Control Conditions

### 5.1 Positive Control

**Barnase-barstar (EXP-13):** kon = 3.7 × 10⁸ M⁻¹s⁻¹ — a fast, electrostatically steered association. If the PMF-based method works for barnase-barstar, it validates the methodology.

### 5.2 Consistency Control

**Thermodynamic consistency:** Kd(kinetic) = koff/kon must agree with Kd(thermodynamic) from EXP-04 within 2 kcal/mol (via ΔG = −RT ln Kd).

### 5.3 Negative Control

**Diffusion limit check:** kon must not exceed the Smoluchowski diffusion limit (~10⁹–10¹⁰ M⁻¹s⁻¹). Values exceeding this indicate methodology errors.

---

## 6. Expected Outcomes

| Metric | Experimental Value | Pipeline Target | Tolerance |
|--------|-------------------|-----------------|-----------|
| kon | ~10⁶ M⁻¹s⁻¹ | 10⁴–10⁸ M⁻¹s⁻¹ | ±2 orders |
| koff | ~6 × 10⁻⁸ s⁻¹ | 10⁻¹⁰–10⁻⁶ s⁻¹ | ±2 orders |
| Kd(kinetic) | 6 × 10⁻¹⁴ M | Within 2 kcal/mol of EXP-04 ΔG | Consistency |
| PMF barrier | ~10.5 kcal/mol | See EXP-10 | Cross-ref |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **1D reaction coordinate inadequacy** | PMF barrier does not capture true TS | Multi-dimensional unbinding path | High |
| **Diffusion coefficient uncertainty** | D(ξ) varies by orders of magnitude | Autocorrelation convergence | High |
| **MSM insufficient sampling** | Unbinding timescale > total simulation | 100 ns too short for koff ≈ 10⁻⁸ s⁻¹ | High |
| **Kramers vs TST disagreement** | Methods differ by >2 orders | Friction regime uncertainty | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| PMF quality | Smooth PMF with clear minimum and barrier | No discontinuities, bootstrap error < 2 kcal/mol |
| Barrier identification | Single, well-defined TS along ξ | Clear maximum between bound and unbound |
| D(ξ) convergence | Block averaging of autocorrelation | CV < 50% across blocks |
| Detailed balance | koff/kon consistent with Kd | Within 2 kcal/mol of EXP-04 |
| Kramers vs TST | Both methods agree within 1 order | Validates rate estimate |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~2–4 GB (PMF-based rate constant calculation; primarily CPU-bound analysis on pre-computed PMF data) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Load EXP-04 PMF (ξ, W(ξ), bootstrap samples) | < 0.1 | < 0.1 | I/O from Drive |
| Diffusion coefficient estimation | 0.5 | 0.3 | Short MD runs for D(r) if needed |
| Kramers/TST rate constant calculation | 0.5 | 0.5 | CPU-bound numerical integration |
| Smoluchowski integration | 0.5 | 0.5 | Full PMF integration |
| Bootstrap uncertainty propagation | 0.5 | 0.5 | Monte Carlo over PMF bootstrap samples |
| **Total** | **~2–3** | **~1.5–2** | §10.17: Tier 2 ≈ 2–3 GPU-hrs |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After diffusion coefficient estimation; after each rate method calculation |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-09/` |
| Restart procedure | Mount Drive → load partial kinetics results → resume from next incomplete method |
| Estimated sessions needed | 1 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Diffusion coefficient D(r) | NumPy `.npy` | `diffusion_coefficient.npy` |
| Rate constants (Kramers, TST, Smoluchowski) | JSON | `rate_constants.json` |
| Bootstrap rate distributions | NumPy `.npz` | `rate_bootstrap_distributions.npz` |

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
| PMF from US/WHAM (ξ, W(ξ)) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/pmf_wham.npz` | Full PMF profile |
| PMF bootstrap samples (200 replicates) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/pmf_bootstrap.npz` | For uncertainty propagation |
| Equilibrated complex structure | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/equilibrated_complex.pdb` | For contact radius R* determination |

**Pre-execution check:** Verify all EXP-04 PMF outputs exist on Drive before starting EXP-09. If missing, EXP-04 must be completed first (including full US/WHAM analysis).

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
