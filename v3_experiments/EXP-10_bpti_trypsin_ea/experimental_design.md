# EXP-10: BPTI–Trypsin Activation Energy (Ea) from PMF Barrier Height

**Experiment ID:** EXP-10  
**Feature ID:** F-10 (benchmarks.md)  
**Category:** Kinetic  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment measures the activation energy for BPTI-trypsin dissociation from the height of the free-energy barrier in the PMF along the unbinding coordinate. The experimental Arrhenius activation energy is Ea = 10.5 kcal/mol (Quast et al. 1974, Vincent & Bhatt 2007), with a pipeline 95% CI of [4.6, 16.4] kcal/mol. The PMF barrier height from umbrella sampling directly estimates this quantity, providing a key validation of the energy landscape along the dissociation pathway.

---

## 2. Hypothesis

**H₁:** The PMF barrier height (ΔW‡ = W(ξ‡) − W(ξ_min)) along the COM unbinding coordinate is 10.5 ± 5.9 kcal/mol (95% CI: [4.6, 16.4] kcal/mol).

**H₂:** The barrier height is consistent with the kinetic koff from EXP-09 via Kramers'/TST relationship.

---

## 3. Background and Rationale

The Arrhenius activation energy for BPTI-trypsin dissociation (10.5 kcal/mol) was measured via temperature-dependent kinetics. In the context of a 1D PMF along the unbinding coordinate ξ, this corresponds to the barrier height between the bound-state minimum and the transition state. A correct barrier height validates both the thermodynamic landscape (EXP-04) and the kinetic predictions (EXP-09) — it is the key link between equilibrium and kinetic properties.

---

## 4. Experimental Protocol

### 4.1 PMF from Umbrella Sampling

Use the PMF computed in EXP-04 (shared data dependency):
- Reaction coordinate ξ: COM distance between BPTI and trypsin
- ξ range: 1.5 nm to 4.0 nm
- Window spacing: 0.05 nm (51 windows)
- Spring constant: k = 1000 kJ/mol/nm²
- Per-window sampling: 10 ns
- Pre-positioning velocity: 0.01 nm/ps
- Equilibration per window: 200 ps (discarded)
- Equilibration detection: enabled
- Save interval: 1.0 ps
- WHAM: tolerance = 10⁻⁶, max_iterations = 100,000, bins = 200
- Bootstrap: n_bootstrap = 200

### 4.2 Barrier Height Extraction

1. Identify ξ_min: the global minimum of W(ξ) (bound state, expected ~2.0–2.5 nm).
2. Identify ξ‡: the global maximum of W(ξ) between ξ_min and the unbound plateau.
3. Compute ΔW‡ = W(ξ‡) − W(ξ_min).
4. Uncertainty from bootstrap: σ(ΔW‡) from the distribution of ΔW‡ across 200 bootstrap samples.

### 4.3 Cross-Validation with MBAR

Independently compute PMF using MBAR:
- Solver: "robust"
- Tolerance: 10⁻⁷
- Max iterations: 10,000
- n_bootstrap: 200
- PMF bins: 200

Compare barrier heights from WHAM and MBAR; they should agree within bootstrap uncertainty.

### 4.4 Acceptance Criteria (from benchmarks.md)

| Classification | Barrier Height (kcal/mol) |
|---------------|--------------------------|
| **PASS** | [4.6, 16.4] |
| **MARGINAL** | [2.0, 20.0] |
| **FAIL** | Outside marginal range or no clear barrier |

---

## 5. Control Conditions

### 5.1 Positive Control

**Flat PMF at large ξ:** Beyond the transition state, the PMF should plateau to a constant value (bulk solvent, no interaction). If the PMF does not plateau, the pulling distance is insufficient.

### 5.2 Consistency Controls

- **EXP-09 cross-check:** ΔW‡ must be consistent with koff via Kramers' theory.
- **EXP-04 cross-check:** The well depth W(ξ_unbound) − W(ξ_min) should approximate ΔG_bind after volume/entropy correction.

### 5.3 Negative Control

**Pure solvent run:** Pulling two non-interacting species (e.g., two water boxes) through the same ξ range should yield a flat PMF (ΔW‡ ≈ 0).

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Barrier height ΔW‡ | 10.5 kcal/mol | Quast 1974, Vincent & Bhatt 2007 |
| 95% CI | [4.6, 16.4] kcal/mol | benchmarks.md |
| ξ_min | ~2.0–2.5 nm | Structural (PDB 2PTC) |
| ξ‡ | ~2.8–3.2 nm | Expected from similar systems |
| WHAM-MBAR agreement | Within bootstrap σ | Internal consistency |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **No clear barrier** | Monotonic PMF (no maximum) | 1D coordinate misses true TS | Critical |
| **Multiple barriers** | Several local maxima | Multi-step dissociation mechanism | Medium |
| **Bootstrap error too large** | σ(ΔW‡) > 5 kcal/mol | Insufficient sampling per window | Medium |
| **WHAM-MBAR disagreement** | Barrier heights differ by >2 kcal/mol | Convergence issues | Medium |
| **ξ range insufficient** | PMF does not plateau at large ξ | Need larger pull_distance | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Umbrella overlap | Histogram overlap between adjacent windows | >10% overlap all pairs |
| PMF convergence | Block averaging (first/second half agree) | ΔW‡ within 2 kcal/mol between halves |
| Plateau check | W(ξ > 3.5 nm) approximately constant | Slope < 0.5 kcal/mol/nm |
| Bootstrap distribution | ΔW‡ bootstrap histogram approximately normal | Shapiro-Wilk p > 0.05 |
| WHAM vs MBAR | Independent barrier height estimates | Agreement within bootstrap σ |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~2–4 GB (PMF barrier analysis on pre-computed data; primarily CPU-bound) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Load EXP-04 PMF (ξ, W(ξ), bootstrap samples) | < 0.1 | < 0.1 | I/O from Drive |
| Barrier height extraction + fitting | 0.3 | 0.3 | CPU-bound curve analysis |
| Temperature-dependent PMF (optional multi-T) | 0.5 | 0.3 | Short re-equilibration if needed |
| Arrhenius/Eyring analysis | 0.2 | 0.2 | CPU-bound fitting |
| Bootstrap uncertainty propagation | 0.3 | 0.3 | Monte Carlo over PMF bootstrap samples |
| **Total** | **~1–2** | **~1** | §10.17: Tier 2 ≈ 1–2 GPU-hrs |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After barrier extraction; after each analysis step |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-10/` |
| Restart procedure | Mount Drive → load partial analysis → resume from next incomplete step |
| Estimated sessions needed | 1 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Barrier height + position | JSON | `barrier_analysis.json` |
| Activation energy (Arrhenius/Eyring fits) | JSON | `activation_energy.json` |
| Bootstrap distributions | NumPy `.npz` | `ea_bootstrap_distributions.npz` |

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
| PMF from US/WHAM (ξ, W(ξ)) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/pmf_wham.npz` | Full PMF profile for barrier identification |
| PMF bootstrap samples (200 replicates) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/pmf_bootstrap.npz` | For uncertainty propagation on E_a |
| Per-window equilibrated states (optional) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/umbrella/` | For optional multi-temperature re-runs |

**Pre-execution check:** Verify all EXP-04 PMF outputs exist on Drive before starting EXP-10. If missing, EXP-04 must be completed first (including full US/WHAM analysis).

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
