# EXP-26: BSA–ΔG Binding Correlation

**Experiment ID:** EXP-26  
**Feature ID:** F-26 (benchmarks.md)  
**Category:** Biophysical  
**Status:** SEMI-QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment tests whether the pipeline reproduces the empirical correlation between buried surface area (BSA) and binding free energy (ΔG_bind) observed in protein-protein complexes. The established correlation predicts ΔG ≈ −0.01 × BSA (kcal/mol per Å²) (Horton & Lewis 1992), or approximately −15 kcal/mol for the ~1530 Å² BPTI-trypsin interface. The pipeline should demonstrate this correlation across its benchmark systems (EXP-01 through EXP-06), validating internal thermodynamic consistency.

---

## 2. Hypothesis

**H₁:** The correlation between BSA (EXP-16) and ΔG_bind (EXP-01 through EXP-06) across benchmark systems will show r² > 0.5 (moderate-to-strong correlation).

**H₂:** The slope of the BSA-ΔG regression will be in the range −0.005 to −0.020 kcal/mol/Å² (literature: ~−0.01).

---

## 3. Background and Rationale

The BSA-ΔG correlation is a well-known empirical relationship in protein-protein interactions. While the correlation is imperfect (R² ~ 0.3–0.7 across diverse complexes), it provides a useful sanity check: larger interfaces should generally bind more tightly. A pipeline that produces ΔG values uncorrelated with BSA likely has systematic errors in either structural or thermodynamic calculations. This experiment is a meta-analysis of EXP-01 through EXP-06 and EXP-16 results.

---

## 4. Experimental Protocol

### 4.1 Data Collection

Compile from completed experiments:

| System | EXP (BSA) | EXP (ΔG) | Expected BSA (Å²) | Expected ΔG (kcal/mol) |
|--------|-----------|-----------|-------------------|------------------------|
| BPTI-trypsin | EXP-16/04 | EXP-04 | 1530 | −18.0 |
| SPINK7-KLK5 | EXP-16/01 | EXP-01 | ~1300 | −9.4 |
| SPINK7-KLK12 | EXP-16/02 | EXP-02 | ~1200 | −7.9 |
| LEKTI-KLK5 | EXP-16/03 | EXP-03 | ~1400 | −11.8 |
| PSTI-chymotrypsinogen | EXP-16/05 | EXP-05 | ~1450 | −14.7 |
| SPINK1-trypsin | EXP-16/06 | EXP-06 | ~1350 | −11.1 |

### 4.2 Correlation Analysis

1. Plot BSA (x-axis) vs. ΔG_bind (y-axis) for all 6 systems.
2. Linear regression: ΔG = m × BSA + b.
3. Report: slope (m), intercept (b), R², p-value.
4. Compare slope to literature value (−0.01 kcal/mol/Å²).

### 4.3 Acceptance Criteria (from benchmarks.md)

| Classification | Criterion |
|---------------|-----------|
| **PASS** | R² > 0.5 AND slope in [−0.005, −0.020] |
| **MARGINAL** | R² > 0.3 OR correct trend (negative slope) |
| **FAIL** | No correlation (R² < 0.1) or positive slope |

---

## 5. Control Conditions

### 5.1 Literature Data Points

Add literature BSA-ΔG data points from non-pipeline systems (e.g., barnase-barstar, SH3-p41 from EXP-29) to increase statistical power.

### 5.2 Rank-Order Test

Even without precise absolute values, the rank order should be preserved: tighter binders should have larger BSA. Spearman rank correlation ρ > 0.5.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Slope | ~−0.01 kcal/mol/Å² | Horton & Lewis 1992 |
| R² | 0.3–0.7 | Literature range |
| Rank order | BPTI > PSTI > LEKTI > SPINK1 > SPINK7 > SPINK7-KLK12 | By ΔG magnitude |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Small sample size** | n=6 too few for robust regression | Limited benchmark set | Medium |
| **BSA measurement error** | Distorts correlation | Different complex geometries | Low |
| **ΔG method inconsistency** | Different methods for different systems | SMD vs US vs FEP | Medium |
| **Non-linear relationship** | Log-linear or saturating relationship | Linear model insufficient | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Data completeness | BSA and ΔG available for all 6 systems | No missing values |
| Negative slope | ΔG decreases (more negative) with increasing BSA | Correct trend |
| No outliers | No system > 3σ from regression line | Physical consistency |
| Rank order preserved | Spearman ρ > 0.5 | Correct ordering |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | Minimal (~1 GB) — pure post-processing; no live simulation |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Collect BSA values from Tier 1 equilibrated structures | 0.1 | 0.1 | MDTraj SASA calculation |
| Collect ΔG values from Tier 1 results | < 0.1 | < 0.1 | File I/O from Drive |
| Linear regression + R² | < 0.1 | < 0.1 | CPU-bound statistics |
| Bootstrap analysis (1000 resamples) | 0.1 | 0.1 | Uncertainty on R² |
| Visualization | < 0.1 | < 0.1 | Matplotlib scatter + regression |
| **Total** | **~0.5** | **~0.5** | §10.17: Tier 4 ≈ 0.5 GPU-hrs (pure analysis) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | Not applicable (single short analysis) |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-26/` |
| Restart procedure | Re-run analysis from start (< 30 min total) |
| Estimated sessions needed | 1 (can share session with other Tier 3/4 analyses) |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| BSA-ΔG data pairs | CSV | `bsa_dg_data.csv` |
| Regression results (R², slope, intercept, CI) | JSON | `correlation_results.json` |
| Scatter plot | PNG | `bsa_dg_correlation.png` |

### 9.5 Platform Configuration

No GPU simulation required. Analysis uses MDTraj (for SASA), NumPy, SciPy, and Matplotlib (CPU).

```python
import mdtraj as md
import numpy as np
from scipy import stats
# SASA calculation on equilibrated structures from each Tier 1 system
```

### 9.6 Dependency Notes

**CRITICAL — Tier 4: requires completed Tier 1 results.**

| Required Upstream Output | Source Experiment | Drive Path |
|-------------------------|-------------------|------------|
| ΔG_bind (BPTI-trypsin) | EXP-04 | `/content/drive/MyDrive/v3_gpu_results/EXP-04/analysis/` |
| ΔG_bind (PSTI-chymotrypsinogen) | EXP-05 | `/content/drive/MyDrive/v3_gpu_results/EXP-05/analysis/` |
| ΔG_bind (SPINK1-trypsin) | EXP-06 | `/content/drive/MyDrive/v3_gpu_results/EXP-06/analysis/` |
| ΔG_bind (SH3-p41) | EXP-29 | `/content/drive/MyDrive/v3_gpu_results/EXP-29/analysis/` |
| Equilibrated complex structures (for BSA) | EXP-04, 05, 06, 29 | Production PDB files from each |

**Pre-execution check:** ALL Tier 1 experiments must be completed before running EXP-26. At minimum 3 data points are required for a meaningful correlation.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
