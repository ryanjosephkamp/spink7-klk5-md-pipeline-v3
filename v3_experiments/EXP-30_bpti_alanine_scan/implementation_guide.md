# EXP-30: BPTI Reactive Region Alanine Scan ΔΔG Panel — Implementation Guide

**Experiment ID:** EXP-30  
**Feature ID:** F-30 (benchmarks.md)  
**Category:** Mutational (Quantitative — Alchemical FEP)  
**Date:** 2026-03-24  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Computes ΔΔG for 15 alanine mutations in the BPTI–trypsin reactive region (T11A through G36A) using Free Energy Perturbation (FEP). Experimental data from Castro & Anderson (1996) provides the gold-standard benchmark with full kinetic (k_on, k_off, Ki) and thermodynamic (ΔG, ΔH, TΔS) characterization. Key mutation: K15A (P1) with ΔΔG ≈ +10 kcal/mol. Acceptance: Spearman ρ > 0.7, RMSE < 2.0 kcal/mol, ≥80% correct sign.

### 2. Hypotheses

- **H₁:** FEP-computed ΔΔG values will achieve a Spearman rank correlation ρ > 0.7 against experimental values.
- **H₂:** Per-mutant |ΔΔG_comp − ΔΔG_exp| < 2.3 kcal/mol for ≥ 12 of 15 mutations (80%).
- **H₃:** K15A (P1) will be identified as the most destabilizing mutation, with computed ΔΔG > +7 kcal/mol.
- **H₄:** G36A will be identified as neutral, with computed |ΔΔG| < 1.5 kcal/mol.

### 3. Classification (§25.1)

- **PASS:** Spearman ρ > 0.7 AND RMSE < 2.0 kcal/mol AND ≥80% correct ΔΔG sign
- **MARGINAL:** Spearman ρ > 0.5 AND RMSE < 3.5 kcal/mol AND ≥60% correct sign
- **FAIL:** ρ < 0.5 OR RMSE > 3.5 OR <60% correct sign

### 4. Mutation Panel

From the canonical `experimental_design.md` and Castro & Anderson (1996):

| # | Mutant | Experimental ΔΔG (kcal/mol) | Ki (M) | Category |
|---|--------|----------------------------|--------|----------|
| 1 | T11A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 2 | G12A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 3 | P13A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 4 | C14A | ≈ +7 | ~10⁻⁹ | Disulfide (Cys14-Cys38) |
| 5 | K15A (P1) | ≈ +10 | 1.4 × 10⁻⁶ | P1 dominant |
| 6 | A16G | varies | — | Framework control |
| 7 | R17A (P2') | ≈ +5 | ~10⁻¹⁰ | Interface contact |
| 8 | I18A | ≈ +3–4 | ~10⁻¹¹ | Hydrophobic packing |
| 9 | I19A | ≈ +3–4 | ~10⁻¹¹ | Hydrophobic packing |
| 10 | R20A | ≈ +2–3 | ~10⁻¹² | Salt bridge |
| 11 | Y21A | ≈ +2–3 | ~10⁻¹² | Aromatic contact |
| 12 | F22A | ≈ +2–3 | ~10⁻¹¹ | Hydrophobic packing |
| 13 | Y23A | ≈ +1–2 | ~10⁻¹³ | Peripheral aromatic |
| 14 | C38A | ≈ +7 | ~10⁻⁹ | Disulfide (Cys14-Cys38) |
| 15 | G36A | ≈ 0 | ~10⁻¹⁴ | Null control |

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json, time
import numpy as np
from pathlib import Path
from scipy import stats

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-30_bpti_alanine_scan"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Correct imports — FEP analysis from src.analyze.fep (NOT src.analyze.mbar)
from src.simulate.fep import run_fep_campaign
from src.analyze.fep import compute_delta_g_mbar, compute_delta_delta_g
from src.config import FEPConfig, KCAL_TO_KJ, BOLTZMANN_KJ

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("All imports successful.")
```

### Step 2: Configure FEP Parameters

```python
TEMPERATURE_K = 310.0

# Optimized: 16 lambda windows (from default 20)
# Justification: Shirts & Chodera (2008) — MBAR achieves optimal efficiency
# at 12–20 windows for single-residue X→Ala mutations
fep_config = FEPConfig(
    n_lambda_windows=16,
    per_window_duration_ns=2.0,
    temperature_k=TEMPERATURE_K,
)
```

### Step 3: Define Mutation Panel

```python
mutations = [
    {"id": 1,  "name": "T11A", "resid": 11, "ddg_exp": 1.5},
    {"id": 2,  "name": "G12A", "resid": 12, "ddg_exp": 1.0},
    {"id": 3,  "name": "P13A", "resid": 13, "ddg_exp": 1.5},
    {"id": 4,  "name": "C14A", "resid": 14, "ddg_exp": 7.0},
    {"id": 5,  "name": "K15A", "resid": 15, "ddg_exp": 10.0},
    {"id": 6,  "name": "A16G", "resid": 16, "ddg_exp": 0.0},   # Framework control
    {"id": 7,  "name": "R17A", "resid": 17, "ddg_exp": 5.0},
    {"id": 8,  "name": "I18A", "resid": 18, "ddg_exp": 3.5},
    {"id": 9,  "name": "I19A", "resid": 19, "ddg_exp": 3.5},
    {"id": 10, "name": "R20A", "resid": 20, "ddg_exp": 2.5},
    {"id": 11, "name": "Y21A", "resid": 21, "ddg_exp": 2.5},
    {"id": 12, "name": "F22A", "resid": 22, "ddg_exp": 2.5},
    {"id": 13, "name": "Y23A", "resid": 23, "ddg_exp": 1.5},
    {"id": 14, "name": "C38A", "resid": 38, "ddg_exp": 7.0},
    {"id": 15, "name": "G36A", "resid": 36, "ddg_exp": 0.0},   # Null control
]
```

### Step 4: System Preparation

Load pre-equilibrated BPTI-trypsin complex from EXP-04 if available, otherwise prepare from PDB 2PTC:

```python
# Complex system: BPTI-trypsin from PDB 2PTC
# If EXP-04 has been run, load equilibrated state
# Otherwise, run full prep: fetch → clean → topology → solvate → minimize → equilibrate

from src.prep.pdb_fetch import fetch_pdb
from src.prep.pdb_clean import clean_pdb
from src.prep.topology import build_topology
from src.prep.solvate import solvate_system
from src.simulate.minimizer import minimize_energy
from src.simulate.equilibrate import run_nvt, run_npt
from src.config import SystemConfig, MinimizationConfig, EquilibrationConfig
```

### Step 5: FEP Execution (Per Mutation, Two Legs)

For each mutation, the FEP thermodynamic cycle requires two legs:

```
ΔΔG_bind = ΔG_complex(WT→mut) − ΔG_free(WT→mut)
```

```python
for mut in mutations:
    # Identify sidechain atom indices beyond backbone (atoms to annihilate)
    sc_complex = get_sidechain_indices(complex_topology, mut["resid"], chain_id=bpti_chain)
    sc_free = get_sidechain_indices(free_topology, mut["resid"], chain_id=0)

    # Complex leg
    fep_result_complex = run_fep_campaign(
        system=complex_system,
        positions=complex_positions,
        mutant_atom_indices=sc_complex,
        config=fep_config,
        output_dir=output_dir / "complex" / mut["name"]
    )
    dg_complex = compute_delta_g_mbar(
        fep_result_complex['energy_matrix'],
        fep_result_complex['n_samples_per_state'],
        TEMPERATURE_K
    )

    # Free leg
    fep_result_free = run_fep_campaign(
        system=free_system,
        positions=free_positions,
        mutant_atom_indices=sc_free,
        config=fep_config,
        output_dir=output_dir / "free" / mut["name"]
    )
    dg_free = compute_delta_g_mbar(
        fep_result_free['energy_matrix'],
        fep_result_free['n_samples_per_state'],
        TEMPERATURE_K
    )

    # ΔΔG from thermodynamic cycle
    ddg = compute_delta_delta_g(dg_complex, dg_free)
    # ddg contains: delta_delta_g_kcal_mol, delta_delta_g_std_kcal_mol
```

### Step 6: Statistical Analysis

```python
from scipy.stats import spearmanr, pearsonr

sim_ddg = np.array([results[m]["ddg_sim"] for m in mutation_names])
exp_ddg = np.array([results[m]["ddg_exp"] for m in mutation_names])

rho, p_rho = spearmanr(exp_ddg, sim_ddg)
r, p_r = pearsonr(exp_ddg, sim_ddg)
rmse = np.sqrt(np.mean((sim_ddg - exp_ddg) ** 2))
sign_correct = np.sum(np.sign(sim_ddg) == np.sign(exp_ddg))
sign_pct = sign_correct / len(sim_ddg) * 100

# Classification
if rho > 0.7 and rmse < 2.0 and sign_pct >= 80:
    verdict = "PASS"
elif rho > 0.5 and rmse < 3.5 and sign_pct >= 60:
    verdict = "MARGINAL"
else:
    verdict = "FAIL"
```

### Step 7: Figure Generation

**Figure 1:** Scatter plot of computed vs. experimental ΔΔG with error bars, identity line, and linear fit.

**Figure 2:** Bar chart showing per-mutation ΔΔG with experimental values overlaid.

**Figure 3:** Heatmap of FEP convergence (forward/reverse overlap per mutation).

---

## Part 3 — Expected Outputs

| File | Description |
|------|-------------|
| `outputs/results.json` | Complete results dictionary with per-mutation ΔΔG |
| `outputs/checkpoint_ddg.json` | Incremental checkpoint (mutation-by-mutation) |
| `figures/ddg_scatter.png` | Computed vs. experimental ΔΔG scatter |
| `figures/ddg_bar_chart.png` | Per-mutation bar chart |
| `figures/fep_convergence.png` | Convergence heatmap |

---

## Part 4 — Verification Criteria

| Criterion | Target | Method |
|-----------|--------|--------|
| Spearman ρ | > 0.7 | `scipy.stats.spearmanr` |
| RMSE | < 2.0 kcal/mol | `np.sqrt(np.mean((sim - exp)²))` |
| Sign agreement | ≥ 80% | Count correct signs / total |
| K15A is most destabilizing | ΔΔG > 7 | Direct comparison |
| G36A is neutral | \|ΔΔG\| < 1.5 | Direct comparison |
| All ΔΔG finite | No NaN or Inf | `np.isfinite()` check |

---

## Part 5 — GPU/Colab Execution

### 5.1 Resource Requirements

| Parameter | Value |
|-----------|-------|
| GPU | NVIDIA A100 or H100 (CUDA) |
| Estimated GPU-hours | ~82h (A100), ~55h (H100) |
| RAM | ≥ 40 GB |
| Disk | ≥ 20 GB |
| Lambda windows | 16 (optimized from default 20) |
| Total simulation | 15 mutations × 16 windows × 2 ns × 2 legs = 960 ns |

### 5.2 Colab Notebook

The experiment is implemented in `EXP-30_colab.ipynb` in this directory. The notebook:

1. Installs OpenMM and dependencies
2. Mounts Google Drive and loads pipeline source from `MyDrive/spink7_pipeline/`
3. Verifies CUDA availability
4. Loads or prepares BPTI-trypsin complex and free BPTI systems
5. Iterates over 15 mutations with checkpoint/resume support
6. Performs correlation analysis and classification
7. Saves results and figures to Google Drive

### 5.3 Checkpoint/Resume

The notebook saves `checkpoint_ddg.json` to Google Drive after each mutation, enabling resume across Colab sessions. This is critical given the ~82h runtime requires 2–3 sessions.

### 5.4 Dependencies

- **EXP-04:** Optionally loads pre-equilibrated complex from `v3_gpu_results/EXP-04/`. Falls back to independent preparation if EXP-04 not yet completed.

---

## References

1. Castro, M. J. & Anderson, S. (1996). Alanine point-mutations in the reactive region of bovine pancreatic trypsin inhibitor. *Biochemistry*, 35(35), 11435–11446.
2. Shirts, M. R. & Chodera, J. D. (2008). Statistically optimal analysis of samples from multiple equilibrium states. *Journal of Chemical Physics*, 129, 124105.
3. Krowarsch, D., Zakrzewska, M., Skowron, P., & Otlewski, J. (2003). Structure–function relationships in serine protease–inhibitor interactions. *Acta Biochimica Polonica*, 50(2), 367–381.

---

**Author:** Ryan Kamp  
**Affiliation:** Dept. of Computer Science, University of Cincinnati  
**Email:** kamprj@mail.uc.edu  
**GitHub:** ryanjosephkamp
