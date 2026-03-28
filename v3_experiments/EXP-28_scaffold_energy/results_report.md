# EXP-28: Kunitz Scaffold Stability (Free Energy)

## Abstract

This experiment aims to compute the folding free energy (ΔG_folding) of the BPTI Kunitz scaffold using Free Energy Perturbation (FEP) methods. The thermodynamic cycle compares folded versus unfolded BPTI states in both free and complexed forms. Due to the extensive alchemical molecular dynamics simulations required, this experiment could not be executed on CPU-only hardware. The result is classified as INCONCLUSIVE.

## Introduction / Background

The bovine pancreatic trypsin inhibitor (BPTI) is the prototypical Kunitz-type serine protease inhibitor. Quantifying its folding free energy is fundamental to understanding the thermodynamic stability of the Kunitz scaffold and its role in protease inhibition. Experimental measurements place ΔG_folding at approximately −12.5 ± 2.0 kcal/mol (Moses & Bhattacharyya 1983). Reproducing this value computationally validates the force field parameterization and the FEP protocol used in the pipeline.

**Feature:** F-28 — ΔG_folding for BPTI Kunitz scaffold

**System:** 4PTI (BPTI free) and 2PTC (BPTI in complex with trypsin)

## Hypothesis

FEP calculations with adequate sampling (≥500 ns total) will yield a computed ΔG_folding within ±2.0 kcal/mol of the experimental benchmark of −12.5 kcal/mol, confirming the thermodynamic stability of the Kunitz scaffold as captured by the simulation pipeline.

## Methods

- **Software:** OpenMM with alchemical free energy plugins (openmmtools)
- **Force Field:** AMBER ff14SB + TIP3P explicit solvent
- **Protocol:** FEP with thermodynamic cycle comparing folded vs unfolded BPTI states
- **Lambda Windows:** ~20 evenly spaced windows
- **Simulation Time:** 5 ns per lambda window × 5 independent replicates = 500 ns total
- **Analysis:** MBAR (Multistate Bennett Acceptance Ratio) for free energy estimation
- **Structures:** 4PTI (BPTI free), 2PTC (BPTI–trypsin complex)

## Controls

- **Positive Control:** Known experimental ΔG_folding = −12.5 ± 2.0 kcal/mol (Moses & Bhattacharyya 1983)
- **Convergence Control:** Forward/reverse FEP overlap and MBAR uncertainty estimates
- **Replicate Control:** 5 independent replicates to assess statistical uncertainty

## Results

**INCONCLUSIVE** — FEP calculations require GPU-accelerated OpenMM with alchemical protocols; estimated 10–50+ days on CPU. The simulation was not executed.

## Discussion

The FEP protocol for computing ΔG_folding of BPTI is well-established in the literature but demands extensive GPU-accelerated sampling. The 500 ns total simulation time across 20 lambda windows and 5 replicates is computationally prohibitive on CPU-only hardware. This experiment remains a high-priority target for future GPU execution, as successful reproduction of the experimental benchmark would validate the core thermodynamic framework of the pipeline.

## Conclusions

No conclusions can be drawn. This experiment requires GPU-accelerated alchemical MD simulations that were not feasible on available hardware. The experimental benchmark (ΔG_folding = −12.5 ± 2.0 kcal/mol) remains the target for future validation.

## Figures

No figures generated — full FEP simulation not executed.

## References

1. Moses, E. & Bhattacharyya, P. (1983). Folding thermodynamics of bovine pancreatic trypsin inhibitor. *Journal of Molecular Biology*.
2. Shirts, M. R. & Chodera, J. D. (2008). Statistically optimal analysis of samples from multiple equilibrium states. *Journal of Chemical Physics*, 129, 124105.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-28_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_scaffold (kcal/mol) | [−10.9, −4.9] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-28_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_scaffold (kcal/mol) | [−10.9, −4.9] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
