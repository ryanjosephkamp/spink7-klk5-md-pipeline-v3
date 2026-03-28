# EXP-33: Barnase–Barstar Double Mutant Cycles

## Abstract

This experiment applies Free Energy Perturbation (FEP) to compute interaction energies (ΔΔG_int) from double mutant cycle analysis at the barnase–barstar interface. With ~45 experimentally characterized double mutant cycles requiring ~90+ individual mutations and an estimated 18 μs of total simulation time, this represents the most computationally demanding experiment in the validation suite. Due to the massive alchemical MD requirements, this experiment could not be executed on CPU-only hardware. The result is classified as INCONCLUSIVE.

## Introduction / Background

Double mutant cycle analysis is the gold-standard experimental method for quantifying the thermodynamic coupling between pairs of residues at protein–protein interfaces. The barnase–barstar system (PDB: 1BRS) has been exhaustively characterized by Schreiber & Fersht (1995), who measured ΔΔG_int values for ~45 residue pairs. These interaction energies reveal cooperative networks at the interface that cannot be captured by single-mutation analysis alone. Computationally reproducing these coupling energies via FEP would provide the most stringent possible validation of the pipeline's accuracy for predicting mutational effects at protein–protein interfaces.

**Feature:** F-33 — ΔΔG_int (interaction energies) from double mutant cycle analysis

**System:** 1BRS (barnase–barstar complex)

## Hypothesis

FEP-based double mutant cycle analysis will yield ΔΔG_int values that correlate (R² > 0.5) with the experimental measurements of Schreiber & Fersht (1995), correctly identifying strongly coupled residue pairs at the barnase–barstar interface.

## Methods

- **Software:** OpenMM with alchemical free energy plugins (openmmtools)
- **Force Field:** AMBER ff14SB + TIP3P explicit solvent
- **Protocol:** FEP for each single and double mutant using dual-topology approach (bound complex + free protein legs); double mutant cycle: ΔΔG_int = ΔΔG_AB − ΔΔG_A − ΔΔG_B
- **Lambda Windows:** ~20 evenly spaced windows per mutation
- **Simulation Time:** ~200 ns per mutation × ~90+ mutations = ~18,000 ns (18 μs) total
- **Analysis:** MBAR for free energy estimation; double mutant cycle thermodynamic analysis for ΔΔG_int
- **Structure:** 1BRS (barnase–barstar complex)
- **Scope:** ~45 double mutant cycles covering key interface residue pairs

## Controls

- **Literature Benchmark:** ~45 experimentally determined ΔΔG_int values (Schreiber & Fersht 1995)
- **Internal Consistency:** Thermodynamic cycle closure for each double mutant cycle
- **Convergence Control:** Forward/reverse overlap and MBAR uncertainty estimates per mutation
- **Single Mutation Validation:** Individual ΔΔG values compared to experimental single-mutant data

## Results

**INCONCLUSIVE** — FEP calculations require GPU-accelerated OpenMM with alchemical protocols; estimated 10–50+ days on CPU. The simulation was not executed.

## Discussion

The barnase–barstar double mutant cycle experiment is the most computationally intensive in the entire validation suite, requiring approximately 18 μs of aggregate simulation time across ~90+ individual alchemical transformations. This is well beyond the capability of CPU-only hardware. However, it is also the most definitive test of the pipeline's FEP accuracy: successfully reproducing the Schreiber & Fersht coupling energies would demonstrate that the pipeline captures not just individual mutational effects but also the cooperative thermodynamic networks that govern protein–protein recognition. This experiment is the highest-priority target for future GPU execution, as it would serve as the most rigorous benchmark for the pipeline's predictive capability at protein–protein interfaces.

## Conclusions

No conclusions can be drawn. This experiment requires GPU-accelerated alchemical MD simulations that were not feasible on available hardware. The ~45 experimental ΔΔG_int benchmarks from Schreiber & Fersht (1995) remain targets for future validation.

## Figures

No figures generated — full FEP simulation not executed.

## References

1. Schreiber, G. & Fersht, A. R. (1995). Energetics of protein–protein interactions: analysis of the barnase–barstar interface by single mutations and double mutant cycles. *Journal of Molecular Biology*, 248(2), 478–486.
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
**Notebook:** `EXP-33_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Correlation R | >0.7 | (observed) | (status) |
| RMSE (kcal/mol) | <2.0 | (observed) | (status) |
| Top 5 coupled pairs | correct | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-33_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Correlation R | >0.7 | (observed) | (status) |
| RMSE (kcal/mol) | <2.0 | (observed) | (status) |
| Top 5 coupled pairs | correct | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
