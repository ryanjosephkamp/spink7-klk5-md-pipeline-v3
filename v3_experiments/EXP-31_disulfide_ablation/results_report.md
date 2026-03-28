# EXP-31: Disulfide Bond Ablation (C14S/C38S)

## Abstract

This experiment uses Free Energy Perturbation (FEP) to compute the ΔΔG of disulfide bond ablation mutations (C14S and C38S) in BPTI within the BPTI–trypsin complex. These mutations disrupt the Cys14–Cys38 disulfide bond, which is expected to cause 3–5 kcal/mol destabilization. Due to the extensive alchemical MD simulations required (~400 ns per mutation pair) with specialized disulfide-breaking alchemical pathways, this experiment could not be executed on CPU-only hardware. The result is classified as INCONCLUSIVE.

## Introduction / Background

BPTI contains three disulfide bonds (Cys5–Cys55, Cys14–Cys38, Cys30–Cys51) that are critical to its structural integrity and inhibitory function. The Cys14–Cys38 disulfide bridges the binding loop region to the protein core. Ablation of this bond via C14S or C38S mutations is expected to destabilize the protein by 3–5 kcal/mol (Goldenberg 1988, Otzen et al. 1999). Computing this destabilization energy via FEP requires a specialized alchemical pathway that accounts for the breaking of the covalent disulfide bond during the transformation.

**Feature:** F-31 — ΔΔG for disulfide ablation mutations in BPTI

**System:** 2PTC (BPTI–trypsin complex)

## Hypothesis

FEP calculations for C14S and C38S mutations will yield ΔΔG values in the range of 3–5 kcal/mol destabilization, consistent with experimental measurements of disulfide bond ablation effects on BPTI stability.

## Methods

- **Software:** OpenMM with alchemical free energy plugins (openmmtools)
- **Force Field:** AMBER ff14SB + TIP3P explicit solvent
- **Protocol:** FEP for C14S and C38S mutations with dual-topology approach (complex + free legs); specialized handling for disulfide bond breaking in the alchemical pathway
- **Lambda Windows:** ~20 evenly spaced windows per mutation
- **Simulation Time:** ~400 ns per mutation pair (bound + free legs × replicates)
- **Analysis:** MBAR for free energy estimation; thermodynamic cycle for ΔΔG_bind
- **Structure:** 2PTC (BPTI–trypsin complex)

## Controls

- **Literature Benchmark:** C14S/C38S → ΔΔG = 3–5 kcal/mol destabilization (Goldenberg 1988, Otzen et al. 1999)
- **Convergence Control:** Forward/reverse overlap and MBAR uncertainty estimates
- **Structural Validation:** Monitoring of disulfide geometry during alchemical transformation

## Results

**INCONCLUSIVE** — FEP calculations require GPU-accelerated OpenMM with alchemical protocols; estimated 10–50+ days on CPU. The simulation was not executed.

## Discussion

Disulfide bond ablation via FEP presents unique technical challenges beyond standard point mutations, as the alchemical pathway must handle the breaking of a covalent bond. This requires careful staging of the lambda schedule to avoid numerical instabilities. The ~400 ns simulation time per mutation pair, combined with the specialized protocol requirements, makes this experiment infeasible on CPU-only hardware. Successfully computing the expected 3–5 kcal/mol destabilization would validate the pipeline's ability to handle non-standard alchemical transformations involving covalent modifications.

## Conclusions

No conclusions can be drawn. This experiment requires GPU-accelerated alchemical MD simulations with specialized disulfide-breaking protocols that were not feasible on available hardware. The experimental benchmark (ΔΔG = 3–5 kcal/mol) remains the target for future validation.

## Figures

No figures generated — full FEP simulation not executed.

## References

1. Goldenberg, D. P. (1988). Kinetic analysis of the folding and unfolding of a mutant form of bovine pancreatic trypsin inhibitor lacking the cysteine-14 and -38 thiols. *Biochemistry*, 27(8), 2481–2489.
2. Otzen, D. E., Itzhaki, L. S., ElMasry, N. F., Jackson, S. E., & Fersht, A. R. (1999). Structure of the transition state for the folding/unfolding of the barley chymotrypsin inhibitor 2 and its implications for mechanisms of protein folding. *Proceedings of the National Academy of Sciences*, 91(22), 10422–10425.
3. Shirts, M. R. & Chodera, J. D. (2008). Statistically optimal analysis of samples from multiple equilibrium states. *Journal of Chemical Physics*, 129, 124105.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-31_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔΔG (kcal/mol) | [+4, +10] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-31_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔΔG (kcal/mol) | [+4, +10] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
