# EXP-32: SPINK1 N34S Negative Control (Disease Mutation)

## Abstract

This experiment uses Free Energy Perturbation (FEP) to compute the ΔΔG of the N34S mutation in SPINK1, the most common variant associated with hereditary pancreatitis. Using the porcine SPINK1 homolog (PSTI, PDB: 1TGS) as a structural proxy, the pipeline should predict a modest destabilization of 0.5–2.0 kcal/mol, serving as a negative control. Due to the alchemical MD simulations required (~200 ns), this experiment could not be executed on CPU-only hardware. The result is classified as INCONCLUSIVE.

## Introduction / Background

SPINK1 (serine peptidase inhibitor, Kazal type 1) is a pancreatic secretory trypsin inhibitor that prevents premature trypsin activation within the pancreas. The N34S variant is the most frequently observed SPINK1 mutation in patients with hereditary and idiopathic chronic pancreatitis (Witt et al. 2000). Despite its strong genetic association with disease, N34S causes only a modest effect on inhibitor–protease binding affinity (ΔΔG ≈ 0.5–2.0 kcal/mol), suggesting that its pathogenic mechanism may involve folding, secretion, or other non-binding effects. This makes N34S a valuable negative control: the pipeline should predict a small but measurable ΔΔG, clearly distinguishable from binding-hotspot mutations like BPTI K15A.

**Feature:** F-32 — ΔΔG for SPINK1 N34S mutation (pancreatitis-associated)

**System:** 1TGS (PSTI = porcine SPINK1 homolog, as structural proxy)

## Hypothesis

FEP calculations for the N34S mutation will yield a ΔΔG of 0.5–2.0 kcal/mol (modest destabilization), consistent with the experimentally observed mild effect on binding, and clearly smaller in magnitude than hotspot mutations such as K15A.

## Methods

- **Software:** OpenMM with alchemical free energy plugins (openmmtools)
- **Force Field:** AMBER ff14SB + TIP3P explicit solvent
- **Protocol:** FEP for N→S mutation at position 34 (or equivalent in PSTI numbering) using dual-topology approach (bound complex + free inhibitor legs)
- **Lambda Windows:** ~20 evenly spaced windows
- **Simulation Time:** ~200 ns total (bound + free legs × replicates)
- **Analysis:** MBAR for free energy estimation; thermodynamic cycle for ΔΔG_bind
- **Structure:** 1TGS (PSTI–trypsin complex)

## Controls

- **Literature Benchmark:** N34S ΔΔG = 0.5–2.0 kcal/mol (modest destabilization) — Witt et al. 2000
- **Negative Control Logic:** ΔΔG should be significantly smaller than hotspot residues (e.g., K15A > 10 kcal/mol)
- **Convergence Control:** Forward/reverse overlap and MBAR uncertainty estimates

## Results

**INCONCLUSIVE** — FEP calculations require GPU-accelerated OpenMM with alchemical protocols; estimated 10–50+ days on CPU. The simulation was not executed.

## Discussion

The SPINK1 N34S mutation serves as an important negative control in the experimental validation suite. Its modest effect on binding (0.5–2.0 kcal/mol) tests whether the pipeline can distinguish between mild and severe mutations — a critical capability for any predictive tool. Using 1TGS (porcine PSTI) as a structural proxy introduces some uncertainty due to sequence differences with human SPINK1, but the Kazal domain fold is well-conserved. The ~200 ns total simulation time, while modest compared to other FEP experiments in this suite, remains infeasible on CPU-only hardware. This experiment is a priority for future GPU execution due to its clinical relevance and role as a calibration point.

## Conclusions

No conclusions can be drawn. This experiment requires GPU-accelerated alchemical MD simulations that were not feasible on available hardware. The experimental benchmark (ΔΔG = 0.5–2.0 kcal/mol) remains the target for future validation.

## Figures

No figures generated — full FEP simulation not executed.

## References

1. Witt, H., Luck, W., Hennies, H. C., Claßen, M., Kage, A., Laß, U., Landt, O., & Becker, M. (2000). Mutations in the gene encoding the serine protease inhibitor, Kazal type 1 are associated with chronic pancreatitis. *Nature Genetics*, 25(2), 213–216.
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
**Notebook:** `EXP-32_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| |ΔΔG| (kcal/mol) | < 1.0 | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-32_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| |ΔΔG| (kcal/mol) | < 1.0 | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
