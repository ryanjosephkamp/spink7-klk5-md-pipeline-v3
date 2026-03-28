# EXP-26: BSA–ΔG Correlation Across Systems

## Abstract

This experiment aims to establish a linear correlation between buried surface area (BSA) and binding free energy (ΔG) across four protein-protein systems: BPTI-trypsin, PSTI-chymotrypsin, barnase-barstar, and SH3-p41. The experiment is INCONCLUSIVE because the ΔG values from EXP-04, EXP-05, EXP-13, and EXP-29 are all unavailable due to computational resource constraints. Static BSA values are available from crystal structure analysis (e.g., EXP-16: 1608.5 Å² for BPTI-trypsin) but the correlation requires computed ΔG values.

## Introduction / Background

The correlation between buried surface area and binding free energy is a foundational concept in structural biology, providing a simple physical model for protein-protein recognition: larger interfaces bury more hydrophobic surface, forming more contacts, and thus bind more tightly. While the correlation is imperfect — specificity arises from the detailed chemical complementarity — a positive correlation across diverse systems validates that the simulation methodology captures the basic thermodynamic driving forces of association.

Feature F-26 tests whether the SMD/umbrella sampling pipeline produces ΔG values that correlate with BSA across four structurally distinct protein-protein systems, benchmarked against the empirical observation of R² > 0.7 (Janin 1997, Horton & Lewis 1992).

## Hypothesis

A linear regression of BSA vs. computed ΔG across the four systems (BPTI-trypsin, PSTI-chymotrypsin, barnase-barstar, SH3-p41) will yield R² > 0.7, consistent with the empirical correlation reported by Janin (1997). The slope is expected to be approximately −20 to −30 cal/mol/Å², reflecting the hydrophobic contribution to binding.

## Methods

- **BSA values**: Computed from crystal structures using FreeSASA or MDAnalysis SASA (probe radius 1.4 Å). BSA = SASA_A + SASA_B − SASA_AB.
- **ΔG values**: From SMD/umbrella sampling PMF integration in EXP-04 (BPTI-trypsin), EXP-05 (PSTI-chymotrypsin), EXP-13 (barnase-barstar), and EXP-29 (SH3-p41).
- **Linear regression**: Ordinary least squares fit of BSA vs. ΔG, with R², slope, intercept, and p-value.
- **Systems**: BPTI-trypsin (3OTJ), PSTI-chymotrypsin, barnase-barstar (1BRS), SH3-p41 (1SSH or similar).

## Controls

- Comparison with experimental ΔG values from literature as an independent benchmark.
- BSA computed from both crystal structures and MD-averaged conformations (when available).
- Leave-one-out cross-validation of the linear fit to assess robustness.
- Comparison of computed slope with literature values (~25 cal/mol/Å²).

## Results

INCONCLUSIVE — this experiment depends on ΔG values from EXP-04 (BPTI-trypsin), EXP-05 (PSTI-chymotrypsin), EXP-13 (barnase-barstar), and EXP-29 (SH3-p41), all of which have not been executed due to computational resource constraints. Static BSA values are available (e.g., 1608.5 Å² for BPTI-trypsin from EXP-16) but cannot be correlated without computed ΔG counterparts.

## Discussion

The BSA-ΔG correlation is a system-level validation that requires successful execution of the entire SMD/umbrella sampling pipeline across all four protein-protein systems. This is the most resource-intensive experiment in the validation suite, as it aggregates results from four independent free energy calculations. The availability of BSA values from crystal structure analysis demonstrates that the structural component of the pipeline is functional; the bottleneck is exclusively the free energy computation requiring GPU-accelerated molecular dynamics.

The expected correlation, while well-established empirically, has significant scatter. With only four data points, statistical power is limited, and outliers could substantially affect R². Nevertheless, a positive correlation would provide strong evidence that the pipeline captures the fundamental thermodynamic trends across diverse systems.

## Conclusions

The BSA-ΔG correlation analysis remains inconclusive pending execution of the upstream free energy calculations (EXP-04, 05, 13, 29). The experimental benchmark of R² > 0.7 (Janin 1997, Horton & Lewis 1992) provides the validation target. BSA values from crystal structures are available and ready for correlation when ΔG values are computed.

## Figures

No figures generated — upstream experiment dependencies unmet.

## References

1. Janin, J. (1997). The kinetics of protein-protein recognition. *Proteins: Structure, Function, and Genetics*, 28(2), 153–161.
2. Horton, N., & Lewis, M. (1992). Calculation of the free energy of association for protein complexes. *Protein Science*, 1(1), 169–181.
3. Chothia, C., & Janin, J. (1975). Principles of protein-protein recognition. *Nature*, 256(5520), 705–708.

## Author Block

- **Author**: Ryan Kamp
- **Affiliation**: Dept. of Computer Science, University of Cincinnati
- **Email**: kamprj@mail.uc.edu
- **GitHub**: ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-26_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| R² | >0.5 | (observed) | (status) |
| Slope (kcal/(mol·Å²)) | [−0.005, −0.020] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-26_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| R² | >0.5 | (observed) | (status) |
| Slope (kcal/(mol·Å²)) | [−0.005, −0.020] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
