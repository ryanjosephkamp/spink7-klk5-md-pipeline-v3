# EXP-24: BPTI Hydrogen-Deuterium Exchange Protection Factors

## Abstract

This experiment aims to predict backbone NH protection factors from an extended production MD trajectory of BPTI, validating against the well-characterized experimental H/D exchange data. The experiment is INCONCLUSIVE because it requires a 100 ns production trajectory from EXP-23 that is not feasible on CPU-only hardware. The crystal structure analysis component of EXP-23 was completed successfully (PASS), confirming correct structural representation of BPTI in the simulation pipeline.

## Introduction / Background

Hydrogen-deuterium (H/D) exchange experiments measure the rate at which backbone amide protons exchange with solvent deuterium, providing residue-level information about structural dynamics and solvent accessibility. Protection factors (Pf) quantify the ratio of the intrinsic (random coil) exchange rate to the observed rate, with high protection indicating persistent hydrogen bonding and burial from solvent. BPTI is a classic system for H/D exchange studies, with protection factors spanning more than 8 orders of magnitude.

Feature F-24 computes backbone NH protection factors from MD trajectory analysis and compares against the experimental data of Wagner & Wüthrich (1982).

## Hypothesis

Core residues of BPTI (β-sheet and α-helix) will show log(Pf) > 5, reflecting persistent hydrogen bonding and solvent exclusion, while surface loops will show log(Pf) < 2, consistent with rapid exchange. The overall pattern of protection factors will correlate with the Wagner & Wüthrich (1982) experimental data.

## Methods

- **Production MD**: 100 ns trajectory of BPTI in explicit solvent (TIP3P) from EXP-23 equilibrated system.
- **H-bond occupancy**: Backbone NH hydrogen bond occupancy computed using geometric criteria (N-H···O distance < 3.5 Å, angle > 135°) over the production trajectory.
- **SASA calculation**: Per-residue backbone NH solvent-accessible surface area computed for each frame.
- **Protection factor estimation**: Pf estimated from the Vendruscolo model: ln(Pf) = β_c × N_c + β_h × N_h, where N_c is the number of heavy-atom contacts and N_h is the H-bond occupancy.
- **Force field**: CHARMM36m, TIP3P water, 150 mM NaCl.

## Controls

- Crystal structure B-factors as a proxy for mobility (available from EXP-23 PASS result).
- Intrinsic exchange rates from the Bai et al. (1993) model as the unprotected reference.
- Correlation coefficient between predicted and experimental log(Pf) values as the primary metric.

## Results

INCONCLUSIVE — this experiment depends on outputs from EXP-23 (100 ns production MD trajectory of BPTI) which has not been executed due to computational resource constraints. The crystal structure analysis component of EXP-23 was completed successfully (PASS), confirming that the pipeline correctly represents BPTI structure, but the extended dynamics required for H/D exchange prediction are unavailable.

## Discussion

Protection factor prediction from MD simulation requires extensive sampling to accurately capture the fluctuations that govern amide exchange. A 100 ns trajectory represents the minimum timescale needed for reliable H-bond occupancy statistics across all residues. While the successful crystal structure validation in EXP-23 confirms the structural fidelity of the pipeline, the dynamic information essential for H/D exchange prediction cannot be obtained without the production trajectory. Short equilibration trajectories are insufficient for this analysis due to inadequate sampling of rare opening events in the protein core.

## Conclusions

The H/D exchange protection factor prediction remains inconclusive pending execution of the 100 ns production trajectory in EXP-23. The experimental benchmark from Wagner & Wüthrich (1982) — core log(Pf) > 5, surface log(Pf) < 2 — provides clear validation targets. The successful structural validation from EXP-23 crystal structure analysis supports the pipeline's readiness for this calculation when GPU resources become available.

## Figures

No figures generated — upstream experiment dependencies unmet.

## References

1. Wagner, G., & Wüthrich, K. (1982). Amide proton exchange and surface conformation of the basic pancreatic trypsin inhibitor in solution. *Journal of Molecular Biology*, 160(2), 343–361.
2. Vendruscolo, M., Paci, E., Dobson, C. M., & Karplus, M. (2003). Rare fluctuations of native proteins sampled by equilibrium hydrogen exchange. *Journal of the American Chemical Society*, 125(51), 15686–15687.
3. Bai, Y., Milne, J. S., Mayne, L., & Englander, S. W. (1993). Primary structure effects on peptide group hydrogen exchange. *Proteins: Structure, Function, and Genetics*, 17(1), 75–86.

## Author Block

- **Author**: Ryan Kamp
- **Affiliation**: Dept. of Computer Science, University of Cincinnati
- **Email**: kamprj@mail.uc.edu
- **GitHub**: ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-24_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| True positives | ≥8/11 | (observed) | (status) |
| False positives | ≤3 | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-24_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| True positives | ≥8/11 | (observed) | (status) |
| False positives | ≤3 | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
