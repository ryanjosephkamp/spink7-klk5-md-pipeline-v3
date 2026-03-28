# EXP-09: BPTI-Trypsin Association Kinetics

## Abstract

This experiment aims to compute the association rate constant (k_on) for BPTI-trypsin binding by applying Kramers theory to a 1D potential of mean force (PMF) obtained from steered molecular dynamics and umbrella sampling. The experiment is INCONCLUSIVE because it depends on the PMF from EXP-04, which has not been executed due to computational resource constraints.

## Introduction / Background

The association kinetics of protease-inhibitor complexes are fundamental to understanding the mechanism of protease inhibition in vivo. BPTI is one of the fastest-associating protein-protein systems known, with k_on on the order of 10⁶ M⁻¹s⁻¹. Reproducing this rate from molecular simulation provides a stringent test of both the underlying free energy surface and the diffusional dynamics captured by the simulation methodology.

Feature F-09 computes k_on for BPTI-trypsin association from the PMF using Kramers theory, connecting the thermodynamic landscape to kinetic observables.

## Hypothesis

The computed k_on from Kramers theory applied to the EXP-04 PMF will be on the order of 10⁶ M⁻¹s⁻¹, consistent with the experimental measurement by Vincent & Bhattacharyya (1982). The rate is expected to be near the diffusion-controlled limit, modulated by the electrostatic steering between the positively charged Lys15 and the negatively charged trypsin specificity pocket.

## Methods

- **PMF source**: 1D PMF along the BPTI-trypsin center-of-mass separation coordinate from EXP-04 umbrella sampling.
- **Kramers theory**: k_on = (D_eff × exp(−ΔG‡/kT)) / ∫ exp(W(r)/kT) dr, where W(r) is the PMF, D_eff is the effective diffusion coefficient, and the integral runs over the reaction coordinate.
- **Diffusion coefficient**: Estimated from the mean-square displacement of the BPTI center-of-mass in bulk solvent windows, or from the friction coefficient via the fluctuation-dissipation theorem.
- **Corrections**: Rotational averaging, hydrodynamic corrections for finite box size.

## Controls

- Comparison with experimental k_on ≈ 10⁶ M⁻¹s⁻¹ (Vincent & Bhattacharyya 1982).
- Sensitivity analysis: variation of D_eff by ±50% to assess impact on computed k_on.
- Smoluchowski diffusion-limited rate as an upper bound.

## Results

INCONCLUSIVE — this experiment depends on outputs from EXP-04 (BPTI-trypsin PMF from umbrella sampling) which has not been executed due to computational resource constraints. No PMF data is available for Kramers theory rate calculation.

## Discussion

Kramers theory provides a rigorous framework for extracting kinetic rates from equilibrium free energy profiles, but requires a well-converged PMF as input. Without the EXP-04 PMF, the rate calculation cannot be performed. The expected result — k_on ≈ 10⁶ M⁻¹s⁻¹ — is consistent with the near-diffusion-limited association expected for electrostatically steered protein-protein interactions. Future execution of EXP-04 would enable this kinetic analysis and provide a direct connection between the simulated free energy landscape and experimentally measurable kinetics.

## Conclusions

The BPTI-trypsin association kinetics calculation remains inconclusive pending generation of the PMF from upstream EXP-04. The experimental benchmark of k_on ≈ 10⁶ M⁻¹s⁻¹ (Vincent & Bhattacharyya 1982) provides a quantitative validation target.

## Figures

No figures generated — upstream experiment dependencies unmet.

## References

1. Vincent, J. P., & Bhattacharyya, A. (1982). Kinetics of the interaction of bovine pancreatic trypsin inhibitor with trypsin. *Biochemistry*, 21(12), 2844–2849.
2. Kramers, H. A. (1940). Brownian motion in a field of force and the diffusion model of chemical reactions. *Physica*, 7(4), 284–304.
3. Northrup, S. H., & Erickson, H. P. (1992). Kinetics of protein-protein association explained by Brownian dynamics computer simulation. *Proceedings of the National Academy of Sciences*, 89(8), 3338–3342.

## Author Block

- **Author**: Ryan Kamp
- **Affiliation**: Dept. of Computer Science, University of Cincinnati
- **Email**: kamprj@mail.uc.edu
- **GitHub**: ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-09_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| k_on agreement | within 2 orders of ~10⁶ M⁻¹s⁻¹ | (observed) | (status) |
| k_off agreement | within 2 orders of ~10⁻⁸ s⁻¹ | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-09_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| k_on agreement | within 2 orders of ~10⁶ M⁻¹s⁻¹ | (observed) | (status) |
| k_off agreement | within 2 orders of ~10⁻⁸ s⁻¹ | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
