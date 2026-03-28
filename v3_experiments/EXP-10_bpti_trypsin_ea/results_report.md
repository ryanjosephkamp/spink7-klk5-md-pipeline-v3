# EXP-10: BPTI-Trypsin Activation Energy

## Abstract

This experiment aims to determine the activation energy (ΔG‡) for BPTI-trypsin dissociation from the PMF barrier height obtained via umbrella sampling, with friction corrections via Grote-Hynes theory. The experiment is INCONCLUSIVE because it depends on the PMF from EXP-04, which has not been executed due to computational resource constraints.

## Introduction / Background

The extraordinarily tight binding of BPTI to trypsin (K_i ≈ 6 × 10⁻¹⁴ M) implies a very large activation barrier for dissociation. Estimating ΔG‡ from the PMF barrier height provides a direct test of whether the simulated free energy landscape is consistent with the experimentally observed kinetic stability of the complex. The activation energy is a key parameter connecting thermodynamic profiles to kinetic observables via transition state theory (TST) and its corrections.

Feature F-10 extracts ΔG‡ from the EXP-04 PMF and applies Grote-Hynes friction corrections to account for recrossing effects in the condensed-phase dissociation process.

## Hypothesis

The activation energy for BPTI-trypsin dissociation is ΔG‡ ≈ 20–25 kcal/mol, estimated from the experimental k_off derived from K_i × k_on (K_i ≈ 6 × 10⁻¹⁴ M, k_on ≈ 10⁶ M⁻¹s⁻¹, giving k_off ≈ 6 × 10⁻⁸ s⁻¹) and transition state theory. The PMF barrier height from umbrella sampling should recover this value within ±3 kcal/mol.

## Methods

- **PMF barrier height**: Extract the maximum barrier along the dissociation coordinate from the EXP-04 umbrella sampling PMF (WHAM-reconstructed).
- **TST estimate**: ΔG‡_TST = −kT × ln(k_off × h / kT), where k_off is derived from experimental K_i and k_on.
- **Grote-Hynes correction**: Apply frequency-dependent friction correction to account for barrier recrossing: κ_GH = ω_b / ω_eff, where ω_b is the barrier frequency and ω_eff is the effective frequency including solvent friction.
- **Barrier shape analysis**: Fit the PMF barrier region to a parabola to extract the barrier frequency ω_b.

## Controls

- TST-derived ΔG‡ from experimental k_off as independent benchmark.
- Sensitivity to umbrella sampling window spacing and simulation length.
- Comparison of uncorrected (TST) and corrected (Grote-Hynes) barrier heights.

## Results

INCONCLUSIVE — this experiment depends on outputs from EXP-04 (BPTI-trypsin PMF barrier from umbrella sampling) which has not been executed due to computational resource constraints. No PMF data is available for barrier height extraction or Grote-Hynes correction.

## Discussion

The activation energy for BPTI-trypsin dissociation is among the largest known for protein-protein interactions, reflecting the exceptional kinetic stability of this complex. Accurately reproducing ΔG‡ from simulation requires a well-converged PMF, which demands extensive umbrella sampling with sufficient overlap between windows. Without the EXP-04 PMF, neither the barrier height nor the Grote-Hynes friction correction can be computed. Future GPU-enabled execution of EXP-04 would provide the necessary data for this analysis.

## Conclusions

The BPTI-trypsin activation energy determination remains inconclusive pending generation of the PMF from upstream EXP-04. The expected ΔG‡ ≈ 20–25 kcal/mol, derived from experimental kinetic data, provides a quantitative validation target.

## Figures

No figures generated — upstream experiment dependencies unmet.

## References

1. Grote, R. F., & Hynes, J. T. (1980). The stable states picture of chemical reactions. II. Rate constants for condensed and gas phase reaction models. *Journal of Chemical Physics*, 73(6), 2715–2732.
2. Vincent, J. P., & Bhattacharyya, A. (1982). Kinetics of the interaction of bovine pancreatic trypsin inhibitor with trypsin. *Biochemistry*, 21(12), 2844–2849.
3. Szabo, A., Schulten, K., & Schulten, Z. (1980). First passage time approach to diffusion controlled reactions. *Journal of Chemical Physics*, 72(8), 4350–4357.

## Author Block

- **Author**: Ryan Kamp
- **Affiliation**: Dept. of Computer Science, University of Cincinnati
- **Email**: kamprj@mail.uc.edu
- **GitHub**: ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-10_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Ea (kcal/mol) | [4.6, 16.4] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-10_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Ea (kcal/mol) | [4.6, 16.4] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
