# EXP-07: P1 Residue Energetic Contribution

## Abstract

This experiment aims to quantify the free energy contribution of the P1 residue (Lys15) to BPTI-trypsin binding through per-residue free energy decomposition from PMF data or alchemical FEP of the K15A mutation. The experiment is INCONCLUSIVE because it depends on EXP-04 (BPTI-trypsin binding free energy), which has not been executed due to computational resource constraints.

## Introduction / Background

The P1 residue at the primary contact site of canonical serine protease inhibitors is the single most important determinant of binding affinity and specificity. In BPTI, the P1 residue is Lys15, which inserts into the trypsin specificity pocket (S1) and forms a salt bridge with Asp189 at the pocket base. Understanding the energetic contribution of this residue is fundamental to the canonical inhibition mechanism.

Feature F-07 quantifies how much of the total BPTI-trypsin binding free energy is attributable to the P1 position, providing a thermodynamic basis for the dominant role of P1 in inhibitor design and engineering.

## Hypothesis

The P1 residue (Lys15) contributes approximately 50% of the total BPTI-trypsin binding free energy, corresponding to 8–10 kcal/mol of the total ΔG_bind ≈ −18 kcal/mol. This is consistent with experimental mutagenesis data showing that K15A mutation reduces binding affinity by several orders of magnitude.

## Methods

- **Per-residue free energy decomposition**: Decompose the PMF from EXP-04 umbrella sampling into per-residue contributions using MM-PBSA or MM-GBSA post-processing of trajectory snapshots along the dissociation coordinate.
- **Alchemical FEP (alternative)**: Perform free energy perturbation of the K15A mutation in both the bound and unbound states to compute ΔΔG_bind(K15A).
- **Force field**: CHARMM36m with TIP3P water.
- **Analysis tools**: GROMACS energy decomposition, MDAnalysis for trajectory processing.

## Controls

- Total ΔG_bind from EXP-04 PMF serves as the reference for decomposition consistency (sum of per-residue contributions should approximate total).
- Comparison with experimental ΔΔG values for P1 mutations from Krowarsch et al. (2003).

## Results

INCONCLUSIVE — this experiment depends on outputs from EXP-04 (BPTI-trypsin binding free energy via SMD/umbrella sampling) which has not been executed due to computational resource constraints. No PMF data or production trajectories are available for free energy decomposition.

## Discussion

Without the PMF and associated trajectory data from EXP-04, per-residue energy decomposition cannot be performed. The expected result — that P1 contributes ~50% of total ΔG_bind — is well-established in the experimental literature through mutagenesis studies. Future execution of EXP-04 on GPU-enabled hardware would enable this analysis. The alchemical FEP approach (K15A mutation) would require additional simulation setup but could provide a more direct measurement of the P1 contribution.

## Conclusions

The P1 energetic contribution analysis remains inconclusive pending execution of the upstream EXP-04 experiment. The experimental benchmark of 8–10 kcal/mol contribution from Lys15 (Krowarsch et al. 2003) provides a clear validation target for when computational resources become available.

## Figures

No figures generated — upstream experiment dependencies unmet.

## References

1. Krowarsch, D., Cierpicki, T., Jelen, F., & Otlewski, J. (2003). Canonical protein inhibitors of serine proteases. *Cellular and Molecular Life Sciences*, 60(11), 2427–2444.
2. Ardèvol, A., Tribello, G. A., & Parrinello, M. (2015). Identification of reaction coordinates from molecular dynamics simulations. *Journal of Chemical Theory and Computation*, 11(3), 1086–1093.

## Author Block

- **Author**: Ryan Kamp
- **Affiliation**: Dept. of Computer Science, University of Cincinnati
- **Email**: kamprj@mail.uc.edu
- **GitHub**: ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-07_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| K15 contact rank | #1 | (observed) | (status) |
| K15 contact fraction | >40% | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-07_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| K15 contact rank | #1 | (observed) | (status) |
| K15 contact fraction | >40% | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
