# EXP-08: Interfacial H-Bond Energetics

## Abstract

This experiment aims to quantify the per-hydrogen-bond energy contribution at the BPTI-trypsin interface through H-bond occupancy analysis from production MD and energy decomposition from free energy perturbation. The experiment is INCONCLUSIVE because it depends on EXP-04 production trajectories, which have not been generated due to computational resource constraints.

## Introduction / Background

Hydrogen bonds at protein-protein interfaces are key determinants of binding specificity and affinity. The BPTI-trypsin interface features a network of direct and water-mediated hydrogen bonds that stabilize the inhibitor-enzyme complex. Quantifying the energetic contribution of individual H-bonds provides insight into which interfacial contacts are most critical for recognition and could guide inhibitor engineering.

Feature F-08 measures the per-H-bond energy contribution at the BPTI-trypsin interface, benchmarked against the established range of 1–3 kcal/mol per interfacial hydrogen bond.

## Hypothesis

Each interfacial hydrogen bond at the BPTI-trypsin interface contributes 1–3 kcal/mol to the binding free energy, consistent with Fersht's (1987) systematic analysis of hydrogen bond contributions to protein stability and binding. High-occupancy H-bonds (>80%) are expected to contribute toward the upper end of this range.

## Methods

- **H-bond occupancy analysis**: Identify and track all interfacial H-bonds from production MD trajectories (EXP-04) using geometric criteria (donor-acceptor distance < 3.5 Å, angle > 135°).
- **Energy decomposition**: Per-H-bond energy contribution estimated via MM-GBSA decomposition or systematic in silico mutagenesis (donor/acceptor → Ala) with FEP.
- **Tools**: MDAnalysis HydrogenBondAnalysis, GROMACS hbond, VMD for visualization.
- **Force field**: CHARMM36m with TIP3P water.

## Controls

- Comparison of H-bond occupancies with crystal structure contacts (PDB: 3OTJ).
- Total interfacial H-bond energy should be consistent with the overall ΔG_bind from EXP-04.
- Benchmark against Fersht (1987) empirical range of 1–3 kcal/mol per H-bond.

## Results

INCONCLUSIVE — this experiment depends on outputs from EXP-04 (BPTI-trypsin production MD trajectory and free energy data) which has not been executed due to computational resource constraints. No trajectory data is available for H-bond occupancy analysis or energy decomposition.

## Discussion

Interfacial H-bond analysis requires equilibrated production trajectories of sufficient length to compute reliable occupancy statistics. Without EXP-04 trajectories, neither the dynamic H-bond network characterization nor the per-bond energy decomposition can be performed. Crystal structure analysis can identify potential H-bonds but cannot provide the dynamic occupancy or energetic data needed for this experiment. Future GPU-enabled execution of EXP-04 would provide the necessary trajectory data.

## Conclusions

The interfacial H-bond energetics analysis remains inconclusive pending execution of upstream EXP-04. The Fersht (1987) benchmark of 1–3 kcal/mol per H-bond provides a well-established validation target for future analysis.

## Figures

No figures generated — upstream experiment dependencies unmet.

## References

1. Fersht, A. R. (1987). The hydrogen bond in molecular recognition. *Trends in Biochemical Sciences*, 12, 301–304.
2. Fersht, A. R., Shi, J. P., Knill-Jones, J., Lowe, D. M., Wilkinson, A. J., Blow, D. M., ... & Winter, G. (1985). Hydrogen bonding and biological specificity analysed by protein engineering. *Nature*, 314(6008), 235–238.

## Author Block

- **Author**: Ryan Kamp
- **Affiliation**: Dept. of Computer Science, University of Cincinnati
- **Email**: kamprj@mail.uc.edu
- **GitHub**: ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-08_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Per-H-bond energy (kcal/mol) | [0.7, 2.3] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-08_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Per-H-bond energy (kcal/mol) | [0.7, 2.3] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
