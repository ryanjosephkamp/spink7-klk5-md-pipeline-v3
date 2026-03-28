# EXP-25: BPTI Conformational Variability

## Abstract

This experiment aims to compute per-residue RMSF and B-factors from a production MD trajectory of BPTI and compare against crystallographic B-factors. The experiment is INCONCLUSIVE because it requires a 100 ns production trajectory that is not feasible on CPU-only hardware.

## Introduction / Background

Root-mean-square fluctuation (RMSF) analysis from molecular dynamics provides residue-level information about conformational dynamics that can be directly compared with crystallographic B-factors. BPTI is a well-characterized system with high-resolution crystal structures available, making it an ideal benchmark for validating the ability of MD simulations to reproduce experimentally observed flexibility patterns.

Feature F-25 computes per-residue RMSF from a 100 ns production trajectory and correlates the results with crystallographic B-factors, testing whether the simulation force field and protocol correctly capture the dynamics of a small, stable protein.

## Hypothesis

Core residues (β-sheet, α-helix) will show RMSF < 0.5 Å, while surface loops (particularly the binding loop residues 11–18 and the C-terminal region) will show RMSF of 1–3 Å. The per-residue RMSF profile will correlate with crystallographic B-factors (B = 8π²/3 × ⟨Δr²⟩) with R² > 0.6, consistent with Liepinsh et al. (1992).

## Methods

- **Production MD**: 100 ns trajectory of BPTI (PDB: 4PTI or from EXP-23 prepared system) in explicit solvent.
- **RMSF calculation**: Per-residue Cα RMSF computed after alignment to the average structure, using the final 80 ns (discarding 20 ns equilibration).
- **B-factor conversion**: Convert RMSF to computed B-factors via B = 8π²/3 × RMSF².
- **Correlation**: Pearson and Spearman correlation between computed and experimental B-factors.
- **Force field**: CHARMM36m, TIP3P water, 150 mM NaCl.
- **Tools**: MDAnalysis for RMSF computation, matplotlib for visualization.

## Controls

- Crystallographic B-factors from PDB: 4PTI (1.0 Å resolution) as the reference.
- Block averaging to assess convergence of RMSF values.
- Comparison of RMSF from first and second halves of trajectory to verify equilibration.

## Results

INCONCLUSIVE — this experiment depends on a production MD trajectory of BPTI (100 ns, from EXP-23 or direct 4PTI setup) which has not been generated due to computational resource constraints. No trajectory data is available for RMSF analysis.

## Discussion

RMSF analysis is one of the most fundamental validation tests for MD simulations, requiring only a single equilibrated trajectory of sufficient length. A 100 ns trajectory is needed to adequately sample loop motions and ensure converged RMSF values across all residues. While shorter trajectories might provide qualitatively correct results for core residues, loop regions require extended sampling. The CHARMM36m force field has been shown to reproduce BPTI dynamics accurately in the literature, so the expected correlation with crystallographic B-factors is high. Future GPU-enabled execution would readily enable this analysis.

## Conclusions

The BPTI conformational variability analysis remains inconclusive pending generation of a 100 ns production trajectory. The experimental benchmark — core RMSF < 0.5 Å, loop RMSF 1–3 Å (Liepinsh et al. 1992) — provides clear validation targets for future analysis.

## Figures

No figures generated — upstream experiment dependencies unmet.

## References

1. Liepinsh, E., Otting, G., & Wüthrich, K. (1992). NMR spectroscopy of hydroxyl protons in aqueous solutions of peptides and proteins. *Journal of Biomolecular NMR*, 2(5), 447–465.
2. van Gunsteren, W. F., & Mark, A. E. (1998). Validation of molecular dynamics simulation. *Journal of Chemical Physics*, 108(15), 6109–6116.

## Author Block

- **Author**: Ryan Kamp
- **Affiliation**: Dept. of Computer Science, University of Cincinnati
- **Email**: kamprj@mail.uc.edu
- **GitHub**: ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-25_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Cα RMSD fluctuation (Å) | [0.25, 0.55] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-25_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Cα RMSD fluctuation (Å) | [0.25, 0.55] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
