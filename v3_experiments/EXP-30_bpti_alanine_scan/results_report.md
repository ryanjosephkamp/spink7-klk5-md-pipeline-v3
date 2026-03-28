# EXP-30: BPTI Reactive Region Alanine Scan ΔΔG Panel

## Abstract

This experiment applies Free Energy Perturbation (FEP) to perform computational alanine scanning across 15 residues in the BPTI reactive region (T11–G36) against trypsin. The goal is to compute ΔΔG values for each X→Ala mutation and compare against the comprehensive experimental dataset of Castro & Anderson (1996). The dominant P1 residue K15A (ΔΔG ≈ +10 kcal/mol) and the null-effect G36A (ΔΔG ≈ 0) bracket the dynamic range. Due to the extensive alchemical MD simulations required (~960 ns total), this experiment could not be executed on CPU-only hardware. The CPU result is classified as INCONCLUSIVE.

## Introduction / Background

Alanine scanning mutagenesis is the gold-standard experimental technique for mapping binding hotspots at protein–protein interfaces. For the BPTI–trypsin system, Castro & Anderson (1996) performed the most comprehensive alanine scanning study of any protease–inhibitor system, measuring full kinetic and thermodynamic parameters for 15 mutations spanning the reactive region. Key findings include:

1. **P1 dominance:** K15A increases Ki by ~10⁷-fold (5 × 10⁻¹⁴ → 1.4 × 10⁻⁶ M), contributing ΔΔG ≈ +10 kcal/mol.
2. **Association vs. dissociation:** Most mutations affect k_off more than k_on, indicating that binding affinity is determined primarily by complex stability.
3. **Disulfide contributions:** C14A and C38A each destabilize by ~7 kcal/mol due to disruption of the Cys14–Cys38 disulfide bond.

Computationally reproducing these values via FEP validates the pipeline's ability to predict mutational effects on binding, which is essential for understanding which residues drive binding and for guiding rational design.

**Feature:** F-30 — ΔΔG for 15 interface alanine mutations  
**System:** 2PTC (BPTI–trypsin complex, 1.9 Å resolution)

## Hypothesis

**H₁:** FEP-computed ΔΔG values will achieve a Spearman rank correlation ρ > 0.7 against experimental values from Castro & Anderson (1996).

**H₂:** Per-mutant RMSE < 2.0 kcal/mol across the 15-mutation panel.

**H₃:** ≥80% of mutations will have the correct sign of ΔΔG.

**H₄:** K15A (P1) will be correctly identified as the most destabilizing mutation (ΔΔG > +7 kcal/mol), and G36A as neutral (|ΔΔG| < 1.5 kcal/mol).

## Methods

- **Software:** OpenMM with alchemical FEP via the pipeline's `src.simulate.fep.run_fep_campaign`
- **Force Field:** AMBER ff14SB + TIP3P explicit solvent
- **Protocol:** FEP for each X→Ala mutation using sidechain annihilation (bound complex + free protein legs)
- **Lambda Windows:** 16 evenly spaced windows per mutation (optimized from default 20; Shirts & Chodera 2008)
- **Simulation Time:** 2 ns per lambda window × 16 windows × 2 legs = 64 ns per mutation × 15 mutations = 960 ns total
- **Analysis:** MBAR via `src.analyze.fep.compute_delta_g_mbar` for free energy estimation; thermodynamic cycle closure via `src.analyze.fep.compute_delta_delta_g` for ΔΔG_bind
- **Structure:** PDB 2PTC (BPTI–trypsin complex); optionally loads EXP-04 equilibrated state
- **Temperature:** 310 K

### Mutation Panel

| # | Mutant | Experimental ΔΔG (kcal/mol) | Ki (M) | Category |
|---|--------|----------------------------|--------|----------|
| 1 | T11A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 2 | G12A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 3 | P13A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 4 | C14A | ≈ +7 | ~10⁻⁹ | Disulfide |
| 5 | K15A (P1) | ≈ +10 | 1.4 × 10⁻⁶ | P1 dominant |
| 6 | A16G | varies | — | Framework control |
| 7 | R17A (P2') | ≈ +5 | ~10⁻¹⁰ | Interface contact |
| 8 | I18A | ≈ +3–4 | ~10⁻¹¹ | Hydrophobic packing |
| 9 | I19A | ≈ +3–4 | ~10⁻¹¹ | Hydrophobic packing |
| 10 | R20A | ≈ +2–3 | ~10⁻¹² | Salt bridge |
| 11 | Y21A | ≈ +2–3 | ~10⁻¹² | Aromatic contact |
| 12 | F22A | ≈ +2–3 | ~10⁻¹¹ | Hydrophobic packing |
| 13 | Y23A | ≈ +1–2 | ~10⁻¹³ | Peripheral aromatic |
| 14 | C38A | ≈ +7 | ~10⁻⁹ | Disulfide |
| 15 | G36A | ≈ 0 | ~10⁻¹⁴ | Null control |

## Controls

- **Positive Control:** K15A with experimental ΔΔG ≈ +10 kcal/mol (Castro & Anderson 1996)
- **Negative Control:** G36A with experimental ΔΔG ≈ 0 (null effect)
- **Literature Benchmarks:** Experimental ΔΔG values from Castro & Anderson (1996)
- **Convergence Control:** Forward/reverse overlap analysis and MBAR uncertainty estimates
- **Thermodynamic Cycle Closure:** Comparison of bound and free legs for consistency

## Results

**INCONCLUSIVE** — FEP calculations require GPU-accelerated OpenMM with alchemical protocols. The estimated 960 ns of total simulation time (~82 GPU-hours on A100) makes this infeasible on CPU-only hardware.

## Discussion

Computational alanine scanning via FEP is one of the most informative experiments for validating the pipeline's predictive capability at protein–protein interfaces. The 15-mutation panel provides statistical power for correlation analysis, spanning a dynamic range from null effects (G36A) to massive destabilization (K15A). This experiment is a high priority for GPU execution.

## Conclusions

No conclusions can be drawn. This experiment requires GPU-accelerated alchemical MD simulations that were not feasible on available hardware.

## Figures

No figures generated — full FEP simulation not executed.

## References

1. Castro, M. J. & Anderson, S. (1996). Alanine point-mutations in the reactive region of bovine pancreatic trypsin inhibitor. *Biochemistry*, 35(35), 11435–11446.
2. Shirts, M. R. & Chodera, J. D. (2008). Statistically optimal analysis of samples from multiple equilibrium states. *Journal of Chemical Physics*, 129, 124105.
3. Krowarsch, D., Zakrzewska, M., Skowron, P., & Otlewski, J. (2003). Structure–function relationships in serine protease–inhibitor interactions. *Acta Biochimica Polonica*, 50(2), 367–381.

---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-30_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| Spearman ρ | > 0.7 | (to be filled) | — |
| RMSE | < 2.0 kcal/mol | (to be filled) | — |
| Sign agreement | ≥ 80% | (to be filled) | — |
| K15A most destabilizing | ΔΔG > 7 kcal/mol | (to be filled) | — |
| G36A neutral | \|ΔΔG\| < 1.5 | (to be filled) | — |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)

---

**Author:** Ryan Kamp  
**Affiliation:** Dept. of Computer Science, University of Cincinnati  
**Email:** kamprj@mail.uc.edu  
**GitHub:** ryanjosephkamp
