# EXP-13: Barnase-Barstar Association/Dissociation Kinetics

## Abstract

This experiment aims to compute the association rate constant (k_on), dissociation rate constant (k_off), and binding free energy (ΔG_bind) for the barnase-barstar complex using SMD pulling, umbrella sampling, and Kramers theory applied to the resulting PMF. The system is based on PDB 1BRS (2.0 Å resolution). Experimental benchmarks: k_on = 4 × 10⁸ M⁻¹s⁻¹, k_off = 8 × 10⁻⁶ s⁻¹, ΔG = −19.0 ± 1.0 kcal/mol. **Result: INCONCLUSIVE** — the required ~1500 ns of aggregate MD simulation is infeasible on the available CPU-only platform.

## Introduction/Background

Barnase-barstar is the prototypical electrostatically steered protein-protein association system. Barnase is a ribonuclease from *Bacillus amyloliquefaciens*, and barstar is its intracellular inhibitor. Their association is diffusion-limited (k_on ~ 10⁸–10⁹ M⁻¹s⁻¹), driven by long-range electrostatic complementarity, making this system a benchmark for both thermodynamic and kinetic computational methods.

The crystal structure (PDB: 1BRS, 2.0 Å) has been extensively used in computational studies of protein-protein binding. Pipeline validation on this system (solvation + NPT equilibration) has been completed successfully, confirming that the pipeline can process barnase-barstar through the preparatory stages.

Feature F-13 extends the pipeline validation to the full kinetics computation: PMF → ΔG from the well depth, and k_on/k_off from Kramers theory applied to the barrier heights and curvatures.

## Hypothesis

The SMD/umbrella sampling pipeline applied to barnase-barstar (1BRS), with Kramers theory analysis of the PMF, will yield:
1. ΔG_bind within ±3.0 kcal/mol of −19.0 ± 1.0 kcal/mol
2. k_on within one order of magnitude of 4 × 10⁸ M⁻¹s⁻¹
3. k_off within two orders of magnitude of 8 × 10⁻⁶ s⁻¹

## Methods

### System Preparation
- **Structure**: PDB 1BRS (barnase-barstar complex, 2.0 Å resolution)
- **Force field**: CHARMM36m (protein) + TIP3P (water)
- **Solvation**: Rectangular box with 12 Å padding, 0.15 M NaCl neutralization
- **Status**: Solvation and NPT equilibration already completed successfully

### Equilibration Protocol
1. Energy minimization (steepest descent, 5000 steps) — **completed**
2. NVT equilibration: 500 ps, 300 K — **completed**
3. NPT equilibration: 1 ns, 300 K, 1 bar — **completed**

### Production — SMD Pulling
- 50 independent replicate SMD pulls along the center-of-mass separation vector
- Pulling rate: 0.01 nm/ps over 10 ns per replicate
- Spring constant: 1000 kJ/mol/nm²
- Total SMD simulation: 500 ns

### Production — Umbrella Sampling
- 50 umbrella windows extracted from SMD trajectories
- Window spacing: ~0.5 Å along the reaction coordinate
- 10 ns per window with harmonic biasing potentials
- Total umbrella sampling: 500 ns

### Kinetics Analysis
- PMF reconstruction via WHAM and MBAR
- **ΔG_bind**: Extracted from PMF well depth with volume/orientation corrections
- **k_on**: Kramers theory — diffusion-limited on-rate from PMF barrier height and diffusion coefficient along the reaction coordinate
- **k_off**: Kramers theory — off-rate from barrier height in the bound-to-unbound direction
- **Consistency check**: ΔG = −RT ln(k_on/k_off)
- Bootstrap error estimation (200 bootstrap samples)

### Total Computational Requirement
- ~1500 ns aggregate simulation time (equilibration completed; SMD + umbrella remaining)
- **Estimated CPU time: 15–30 days** on current hardware (MacBook Pro, CPU/OpenCL only)

## Controls

- **Pipeline validation control**: Solvation and equilibration stages completed successfully on 1BRS, confirming structural integrity through the preparatory pipeline.
- **Thermodynamic-kinetic consistency**: ΔG computed from PMF depth must be consistent with ΔG = −RT ln(k_on/k_off) from Kramers-derived rates.
- **Electrostatic steering control**: The strong electrostatic complementarity of barnase-barstar should produce an asymmetric PMF with a shallow long-range attraction funnel.
- **Replicate control**: 50 independent SMD pulls for pathway sampling robustness.

## Results

**INCONCLUSIVE** — experiment requires GPU-accelerated MD simulations not available on the current platform (MacBook Pro, CPU/OpenCL only).

Pipeline validation on 1BRS (solvation + NVT/NPT equilibration) completed successfully, confirming the pipeline can process this system through preparatory stages. The full kinetics computation (SMD pulling + umbrella sampling + Kramers analysis) awaits GPU resources.

The required ~1500 ns of aggregate production simulation would require an estimated 15–30 days of continuous CPU computation. With a single CUDA-enabled GPU (e.g., NVIDIA A100), the same computation would complete in approximately 1–3 days.

## Discussion

Barnase-barstar is uniquely suited for kinetics validation because both the thermodynamic (ΔG) and kinetic (k_on, k_off) parameters are precisely measured. The diffusion-limited association rate provides a stringent test of the PMF shape at long range, where electrostatic steering dominates. The extremely slow dissociation (k_off ~ 10⁻⁶ s⁻¹) reflects a deep, narrow binding well that must be accurately captured by the umbrella sampling.

The successful completion of pipeline validation (solvation + equilibration) on this system is an important milestone: it confirms that the 1BRS structure can be processed through the pipeline without errors, and the equilibrated system is ready for production simulations immediately upon GPU availability.

Kramers theory provides a rigorous framework for extracting rate constants from a 1D PMF, though the reduction to a single reaction coordinate introduces systematic approximations. For barnase-barstar, the center-of-mass separation is expected to be a reasonable reaction coordinate given the funnel-like binding geometry.

**Expected outcome with GPU resources**: ΔG_bind should be reproduced within ±3–5 kcal/mol. Rate constants are more challenging: k_on predictions from PMF-based methods are typically accurate within 1–2 orders of magnitude; k_off predictions are less reliable due to sensitivity to barrier height errors (a 1.4 kcal/mol error in barrier height corresponds to a 10-fold error in rate).

## Conclusions

EXP-13 is **INCONCLUSIVE** due to computational resource constraints. This experiment is notable because pipeline validation (solvation + equilibration) has been completed, and the system is ready for immediate production upon GPU availability. The kinetics analysis (Kramers theory applied to PMF) adds novel methodology beyond the thermodynamic-only experiments (EXP-04/05/06), making this a high-priority target for GPU execution.

## Figures

No figures generated — full simulation not executed due to computational constraints.

## References

1. Schreiber, G. & Fersht, A. R. (1993). Interaction of barnase with its polypeptide inhibitor barstar studied by protein engineering. *Biochemistry*, 32, 5145–5150.
2. Buckle, A. M., Schreiber, G. & Fersht, A. R. (1994). Protein-protein recognition: Crystal structural analysis of a barnase-barstar complex at 2.0 Å resolution. *Biochemistry*, 33, 8878–8889. (PDB: 1BRS)
3. Gabdoulline, R. R. & Wade, R. C. (2002). Biomolecular diffusional association. *Curr. Opin. Struct. Biol.*, 12, 204–213.
4. Kramers, H. A. (1940). Brownian motion in a field of force and the diffusion model of chemical reactions. *Physica*, 7, 284–304.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-13_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| k_on agreement | within 2 orders of experimental | (observed) | (status) |
| Electrostatic enhancement | >100× | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-13_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| k_on agreement | within 2 orders of experimental | (observed) | (status) |
| Electrostatic enhancement | >100× | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
