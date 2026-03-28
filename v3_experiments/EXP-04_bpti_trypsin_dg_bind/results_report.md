# EXP-04: BPTI-Trypsin Binding Free Energy (ΔG_bind)

## Abstract

This experiment aims to compute the absolute binding free energy (ΔG_bind) of the BPTI-trypsin complex using steered molecular dynamics (SMD) pulling and umbrella sampling with WHAM/MBAR analysis. The system is based on PDB 2PTC (1.9 Å resolution). The experimental benchmark is ΔG_bind = −18.2 ± 1.5 kcal/mol. **Result: INCONCLUSIVE** — the required ~1500 ns of aggregate MD simulation is infeasible on the available CPU-only platform within a reasonable timeframe.

## Introduction/Background

BPTI (bovine pancreatic trypsin inhibitor) is the canonical Kunitz-type serine protease inhibitor and one of the most extensively characterized protein-protein interactions in biophysics. The BPTI-trypsin complex (PDB: 2PTC) has been resolved at 1.9 Å and serves as a gold-standard benchmark for computational binding free energy methods. The extraordinarily tight binding (Ki = 6 × 10⁻¹⁴ M) makes this system a rigorous test of any free energy pipeline.

Feature F-04 targets the reproduction of the experimentally measured ΔG_bind using the SMD → umbrella sampling → PMF → ΔG computational pathway implemented in this project's pipeline.

## Hypothesis

The SMD/umbrella sampling pipeline, applied to the BPTI-trypsin system (2PTC), will yield a computed ΔG_bind within ±3.0 kcal/mol of the experimental value of −18.2 ± 1.5 kcal/mol (Ardelt & Bhattacharyya 1988), corresponding to a sub-picomolar inhibition constant.

## Methods

### System Preparation
- **Structure**: PDB 2PTC (BPTI-trypsin complex, 1.9 Å resolution)
- **Force field**: CHARMM36m (protein) + TIP3P (water)
- **Solvation**: Rectangular box with 12 Å padding, 0.15 M NaCl neutralization

### Equilibration Protocol
1. Energy minimization (steepest descent, 5000 steps)
2. NVT equilibration: 500 ps, 300 K, position restraints on heavy atoms
3. NPT equilibration: 1 ns, 300 K, 1 bar, Parrinello-Rahman barostat

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

### Analysis
- WHAM and MBAR for PMF reconstruction
- Bootstrap error estimation (200 bootstrap samples)
- ΔG_bind extraction from PMF with volume/orientation corrections

### Total Computational Requirement
- ~1500 ns aggregate simulation time
- **Estimated CPU time: 15–30 days** on current hardware (MacBook Pro, CPU/OpenCL only)

## Controls

- **Positive control**: BPTI-trypsin is among the tightest known protein-protein interactions; the PMF should show a deep, well-defined minimum.
- **Convergence control**: Block averaging of umbrella windows to assess sampling convergence.
- **Force field control**: CHARMM36m is well-validated for protein-protein interactions.
- **Replicate control**: 50 independent SMD pulls provide statistical robustness for initial pathway sampling.

## Results

**INCONCLUSIVE** — experiment requires GPU-accelerated MD simulations not available on the current platform (MacBook Pro, CPU/OpenCL only).

The required ~1500 ns of aggregate simulation (50 SMD pulls × 10 ns + 50 umbrella windows × 10 ns + equilibration) would require an estimated 15–30 days of continuous computation on CPU, making execution infeasible within the project timeframe. With a single CUDA-enabled GPU (e.g., NVIDIA A100), the same computation would complete in approximately 1–3 days.

## Discussion

The BPTI-trypsin system represents a demanding but well-characterized benchmark for binding free energy calculations. The deep binding well (ΔG ≈ −18 kcal/mol) requires extensive sampling along the unbinding pathway to capture the full PMF profile. The SMD/umbrella sampling approach is well-established for this class of problem, and the pipeline infrastructure (solvation, equilibration, SMD, umbrella sampling, WHAM/MBAR) has been validated on smaller test systems.

The primary limitation is purely computational: the aggregate simulation time of ~1500 ns exceeds what is feasible on a CPU-only platform. This is not a methodological limitation — the pipeline is fully implemented and ready for execution on GPU-equipped hardware.

**Expected outcome with GPU resources**: Based on published studies using similar protocols (e.g., Woo & Roux 2005; Gumbart et al. 2013), we expect the computed ΔG_bind to fall within ±3–5 kcal/mol of the experimental value, with the primary sources of error being force field accuracy and sampling completeness.

## Conclusions

EXP-04 is **INCONCLUSIVE** due to computational resource constraints. The experiment is fully designed, the pipeline is validated, and all input structures are prepared. Execution awaits access to GPU-accelerated hardware, which would reduce the computation time from weeks to days. This experiment remains a high-priority target for future execution.

## Figures

No figures generated — full simulation not executed due to computational constraints.

## References

1. Ardelt, W. & Bhattacharyya, A. (1988). Thermodynamics and kinetics of the association of bovine pancreatic trypsin inhibitor with trypsin. *J. Biol. Chem.*, 263, 3379–3385.
2. Marquart, M., Walter, J., Deisenhofer, J., Bode, W. & Huber, R. (1983). The geometry of the reactive site and of the peptide groups in trypsin, trypsinogen and its complexes with inhibitors. *Acta Crystallogr. B*, 39, 480–490. (PDB: 2PTC)
3. Woo, H.-J. & Roux, B. (2005). Calculation of absolute protein-ligand binding free energy from computer simulations. *Proc. Natl. Acad. Sci. USA*, 102, 6825–6830.
4. Gumbart, J. C., Roux, B. & Chipot, C. (2013). Standard binding free energies from computer simulations: What is the best strategy? *J. Chem. Theory Comput.*, 9, 794–802.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-04_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−24.7, −11.3] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-04_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−24.7, −11.3] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
