# EXP-29: SH3-p41 Binding Validation

## Abstract

This experiment aims to compute the binding free energy (ΔG_bind) of the SH3 domain–p41 peptide complex as a cross-validation test for the SMD/umbrella sampling pipeline on a protein-peptide system. The system is based on PDB 4EIK. The experimental benchmark is ΔG_bind = −8.5 ± 1.5 kcal/mol (estimated from ITC). **Result: INCONCLUSIVE** — the required ~1200 ns of aggregate MD simulation is infeasible on the available CPU-only platform.

## Introduction/Background

SH3 (Src Homology 3) domains are small (~60 residue) protein interaction modules that recognize proline-rich peptide motifs. The SH3-p41 complex (PDB: 4EIK) represents a qualitatively different binding regime from the protease-inhibitor systems tested in EXP-04/05/06: the interaction is weaker (ΔG ≈ −8.5 kcal/mol vs. −14 to −18 kcal/mol), the binding interface is smaller, and the peptide ligand is flexible.

Feature F-29 serves as a cross-validation target to test whether the pipeline, originally developed and validated on tight protein-protein complexes, can also handle moderate-affinity protein-peptide interactions. This tests the generalizability of the computational approach.

## Hypothesis

The SMD/umbrella sampling pipeline applied to the SH3-p41 system (4EIK) will yield a computed ΔG_bind within ±2.0 kcal/mol of the experimental value of −8.5 ± 1.5 kcal/mol (Mayer et al., ITC measurements). The weaker binding and smaller interface should produce a shallower PMF profile compared to the protease-inhibitor systems.

## Methods

### System Preparation
- **Structure**: PDB 4EIK (SH3 domain + p41 peptide)
- **Force field**: CHARMM36m (protein/peptide) + TIP3P (water)
- **Solvation**: Rectangular box with 12 Å padding, 0.15 M NaCl neutralization

### Equilibration Protocol
1. Energy minimization (steepest descent, 5000 steps)
2. NVT equilibration: 500 ps, 300 K, position restraints on heavy atoms
3. NPT equilibration: 1 ns, 300 K, 1 bar, Parrinello-Rahman barostat

### Production — SMD Pulling
- 50 independent replicate SMD pulls along the center-of-mass separation vector
- Pulling rate: 0.01 nm/ps over 8 ns per replicate (shorter system → shorter pulling distance)
- Spring constant: 1000 kJ/mol/nm²
- Total SMD simulation: 400 ns

### Production — Umbrella Sampling
- 40 umbrella windows extracted from SMD trajectories
- Window spacing: ~0.5 Å along the reaction coordinate
- 10 ns per window with harmonic biasing potentials
- Total umbrella sampling: 400 ns

### Analysis
- WHAM and MBAR for PMF reconstruction
- Bootstrap error estimation (200 bootstrap samples)
- ΔG_bind extraction from PMF with volume corrections

### Total Computational Requirement
- ~1200 ns aggregate simulation time
- **Estimated CPU time: 12–25 days** on current hardware (MacBook Pro, CPU/OpenCL only)

## Controls

- **Cross-system control**: Comparison with EXP-04/05/06 (protease-inhibitor systems) to assess pipeline accuracy across different binding affinity regimes.
- **Affinity regime control**: ΔG ≈ −8.5 kcal/mol is ~6–10 kcal/mol weaker than the protease-inhibitor benchmarks, testing the pipeline in a regime where sampling errors are a larger fraction of the signal.
- **Convergence control**: Block averaging of umbrella windows; shorter total pathway (protein-peptide vs. protein-protein) may improve convergence.
- **Replicate control**: 50 independent SMD pulls for robust pathway sampling.

## Results

**INCONCLUSIVE** — experiment requires GPU-accelerated MD simulations not available on the current platform (MacBook Pro, CPU/OpenCL only).

The required ~1200 ns of aggregate simulation (50 SMD pulls × 8 ns + 40 umbrella windows × 10 ns + equilibration) would require an estimated 12–25 days of continuous CPU computation. With a single CUDA-enabled GPU (e.g., NVIDIA A100), the same computation would complete in approximately 1–2 days.

## Discussion

The SH3-p41 system provides a critical cross-validation test for several reasons:

1. **Different affinity regime**: At ΔG ≈ −8.5 kcal/mol, this is a moderate-affinity interaction where the computed PMF must resolve a relatively shallow well. Sampling errors that are acceptable for deep wells (as in BPTI-trypsin) may dominate the signal here.

2. **Protein-peptide vs. protein-protein**: The flexible peptide ligand introduces conformational sampling challenges not present in the rigid protease-inhibitor systems. The reaction coordinate (COM separation) may be less well-defined.

3. **Smaller system size**: The SH3 domain (~60 residues) + p41 peptide is significantly smaller than the protease-inhibitor complexes, which partially compensates for the reduced signal-to-noise by enabling faster per-step computation.

**Expected outcome with GPU resources**: The computed ΔG_bind should fall within ±2–4 kcal/mol of the experimental value. Peptide flexibility may require enhanced sampling techniques (e.g., replica exchange) for convergence, which would increase the computational cost beyond the current estimate.

## Conclusions

EXP-29 is **INCONCLUSIVE** due to computational resource constraints. This cross-validation experiment is important for establishing the generalizability of the pipeline beyond the protease-inhibitor systems that constitute the primary validation set. Execution is prioritized after the core experiments (EXP-04/05/06/13) are completed on GPU hardware.

## Figures

No figures generated — full simulation not executed due to computational constraints.

## References

1. Mayer, B. J. (2001). SH3 domains: complexity in moderation. *J. Cell Sci.*, 114, 1253–1263.
2. PDB: 4EIK — SH3 domain in complex with p41 peptide.
3. Gumbart, J. C., Roux, B. & Chipot, C. (2013). Standard binding free energies from computer simulations: What is the best strategy? *J. Chem. Theory Comput.*, 9, 794–802.
4. Woo, H.-J. & Roux, B. (2005). Calculation of absolute protein-ligand binding free energy from computer simulations. *Proc. Natl. Acad. Sci. USA*, 102, 6825–6830.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-29_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−10.0, −7.0] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-29_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−10.0, −7.0] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
