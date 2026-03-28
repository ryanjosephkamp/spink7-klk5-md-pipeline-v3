# EXP-05: PSTI-Chymotrypsinogen Binding Free Energy (ΔG_bind)

## Abstract

This experiment aims to compute the absolute binding free energy (ΔG_bind) of the PSTI-chymotrypsinogen complex using steered molecular dynamics (SMD) and umbrella sampling with WHAM/MBAR analysis. The system is based on PDB 1TGS (1.8 Å resolution). The experimental benchmark is ΔG_bind = −14.5 ± 2.0 kcal/mol. **Result: INCONCLUSIVE** — the required ~1500 ns of aggregate MD simulation is infeasible on the available CPU-only platform.

## Introduction/Background

PSTI (pancreatic secretory trypsin inhibitor), also known as SPINK1 in the human homolog, is a Kazal-type serine protease inhibitor. The crystal structure of the PSTI-chymotrypsinogen complex (PDB: 1TGS) at 1.8 Å resolution provides a well-characterized system for binding free energy benchmarking. The inhibition constant (Ki ~ 10⁻¹¹ M) reflects tight but not ultra-tight binding, placing this system in a complementary regime to BPTI-trypsin (EXP-04).

Feature F-05 targets reproduction of the experimentally measured ΔG_bind to validate the pipeline across different inhibitor-protease families (Kazal-type vs. Kunitz-type).

## Hypothesis

The SMD/umbrella sampling pipeline applied to the PSTI-chymotrypsinogen system (1TGS) will yield a computed ΔG_bind within ±3.0 kcal/mol of the experimental value of −14.5 ± 2.0 kcal/mol (Laskowski & Kato 1980).

## Methods

### System Preparation
- **Structure**: PDB 1TGS (PSTI-chymotrypsinogen complex, 1.8 Å resolution)
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

- **Cross-system control**: Comparison with EXP-04 (BPTI-trypsin) and EXP-06 (SPINK1-trypsin) to validate pipeline consistency across inhibitor-protease families.
- **Convergence control**: Block averaging of umbrella windows to assess sampling convergence.
- **Force field control**: CHARMM36m validated for Kazal-type inhibitor systems.
- **Replicate control**: 50 independent SMD pulls for statistical robustness.

## Results

**INCONCLUSIVE** — experiment requires GPU-accelerated MD simulations not available on the current platform (MacBook Pro, CPU/OpenCL only).

The required ~1500 ns of aggregate simulation (50 SMD pulls × 10 ns + 50 umbrella windows × 10 ns + equilibration) would require an estimated 15–30 days of continuous computation on CPU. With a single CUDA-enabled GPU (e.g., NVIDIA A100), the same computation would complete in approximately 1–3 days.

## Discussion

The PSTI-chymotrypsinogen system provides a complementary benchmark to BPTI-trypsin (EXP-04), testing the pipeline on a Kazal-type inhibitor rather than a Kunitz-type. The moderately tight binding (ΔG ≈ −14.5 kcal/mol, ~4 kcal/mol weaker than BPTI-trypsin) probes a different region of the free energy landscape, which is important for validating that the pipeline can resolve differences in binding affinity across systems.

The 1TGS structure at 1.8 Å provides excellent starting coordinates. The PSTI-chymotrypsinogen interface involves a canonical Kazal-type reactive loop insertion into the protease active site, a structurally distinct binding mode from the Kunitz-type BPTI.

**Expected outcome with GPU resources**: The computed ΔG_bind is expected to fall within ±3–5 kcal/mol of the experimental value, consistent with the accuracy typically achieved by SMD/umbrella sampling methods for protein-protein systems.

## Conclusions

EXP-05 is **INCONCLUSIVE** due to computational resource constraints. The experiment design is complete and the pipeline is ready for execution. This system is prioritized for GPU execution alongside EXP-04 and EXP-06 to enable cross-system validation of the binding free energy pipeline.

## Figures

No figures generated — full simulation not executed due to computational constraints.

## References

1. Laskowski, M. Jr. & Kato, I. (1980). Protein inhibitors of proteinases. *Annu. Rev. Biochem.*, 49, 593–626.
2. Bolognesi, M., Gatti, G., Menegatti, E., Guarneri, M., Marquart, M., Papamokos, E. & Huber, R. (1982). Three-dimensional structure of the complex between pancreatic secretory trypsin inhibitor (Kazal type) and trypsinogen at 1.8 Å resolution. *J. Mol. Biol.*, 162, 839–868. (PDB: 1TGS)
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
**Notebook:** `EXP-05_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−16.5, −12.5] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-05_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−16.5, −12.5] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
