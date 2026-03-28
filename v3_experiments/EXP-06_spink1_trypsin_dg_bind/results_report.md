# EXP-06: SPINK1-Trypsin Binding Free Energy (ΔG_bind)

## Abstract

This experiment aims to compute the absolute binding free energy (ΔG_bind) of the SPINK1-trypsin complex using steered molecular dynamics (SMD) and umbrella sampling with WHAM/MBAR analysis. The structural proxy is PDB 1TGS (PSTI, the porcine homolog of human SPINK1, ~70% sequence identity). The experimental benchmark is ΔG_bind = −16.0 ± 2.0 kcal/mol. **Result: INCONCLUSIVE** — the required ~1500 ns of aggregate MD simulation is infeasible on the available CPU-only platform.

## Introduction/Background

SPINK1 (serine peptidase inhibitor, Kazal type 1) is the human pancreatic secretory trypsin inhibitor, a key physiological regulator of premature trypsin activation in the pancreas. Mutations in SPINK1 (notably N34S) are strongly associated with chronic pancreatitis, making this system of direct biomedical relevance.

No crystal structure of the human SPINK1-trypsin complex is available. PDB 1TGS (PSTI-trypsinogen, 1.8 Å) serves as the structural proxy; PSTI is the porcine homolog of SPINK1 with ~70% sequence identity and a conserved reactive loop. The experimental Ki = 2 × 10⁻¹² M for SPINK1-trypsin (Pubols et al. 1974) corresponds to ΔG_bind ≈ −16.0 kcal/mol, intermediate between BPTI-trypsin (EXP-04) and PSTI-chymotrypsinogen (EXP-05).

Feature F-06 targets this binding free energy to complete the three-system validation triangle for the SMD/umbrella sampling pipeline.

## Hypothesis

The SMD/umbrella sampling pipeline applied to the 1TGS-based SPINK1-trypsin structural model will yield a computed ΔG_bind within ±3.0 kcal/mol of the experimental value of −16.0 ± 2.0 kcal/mol (Pubols et al. 1974), despite the use of the porcine PSTI structure as a proxy for human SPINK1.

## Methods

### System Preparation
- **Structure**: PDB 1TGS (PSTI-trypsinogen complex, 1.8 Å resolution), used as structural proxy for SPINK1-trypsin
- **Homology note**: PSTI (porcine) → SPINK1 (human), ~70% sequence identity, conserved Kazal-type reactive loop
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

- **Homology control**: Comparison of computed ΔG_bind with EXP-05 (same 1TGS structure, different benchmark) to assess the impact of using a porcine proxy for the human system.
- **Cross-system control**: Three-way comparison with EXP-04 (BPTI-trypsin) and EXP-05 (PSTI-chymotrypsinogen) to validate pipeline accuracy across binding affinities spanning −14.5 to −18.2 kcal/mol.
- **Convergence control**: Block averaging of umbrella windows.
- **Replicate control**: 50 independent SMD pulls.

## Results

**INCONCLUSIVE** — experiment requires GPU-accelerated MD simulations not available on the current platform (MacBook Pro, CPU/OpenCL only).

The required ~1500 ns of aggregate simulation (50 SMD pulls × 10 ns + 50 umbrella windows × 10 ns + equilibration) would require an estimated 15–30 days of continuous computation on CPU. With a single CUDA-enabled GPU (e.g., NVIDIA A100), the same computation would complete in approximately 1–3 days.

## Discussion

The SPINK1-trypsin system is the most biomedically relevant of the three binding free energy benchmarks (EXP-04, -05, -06), given the direct link between SPINK1 mutations and chronic pancreatitis. The use of the porcine PSTI structure (1TGS) as a proxy introduces a systematic uncertainty: the ~30% sequence divergence may affect peripheral contacts while the reactive loop interaction (which dominates the binding energy) is highly conserved.

The expected ΔG_bind of −16.0 kcal/mol places this system between BPTI-trypsin (−18.2) and PSTI-chymotrypsinogen (−14.5), providing a critical intermediate data point for assessing whether the pipeline can resolve quantitative differences in binding affinity across related systems.

**Expected outcome with GPU resources**: We anticipate the computed ΔG_bind to fall within ±3–5 kcal/mol of the experimental value. The use of a homologous structure rather than the exact complex may introduce an additional ~1–2 kcal/mol systematic error, which should be considered when interpreting results.

## Conclusions

EXP-06 is **INCONCLUSIVE** due to computational resource constraints. This experiment completes the three-system validation set (EXP-04/05/06) and is prioritized for GPU execution alongside the other binding free energy experiments. The homology-based approach (PSTI as proxy for SPINK1) adds interpretive complexity but is well-justified given the high sequence conservation in the reactive loop.

## Figures

No figures generated — full simulation not executed due to computational constraints.

## References

1. Pubols, M. H., Bartelt, D. C. & Greene, L. J. (1974). Trypsin inhibitor from human pancreatic juice. *J. Biol. Chem.*, 249, 2235–2242.
2. Bolognesi, M., Gatti, G., Menegatti, E., Guarneri, M., Marquart, M., Papamokos, E. & Huber, R. (1982). Three-dimensional structure of the complex between pancreatic secretory trypsin inhibitor (Kazal type) and trypsinogen at 1.8 Å resolution. *J. Mol. Biol.*, 162, 839–868. (PDB: 1TGS)
3. Witt, H., Luck, W. & Becker, M. (2000). A signal peptide cleavage site mutation in the cationic trypsinogen gene is strongly associated with chronic pancreatitis. *Gastroenterology*, 117, 7–10.
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
**Notebook:** `EXP-06_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−18.0, −14.0] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)


---

## GPU Results

**Execution Platform:** Google Colab — NVIDIA A100/H100 (to be filled)  
**Execution Date:** (to be filled)  
**Notebook:** `EXP-06_colab.ipynb`  
**Runtime:** (to be filled)

### GPU Quantitative Results

(To be filled after GPU execution)

### GPU Classification

| Criterion | Target | Observed | Status |
|-----------|--------|----------|--------|
| ΔG_bind (kcal/mol) | [−18.0, −14.0] | (observed) | (status) |

**Overall GPU Classification:** (PASS/MARGINAL/FAIL — to be filled)

### GPU Figures

(Figure references to be added after execution)

### GPU Discussion

(To be filled after GPU execution. Compare CPU and GPU results. Discuss convergence, sampling quality, and agreement with experimental benchmarks.)
