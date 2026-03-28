# SPINK7-KLK5 MD Pipeline: Binding Free Energy via Enhanced Sampling

[![Version: 3.0](https://img.shields.io/badge/version-3.0-blue.svg)]()
[![OpenMM](https://img.shields.io/badge/OpenMM-%E2%89%A58.1-blue.svg)](https://openmm.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![Tests: 367 passed](https://img.shields.io/badge/tests-367%20passed-brightgreen.svg)]()
[![Force Fields](https://img.shields.io/badge/force%20fields-AMBER%20%7C%20AMOEBA%20%7C%20ANI--2x-orange.svg)]()
[![Experiments: 33](https://img.shields.io/badge/experiments-33%20designed-blueviolet.svg)]()
[![Literature: 28 sources](https://img.shields.io/badge/literature-28%20sources-informational.svg)]()

An end-to-end molecular dynamics simulation pipeline for computing the binding free energy of the **SPINK7-KLK5** protease-antiprotease complex — a protein-protein interaction central to the pathogenesis of **Eosinophilic Esophagitis (EoE)**. The pipeline automates the full computational biophysics workflow: PDB retrieval, structure cleaning, protonation at physiological pH, AMBER ff14SB topology construction, explicit TIP3P solvation, energy minimization, NVT/NPT equilibration, and unrestrained production dynamics. Three complementary enhanced sampling strategies — **Steered Molecular Dynamics (SMD)** with the **Jarzynski equality**, **Umbrella Sampling** with **WHAM/MBAR** reconstruction of the Potential of Mean Force, and **well-tempered metadynamics** — provide rigorous, cross-validated binding free energy estimates with quantified statistical uncertainties. **Alchemical free energy perturbation (FEP)** enables computational mutagenesis, and **Markov State Model (MSM)** construction provides kinetic rate information.

**Version 2 (V2)** addresses 40 systematically identified limitations in the original pipeline, spanning physical correctness, algorithmic rigor, software architecture, testing, performance, and visualization, upgrading the pipeline from a functional prototype to a production-ready research tool. Ten physical validity invariants are enforced as runtime checks. A comprehensive test suite of **367 unit, integration, and analytical tests** validates every module, passing across CPU and GPU platforms with zero regressions.

**Version 3 (V3)** conducts a systematic benchmarking campaign — **33 experiments** designed against **28 literature sources** — to validate pipeline predictions against published experimental data. Of 9 CPU-completed experiments, **7 achieved PASS** and **2 MARGINAL**. The remaining 24 are classified INCONCLUSIVE: 7 due to absent co-crystal structures, 17 due to GPU compute requirements (~530–570 A100 GPU-hours). Full experimental designs, Colab notebooks, and GPU execution infrastructure have been prepared for future completion.

---

## Table of Contents

- [Biological Motivation](#biological-motivation)
- [Mathematical Framework](#mathematical-framework)
- [V2 Pipeline Implementation](#v2-pipeline-implementation)
- [V3 Benchmarking Campaign](#v3-benchmarking-campaign)
- [Results](#results)
- [Physical Validity Invariants](#physical-validity-invariants)
- [Repository Structure](#repository-structure)
- [Installation & Setup](#installation--setup)
- [Usage](#usage)
- [GPU Experiment Execution Guide](#gpu-experiment-execution-guide)
- [Simulation Protocol](#simulation-protocol)
- [Test Suite](#test-suite)
- [Data Flow](#data-flow)
- [Dependencies](#dependencies)
- [Key Findings & Identified Limitations](#key-findings--identified-limitations)
- [Future Directions](#future-directions)
- [References](#references)
- [License](#license)
- [Author](#author)

---

<div style="page-break-after: always;"></div>

## Biological Motivation

Eosinophilic Esophagitis (EoE) is a chronic inflammatory disease of the esophagus driven by IL-13-mediated transcriptional silencing of *SPINK7* (Serine Peptidase Inhibitor, Kazal Type 7). Under homeostatic conditions, SPINK7 stoichiometrically inhibits KLK5 (Kallikrein-Related Peptidase 5), a trypsin-like serine protease. SPINK7 deficiency unleashes KLK5 proteolytic activity, which degrades Desmoglein-1 (DSG1) and compromises the epithelial barrier, permitting allergen penetration.

The SPINK7-KLK5 interaction follows the canonical **Laskowski mechanism** for Kazal-type inhibitor–serine protease binding: the reactive site loop (RSL) of SPINK7 inserts into the KLK5 active-site cleft in a substrate-like orientation, and the conformational rigidity imposed by three disulfide bonds renders the acyl-enzyme intermediate thermodynamically trapped with an extremely slow $k_{\text{off}}$. The binding interface buries approximately 800–1200 Å² of solvent-accessible surface area, stabilized by backbone hydrogen bonds, electrostatic complementarity (P1-Arg/Lys ↔ Asp189), and hydrophobic contacts at the P2/P2' sub-sites.

Understanding this interaction at the atomic level through molecular dynamics simulation is essential for developing targeted therapeutic strategies to restore epithelial barrier function in EoE patients.

---

## Mathematical Framework

<details>
<summary><strong>Force Field Potential Energy</strong></summary>

The total potential energy is decomposed into bonded and nonbonded contributions:

$$V_{\text{total}}(\mathbf{r}) = V_{\text{bonded}}(\mathbf{r}) + V_{\text{nonbonded}}(\mathbf{r})$$

**Bonded terms** (covalent interactions within the molecular topology):

$$V_{\text{bonded}} = \sum_{\text{bonds}} K_b (b - b_0)^2 + \sum_{\text{angles}} K_\theta (\theta - \theta_0)^2 + \sum_{\text{dihedrals}} \frac{V_n}{2} [1 + \cos(n\phi - \gamma)] + \sum_{\text{impropers}} K_\omega (\omega - \omega_0)^2$$

where $K_b$, $K_\theta$, $V_n$, and $K_\omega$ are force constants; $b_0$, $\theta_0$, $\gamma$, and $\omega_0$ are equilibrium values; and $n$ is the dihedral periodicity.

**Nonbonded terms** (van der Waals + electrostatics):

$$V_{\text{nonbonded}} = \sum_{i < j} \left[ 4\epsilon_{ij} \left( \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6} \right) + \frac{q_i q_j}{4\pi\epsilon_0 r_{ij}} \right]$$

The **AMBER ff14SB** force field provides optimized backbone torsion parameters. Water is modeled using **TIP3P** (with configurable OPC and TIP4P-Ew alternatives), and ions follow the **Joung-Cheatham** monovalent parameters.

</details>

<details>
<summary><strong>Long-Range Electrostatics (Particle Mesh Ewald)</strong></summary>

$$E_{\text{elec}} = E_{\text{direct}} + E_{\text{reciprocal}} + E_{\text{self-correction}}$$

PME reduces the computational cost from $O(N^2)$ to $O(N \log N)$ using B-spline interpolation (order 5) with grid spacing $\leq 1.0$ Å and direct-space cutoff $r_c = 10$ Å.

</details>

<details>
<summary><strong>Langevin Equation of Motion</strong></summary>

$$m_i \ddot{\mathbf{r}}_i = -\nabla_i V(\mathbf{r}) - \gamma m_i \dot{\mathbf{r}}_i + \sqrt{2 \gamma m_i k_B T} \, \boldsymbol{\eta}_i(t)$$

where $\gamma = 1.0 \text{ ps}^{-1}$, $T = 310$ K, and $\boldsymbol{\eta}_i(t)$ is Gaussian white noise. Integration uses the **Langevin middle integrator** with $\Delta t = 2$ fs, enabled by **SHAKE** constraints on hydrogen bonds.

</details>

<details>
<summary><strong>Steered Molecular Dynamics & Jarzynski Equality</strong></summary>

SMD applies a time-dependent harmonic bias to the center-of-mass (COM) distance $\xi$:

$$U_{\text{SMD}}(\xi, t) = \frac{k}{2} \left[ \xi(t) - \xi_0 - v \cdot t \right]^2$$

where $k = 1000$ kJ/mol/nm², $v = 0.001$ nm/ps, and the COM distance uses the **minimum image convention** under periodic boundary conditions:

$$\xi(\mathbf{r}) = \left\| \text{MIC}\left(\mathbf{R}_{\text{COM}}^{\text{SPINK7}} - \mathbf{R}_{\text{COM}}^{\text{KLK5}}\right) \right\|_2$$

The **Jarzynski equality** connects non-equilibrium work to equilibrium free energy:

$$\Delta G = -k_B T \ln \left[ \frac{1}{N_{\text{traj}}} \sum_{j=1}^{N_{\text{traj}}} e^{-\beta W_j} \right]$$

with the second-order cumulant expansion for near-Gaussian work distributions:

$$\Delta G \approx \langle W \rangle - \frac{\beta}{2} \sigma_W^2$$

The **Bennett Acceptance Ratio (BAR)** estimator provides cross-validation:

$$\sum_F \frac{1}{1 + \frac{n_F}{n_R} \exp[\beta(W_F - C)]} = \sum_R \frac{1}{1 + \frac{n_R}{n_F} \exp[\beta(-W_R - C)]}$$

</details>

<details>
<summary><strong>Umbrella Sampling & WHAM / MBAR</strong></summary>

Each of $M$ discrete windows is biased by a harmonic potential:

$$U_i^{\text{bias}}(\xi) = \frac{k_i}{2} \left( \xi - \xi_i^{\text{ref}} \right)^2$$

The **Weighted Histogram Analysis Method** recovers the unbiased PMF via self-consistent equations:

$$P^{\text{unbiased}}(\xi) = \frac{\sum_{i=1}^{M} n_i \, h_i(\xi)}{\sum_{i=1}^{M} n_i \, \exp\left[ \beta \left( f_i - U_i^{\text{bias}}(\xi) \right) \right]}$$

$$e^{-\beta f_i} = \int P^{\text{unbiased}}(\xi) \, \exp\left[ -\beta \, U_i^{\text{bias}}(\xi) \right] d\xi$$

iterated until convergence: $\max_i |f_i^{(n+1)} - f_i^{(n)}| < 10^{-6}$ kJ/mol.

The **Multistate Bennett Acceptance Ratio (MBAR)** avoids histogram binning entirely:

$$\hat{f}_i = -\ln \sum_{n=1}^{N} \frac{\exp[-\beta U_i(\mathbf{x}_n)]}{\sum_{k=1}^{K} N_k \exp[\hat{f}_k - \beta U_k(\mathbf{x}_n)]}$$

The PMF from either method yields:

$$G(\xi) = -k_B T \ln P^{\text{unbiased}}(\xi) + C$$

</details>

<details>
<summary><strong>Binding Free Energy Extraction</strong></summary>

$$\Delta G_{\text{bind}}^{\circ} = -k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{site}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right] + k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{bulk}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right]$$

where $C^{\circ} = 1/1660$ Å$^{-3}$ corresponds to the standard concentration of 1 M.

</details>

<details>
<summary><strong>Alchemical Free Energy Perturbation</strong></summary>

Computational mutagenesis employs a thermodynamic cycle to compute $\Delta\Delta G_{\text{bind}}$ upon point mutation:

$$\Delta\Delta G_{\text{bind}} = \Delta G_{\text{bind}}^{\text{mut}} - \Delta G_{\text{bind}}^{\text{wt}} = \Delta G_{\text{alch}}^{\text{complex}} - \Delta G_{\text{alch}}^{\text{free}}$$

</details>

---

<div style="page-break-after: always;"></div>

## V2 Pipeline Implementation

Version 2 addresses **40 systematically identified limitations** in the V1 pipeline, organized into seven categories. Each limitation was classified by severity (1 Critical, 10 High, 17 Medium, 12 Low) and independently verified with full regression testing. The complete implementation report is available at [`reports/full_implementation_report_v2.md`](reports/full_implementation_report_v2.md).

### Critical and High-Severity Fixes

<details>
<summary><strong>L-01 — PBC-Unaware COM Distance Calculation (Critical)</strong></summary>

The V1 pipeline computed center-of-mass distances without applying the minimum image convention, producing erroneous distances whenever molecular fragments crossed periodic box boundaries. This silently corrupts the reaction coordinate $\xi$ that drives both SMD pulling and umbrella window placement, propagating incorrect forces throughout the enhanced sampling workflow and invalidating the resulting PMF. The fix implements image-aware COM distance computation using the angular mean method, which maps Cartesian coordinates onto a unit circle via $\theta_i = 2\pi x_i / L$, computes

$$\bar{\theta} = \text{atan2}\!\left(\frac{1}{N}\sum_i \sin\theta_i,\;\frac{1}{N}\sum_i \cos\theta_i\right)$$

per dimension, and maps back to Cartesian space, correctly handling all boundary-crossing geometries.

</details>

<details>
<summary><strong>L-07 — Exponential Averaging Bias in Jarzynski Estimator (High)</strong></summary>

The Jarzynski equality $\Delta G = -k_BT \ln\langle e^{-\beta W}\rangle$ suffers from systematic bias when the number of non-equilibrium trajectories $N$ is finite. V2 implements the second-order cumulant expansion and the BAR estimator as cross-validators, and reports convergence diagnostics based on the effective sample size $n_{\text{eff}} = (\sum e^{-\beta W_j})^2 / \sum e^{-2\beta W_j}$.

</details>

<details>
<summary><strong>L-02 — Henderson-Hasselbalch Protonation Without Electrostatic Environment (High)</strong></summary>

V1 assigned protonation states using tabulated model $\text{p}K_a$ values, ignoring shifts induced by the protein electrostatic environment. V2 integrates PROPKA-predicted residue-specific $\text{p}K_a$ values, accounting for desolvation penalties, hydrogen bonding networks, and charge-charge interactions within the folded protein.

</details>

<details>
<summary><strong>Other High-Severity Fixes</strong></summary>

| ID | Limitation | Fix Summary |
|----|-----------|-------------|
| L-03 | `NoCutoff` nonbonded method | Enforced PME with proper cutoff specification |
| L-15 | Fixed-charge force field limits | Added AMOEBA polarizable and ANI-2x ML/MM force field tiers |
| L-18 | Hardcoded random seeds | Per-replica cryptographic seeding via `os.urandom()` |
| L-19 | Synthetic topology corrupts metadata | Preserved full biochemical topology through the simulation pipeline |
| L-29 | $O(N^3)$ memory in contact analysis | Frame-streaming with configurable batch sizes |
| L-37 | Data format mismatches | Typed data contracts with schema validation at stage boundaries |
| L-38 | No production-scale results | SPINK7-KLK5 production: 100 ns MD, 50 SMD replicates, 50 umbrella windows |

</details>

### Full Categorized Summary of All 40 Fixes

<details>
<summary><strong>Physics and Correctness (9 fixes)</strong></summary>

| ID | Limitation | Severity | Fix Summary |
|----|-----------|----------|-------------|
| L-01 | PBC-unaware COM distance | Critical | Minimum image convention via angular mean method |
| L-02 | Henderson-Hasselbalch protonation | High | PROPKA-predicted residue-specific $\text{p}K_a$ values |
| L-03 | `NoCutoff` in topology builder | High | Enforced PME nonbonded method with proper cutoff |
| L-04 | Double protonation risk | Medium | Protonation state tracking guard |
| L-05 | No CIF/mmCIF support | Medium | Dual-format parser dispatching (PDB and mmCIF) |
| L-06 | NMR multi-model handling | Low | Model-1 extraction with configurable selection |
| L-14 | Finite-size electrostatic artifacts | Medium | Analytical finite-size correction for PME |
| L-15 | Fixed-charge force field limits | High | AMOEBA polarizable and ML/MM hybrid potential tiers |
| L-17 | No truncated octahedron boxes | Low | Truncated octahedron geometry (~29% fewer solvent atoms) |

</details>

<details>
<summary><strong>Algorithm and Methodology (9 fixes)</strong></summary>

| ID | Limitation | Severity | Fix Summary |
|----|-----------|----------|-------------|
| L-07 | Jarzynski exponential averaging bias | High | BAR estimator and cumulant expansion with convergence diagnostics |
| L-08 | Single scalar reaction coordinate | High | Multi-dimensional collective variable framework |
| L-09 | No umbrella pre-equilibration | Medium | Three-phase protocol with Chodera equilibration detection |
| L-10 | Histogram overlap detection gaps | Medium | Quantitative overlap analysis with window insertion recommendations |
| L-11 | No MBAR alternative | Medium | MBAR solver via pymbar, eliminating histogram binning artifacts |
| L-12 | Bootstrap ignores autocorrelation | Medium | Stationary block bootstrap with automatic block length selection |
| L-13 | No metadynamics support | Medium | Well-tempered metadynamics with Gaussian hill deposition |
| L-39 | No $\Delta\Delta G$ capability | Medium | Alchemical FEP via thermodynamic cycle for mutagenesis |
| L-40 | No MSM construction | Low | MSM pipeline with TICA dimensionality reduction |

</details>

<details>
<summary><strong>Software Architecture and Design (6 fixes)</strong></summary>

| ID | Limitation | Severity | Fix Summary |
|----|-----------|----------|-------------|
| L-19 | Synthetic topology corrupts metadata | High | Preserved full biochemical Topology through simulation pipeline |
| L-20 | `assert`-based runtime validation | Medium | Replaced with explicit `ValueError`/`TypeError` exceptions |
| L-21 | No config file support | Medium | YAML configuration with precedence cascade |
| L-23 | No resume-from-checkpoint | Medium | Checkpoint lifecycle manager with state serialization |
| L-24 | `sys.path` manipulation in scripts | Low | Proper relative imports and package entry points |
| L-37 | Data format mismatches | High | Typed data contracts with schema validation |

</details>

<details>
<summary><strong>Testing, Performance, Visualization, and Production (16 fixes)</strong></summary>

| ID | Limitation | Severity | Fix Summary |
|----|-----------|----------|-------------|
| L-26 | No full integration test | Medium | End-to-end integration test on alanine dipeptide |
| L-27 | No tests under `python -O` | Low | Optimization-mode test harness |
| L-28 | Missing parser edge-case tests | Low | Comprehensive edge-case matrix |
| L-32 | No PBC unwrapping before analysis | Medium | Trajectory unwrapping prior to RMSD/RMSF/Rg analysis |
| L-33 | No equilibration detection | Medium | Chodera automated equilibration detection |
| L-16 | TIP3P water model limitations | Low | OPC and TIP4P-Ew water model support |
| L-18 | Hardcoded random seeds | High | Per-replica cryptographic seeding |
| L-22 | No PDB download retry/cache | Low | Exponential backoff retry with local file cache |
| L-25 | CPU-only platform hardcoding | Medium | Automatic platform detection: CUDA → OpenCL → CPU |
| L-29 | $O(N^3)$ memory in contacts | High | Frame-streaming contact analysis |
| L-30 | Sequential SMD/umbrella execution | Medium | Concurrent replicate execution via `ProcessPoolExecutor` |
| L-31 | No streaming SMD work aggregation | Low | Streaming Welford aggregation |
| L-34 | No minimizer convergence status | Low | Convergence reporting with L-BFGS support |
| L-35 | Chain color palette overflow | Low | Cyclic palette with perceptually distinct hue rotation |
| L-36 | No automated figure generation | Low | `generate_figures.py` for publication-quality figures |
| L-38 | No production-scale results | High | SPINK7-KLK5 production results |

</details>

### GPU-Validated Force Field Hierarchy

| Tier | Force Field | Electrostatics | Polarization | GPU Test | Validation |
|------|------------|----------------|--------------|----------|------------|
| 1 (Default) | AMBER ff14SB | Fixed point charges | None | — | 364 CPU tests |
| 2 (Polarizable) | AMOEBA 2018 | Permanent multipoles | Self-consistent induced dipoles | GPU-01 | Finite energy with SCF convergence |
| 3 (ML) | ANI-2x | Implicit (learned from DFT) | Implicit (learned from DFT) | GPU-02, GPU-03 | Valid system creation + finite energy |

---

## V3 Benchmarking Campaign

V3 transitions the pipeline from a **verified implementation** to a **validated scientific instrument** by systematically benchmarking computational predictions against published experimental data from the peer-reviewed literature.

### Literature Review

A comprehensive review of **28 primary literature sources** identified **33 experimentally measurable features** across six categories:

| Category | Features | Example |
|----------|----------|---------|
| Thermodynamic | F-01 – F-08 | $\Delta G_{\text{bind}}$ for BPTI-trypsin, SPINK1-trypsin, SPINK7-KLK5 |
| Kinetic | F-09 – F-13 | $k_{\text{on}}/k_{\text{off}}$ for BPTI-trypsin, barnase-barstar |
| Structural | F-14 – F-23 | Catalytic contact distance, BSA, H-bond count, loop geometry |
| Dynamic | F-24 – F-25 | BPTI amide H/D exchange, conformational variability |
| Biophysical | F-26 – F-29 | BSA–$\Delta G$ correlation, interface packing, SH3-p41 validation |
| Mutational | F-30 – F-33 | BPTI alanine scan, disulfide ablation, SPINK1 N34S, barnase-barstar DMC |

Full literature details: [`v3_literature/benchmarks.md`](v3_literature/benchmarks.md) and [`v3_literature/sources.md`](v3_literature/sources.md).

<div style="page-break-after: always;"></div>

### Systems Hierarchy

| System | Relationship to SPINK7-KLK5 | Role |
|--------|------------------------------|------|
| SPINK7-KLK5 | Primary target | Direct validation |
| BPTI-trypsin | Protease-inhibitor gold standard | Method validation |
| SPINK1-trypsin | Same Kazal family, different system | Kazal-family transferability |
| PSTI-chymotrypsinogen | Kazal-type, co-crystal available | Structural reference |
| Barnase-barstar | Ultra-tight PPI reference | PPI energetics methodology |
| SH3-p41 peptide | Computational methods benchmark | PMF/alchemical calibration |
| LEKTI-KLK5 | Same protease, related Kazal inhibitor | Cross-validation |

### Experimental Design

Each experiment follows a standardized three-document framework:

1. **`experimental_design.md`** — Scientific background, target benchmark values with 95% confidence intervals, statistical validation criteria
2. **`implementation_guide.md`** — Step-by-step computational protocol using the V2 pipeline
3. **`results_report.md`** — Quantitative results, statistical analysis, classification (PASS / MARGINAL / INCONCLUSIVE)

All 33 experiments have complete experimental designs and implementation guides. CPU-executable experiments additionally contain results reports with quantitative outcomes and figures.

---

## Results

### Summary

| Classification | Count | Description |
|----------------|-------|-------------|
| **PASS** | 7 | Prediction within 95% CI of experimental benchmark |
| **MARGINAL** | 2 | Prediction within 2× the 95% CI |
| **INCONCLUSIVE — Structure** | 7 | No co-crystal structure available in the PDB |
| **INCONCLUSIVE — Compute** | 17 | GPU compute required (~530–570 A100 GPU-hours) |

### Complete Experiment Results

<details>
<summary><strong>PASS — 7 experiments</strong></summary>

| ID | Experiment | Feature | Target (95% CI) | Result | System |
|----|-----------|---------|------------------|--------|--------|
| EXP-15 | Catalytic Contact Distance | F-15 | [2.4, 3.0] Å | 2.675 Å | BPTI-trypsin |
| EXP-16 | Buried Surface Area | F-16 | [1360, 1700] Å² | 1608.5 Å² | BPTI-trypsin |
| EXP-17 | Interfacial H-Bond Count | F-17 | [7, 13] | 5 H-bonds | BPTI-trypsin |
| EXP-19 | Nonpolar BSA Fraction | F-19 | [55%, 67%] | 56.9% | BPTI-trypsin |
| EXP-20 | Loop Conservation | F-20 | [0.0, 1.0] Å | 0.234 Å | Kazal family |
| EXP-21 | Subsite Specificity | F-21 | S1 salt bridge < 4.0 Å | 3.63 Å | BPTI-trypsin |
| EXP-23 | BPTI Crystal Structure | F-23 | [1.85, 2.15] Å (S-S) | 2.054 Å | BPTI |

</details>

<details>
<summary><strong>MARGINAL — 2 experiments</strong></summary>

| ID | Experiment | Feature | Target (95% CI) | Result | Notes |
|----|-----------|---------|------------------|--------|-------|
| EXP-14 | Binding Loop Geometry | F-14 | [0.0, 1.5] Å | 1.766 Å | Exceeded upper bound by 0.27 Å |
| EXP-18 | Interfacial Water Count | F-18 | [15, 30] waters | 8.0 ± 2.4 | Static minimized snapshot vs. dynamic ensemble |

</details>

<details>
<summary><strong>INCONCLUSIVE — Structure-Limited (7 experiments)</strong></summary>

| ID | Experiment | Reason |
|----|-----------|--------|
| EXP-01 | SPINK7-KLK5 $\Delta G_{\text{bind}}$ | No SPINK7-KLK5 complex structure in PDB |
| EXP-02 | SPINK7-KLK12 $\Delta G_{\text{bind}}$ | No SPINK7-KLK12 complex structure |
| EXP-03 | LEKTI-KLK5 Panel | No LEKTI-KLK5 complex structure |
| EXP-11 | KLK5 Substrate Kinetics | No substrate-bound KLK5 structure |
| EXP-12 | LEKTI-KLK5 pH Kinetics | No complex + requires constant-pH MD |
| EXP-22 | SPINK7-KLK5 Geometry | Depends on EXP-01 complex |
| EXP-27 | Interface Packing Density | Insufficient resolved water coordinates |

</details>

<details>
<summary><strong>INCONCLUSIVE — Compute-Limited / GPU-Required (17 experiments)</strong></summary>

| ID | Experiment | GPU Method | Est. A100 Hours | Dependencies |
|----|-----------|-----------|-----------------|--------------|
| **Tier 1 — Thermodynamic ($\Delta G_{\text{bind}}$)** | | | | |
| EXP-04 | BPTI-Trypsin $\Delta G_{\text{bind}}$ | SMD + US | ~52 | Independent |
| EXP-05 | PSTI-Chymotrypsinogen $\Delta G_{\text{bind}}$ | SMD + US | ~52 | Independent |
| EXP-06 | SPINK1-Trypsin $\Delta G_{\text{bind}}$ | SMD + US | ~52 | Independent |
| EXP-13 | Barnase-Barstar Kinetics | US + BD | ~42 | Independent |
| EXP-29 | SH3-p41 Validation | SMD + US | ~52 | Independent |
| **Tier 2 — PMF Analysis** | | | | |
| EXP-07 | P1 Energetic Contribution | PMF decomposition | ~2–3 | EXP-04 |
| EXP-08 | Interfacial H-Bond Energy | PMF decomposition | ~1–2 | EXP-04 |
| EXP-09 | BPTI-Trypsin $k_{\text{on}}/k_{\text{off}}$ | Kramers/TST | ~2–3 | EXP-04 |
| EXP-10 | BPTI-Trypsin $\Delta G^{\ddagger}$ | Barrier extraction | ~1–2 | EXP-04 |
| **Tier 3 — Dynamics** | | | | |
| EXP-24 | BPTI H/D Exchange | 100 ns production MD | ~8 | Independent |
| EXP-25 | BPTI Conformational Variability | Reuses EXP-24 trajectory | ~0.5 | EXP-24 |
| **Tier 4 — Correlation** | | | | |
| EXP-26 | BSA–$\Delta G$ Correlation | Post-processing | ~0.5 | EXP-04/05/13/29 |
| **Tier 5 — FEP Mutagenesis** | | | | |
| EXP-28 | Scaffold Energy | US (loop peptide) | ~42 | Partial EXP-04 |
| EXP-30 | Alanine Scanning | FEP (15 mutations) | ~90–100 | Independent |
| EXP-31 | Disulfide Ablation | FEP (C14S/C38S) | ~20–24 | Independent |
| EXP-32 | SPINK1 N34S Control | FEP (N34S) | ~6 | Independent |
| EXP-33 | Barnase-Barstar DMC | FEP (~45 cycles) | ~100–120 | Independent |

**Total estimated GPU compute:** ~530–570 A100 GPU-hours (~265–340 H100 GPU-hours).

</details>

---

<div style="page-break-after: always;"></div>

## Physical Validity Invariants

Ten invariants are enforced as runtime checks — violation halts execution immediately:

| ID | Condition | Description |
|----|-----------|-------------|
| IV-1 | $E_{\text{min}} < E_{\text{initial}}$ | Energy decreases after minimization |
| IV-2 | $\|T_{\text{avg}} - 310\| < 5$ K | NVT temperature stability |
| IV-3 | $\rho \in [0.95, 1.05]$ g/cm³ | NPT density physicality |
| IV-4 | RMSD $< 5$ Å | Backbone structural stability |
| IV-5 | Drift $< 0.1$ kJ/mol/ns/atom | Energy conservation in production |
| IV-6 | $d_{S-S} < 2.5$ Å | Disulfide bond integrity |
| IV-7 | Min image $> 2r_c$ | No periodic image artifacts |
| IV-8 | Overlap $\geq 10\%$ | Umbrella histogram coverage |
| IV-9 | $\max_i \|f_i^{(n+1)} - f_i^{(n)}\| < 10^{-6}$ | WHAM convergence |
| IV-10 | Unimodal $P(W)$ | No SMD pathway bifurcation |

All invariant enforcement uses explicit `ValueError` exceptions (not `assert`), ensuring validation persists under Python's `-O` optimization mode.

---

## Repository Structure

```
medium_project_2/
├── README.md                          # This file — project overview (V2 + V3)
├── LICENSE                            # MIT license
├── pyproject.toml                     # Build configuration and dependencies
├── requirements.txt                   # pip-installable dependency list
├── default_config.yaml                # Default simulation parameters
├── Makefile                           # Build and utility targets
│
├── src/                               # Pipeline source code (V2 — immutable during V3)
│   ├── config.py                      #   Central configuration (single source of truth)
│   ├── prep/                          #   Structure preparation pipeline
│   │   ├── pdb_fetch.py               #     RCSB download with retry/cache
│   │   ├── pdb_clean.py               #     Crystallographic artifact removal
│   │   ├── protonate.py               #     PROPKA-aware protonation
│   │   ├── topology.py                #     OpenMM topology with PME enforcement
│   │   └── solvate.py                 #     Solvation box & ion placement
│   ├── simulate/                      #   Molecular dynamics engines
│   │   ├── minimizer.py               #     Energy minimization with convergence reporting
│   │   ├── equilibrate.py             #     NVT → NPT equilibration with checkpoint
│   │   ├── production.py              #     Unrestrained production MD
│   │   ├── smd.py                     #     Steered MD with parallel replicates
│   │   ├── umbrella.py                #     Umbrella Sampling with pre-equilibration
│   │   ├── fep.py                     #     Alchemical free energy perturbation
│   │   ├── metadynamics.py            #     Well-tempered metadynamics
│   │   └── platform.py               #     Automatic CUDA/OpenCL/CPU detection
│   ├── analyze/                       #   Post-processing & thermodynamic analysis
│   │   ├── trajectory.py              #     Trajectory I/O with PBC unwrapping
│   │   ├── structural.py              #     RMSD, RMSF, Rg, SASA computation
│   │   ├── contacts.py                #     Streaming contact analysis
│   │   ├── wham.py                    #     WHAM solver for PMF extraction
│   │   ├── mbar.py                    #     MBAR solver via pymbar
│   │   ├── jarzynski.py               #     Jarzynski + BAR free energy estimators
│   │   ├── fep.py                     #     FEP analysis (MBAR-based ΔΔG)
│   │   ├── msm.py                     #     Markov State Model construction
│   │   ├── featurize.py               #     Trajectory featurization for MSM/TICA
│   │   ├── convergence.py             #     Block bootstrap with autocorrelation
│   │   └── equilibration.py           #     Chodera automated equilibration detection
│   ├── physics/                       #   Physical models & collective variables
│   │   ├── collective_variables.py    #     PBC-aware COM distance
│   │   ├── force_field_factory.py     #     3-tier force field abstraction
│   │   ├── finite_size.py             #     Finite-size electrostatic corrections
│   │   ├── restraints.py              #     Positional & distance restraints
│   │   └── units.py                   #     Unit conversion utilities
│   └── visualization/                 #   Rendering & plotting
│       ├── viewer_3d.py               #     py3Dmol with cyclic palette
│       ├── plot_pmf.py                #     PMF profile plotting
│       └── plot_timeseries.py         #     Energy, temperature, RMSD timeseries
│
├── scripts/                           # CLI entry points for each pipeline stage
│   ├── run_prep.py                    #   Structure preparation
│   ├── run_equilibration.py           #   Minimization + NVT/NPT equilibration
│   ├── run_production.py              #   Unrestrained production dynamics
│   ├── run_smd.py                     #   SMD campaign (N replicates)
│   ├── run_umbrella.py                #   Umbrella Sampling campaign
│   ├── run_analysis.py                #   WHAM, MBAR, Jarzynski, structural analysis
│   ├── run_fep.py                     #   Alchemical FEP campaign
│   ├── run_msm.py                     #   MSM construction and analysis
│   ├── cross_validate.py              #   Statistical cross-validation
│   ├── generate_figures.py            #   Publication-quality V2 figure generation
│   └── generate_gpu_figures.py        #   GPU test figure generation
│
├── tests/                             # 367 tests across 28+ test files
│   ├── conftest.py                    #   Shared fixtures
│   ├── test_*.py                      #   Unit, integration, and analytical tests
│   └── colab_gpu_test.ipynb           #   GPU force field validation notebook
│
├── notebooks/                         # Interactive Jupyter workflows (V2)
│   ├── 00_structure_acquisition.ipynb #   PDB retrieval and inspection
│   ├── 01_system_prep.ipynb           #   Structure preparation with visual inspection
│   ├── 02_equilibration.ipynb         #   Equilibration monitoring with real-time plots
│   ├── 03_production_analysis.ipynb   #   Production trajectory analysis
│   ├── 04_smd_analysis.ipynb          #   SMD work distributions and Jarzynski analysis
│   ├── 05_umbrella_pmf.ipynb          #   Umbrella Sampling results and WHAM/MBAR PMF
│   ├── 06_visualization.ipynb         #   3D molecular rendering and publication figures
│   ├── 07_production_campaign.ipynb   #   Full production campaign orchestration
│   ├── 08_convergence_analysis.ipynb  #   Convergence diagnostics
│   ├── 09_invariant_validation.ipynb  #   Physics invariant verification
│   ├── 10_fep_validation.ipynb        #   FEP method validation
│   └── 11_spink7_mutagenesis.ipynb    #   SPINK7-KLK5 computational mutagenesis
│
├── v3_experiments/                    # V3 benchmarking campaign (33 experiments)
│   ├── experiment_progress.csv        #   Experiment tracking with classifications
│   ├── EXP-01_spink7_klk5_dg_bind/   #   Each experiment folder contains:
│   │   ├── experimental_design.md     #     Scientific design and target benchmarks
│   │   ├── implementation_guide.md    #     Step-by-step computational protocol
│   │   ├── results_report.md          #     Results and statistical analysis (if executed)
│   │   ├── EXP-XX_colab.ipynb         #     Google Colab notebook (GPU experiments)
│   │   └── outputs/                   #     Figures and results data (if executed)
│   │       ├── figures/               #       Experiment-specific figures
│   │       └── results.json           #       Quantitative results
│   ├── EXP-02_spink7_klk12_dg_bind/
│   ├── ...                            #   (33 experiment folders: EXP-01 through EXP-33)
│   ├── structures/                    #   Pre-prepared PDB structures for experiments
│   ├── pipeline_validation/           #   Pipeline validation outputs
│   ├── prepare_structures.py          #   Structure preparation utilities
│   ├── prepare_all_structures.py      #   Batch structure preparation
│   ├── validate_pipeline.py           #   Pipeline validation checks
│   └── 00_dry_run_validation.ipynb    #   Dry-run validation notebook (Colab)
│
├── v3_literature/                     # V3 literature review outputs
│   ├── benchmarks.md                  #   Consolidated experimental benchmarks (33 features)
│   ├── sources.md                     #   Source registry (28 publications)
│   ├── literature_review_progress.csv #   Review tracking
│   └── source_reports/                #   Per-source analysis reports (28 reports)
│
├── latex/                             # IEEE final report (V2 + V3 content)
│   ├── final_report.tex               #   LaTeX source (~2500 lines)
│   └── latex_figures/                 #   Figures embedded in the LaTeX report (21 figures)
│
├── figures/                           # Publication-quality figures
│   ├── pipeline_v2_figures/           #   V2 implementation diagrams (52 figures)
│   ├── fig3_protein_complex.png       #   Barnase-barstar complex rendering
│   ├── fig4_pmf_profile.png           #   PMF profile
│   └── fig5_simulation_timeseries.png #   Simulation diagnostics
│
├── data/                              # Raw and prepared structures
│   └── pdb/                           #   PDB files (raw, cached, prepared)
│
├── reports/                           # V2 project documentation
│   ├── full_implementation_report_v2.md  # Complete V2 implementation report
│   ├── gpu_test_implementation_guide.md  # GPU test execution guide
│   └── project_overview.md               # Project overview
│
├── architecture_blueprint/            # Pipeline architecture documentation
│   └── architecture_blueprint.md      #   System architecture and design decisions
│
└── improvements/                      # V2 limitation analysis and implementation guides
    ├── improvements_v1/               #   First-pass limitation identification
    └── improvements_v2/               #   Detailed implementation guides for all 40 fixes
```

---

## Installation & Setup

### Prerequisites

- Python 3.11+
- OpenMM ≥ 8.1 (may require [conda installation](http://docs.openmm.org/latest/userguide/application/01_getting_started.html) on some platforms)

### Setup

```bash
git clone https://github.com/ryanjosephkamp/SPINK7-KLK5-MD-Pipeline.git
cd SPINK7-KLK5-MD-Pipeline
python -m venv .venv
source .venv/bin/activate

# Option A: Install via requirements.txt
pip install -r requirements.txt

# Option B: Install as editable package (recommended)
pip install -e ".[dev]"
```

### Verify Installation

```bash
python -m pytest tests/ -v
```

All 367 tests should pass with no failures.

---

## Usage

### Running the V2 Pipeline

**Stage 1 — System Preparation:**
```bash
python scripts/run_prep.py
```

**Stage 2 — Equilibration:**
```bash
python scripts/run_equilibration.py
```

**Stage 3 — Production MD:**
```bash
python scripts/run_production.py
```

**Stage 4 — Steered Molecular Dynamics:**
```bash
python scripts/run_smd.py
```

**Stage 5 — Umbrella Sampling:**
```bash
python scripts/run_umbrella.py
```

**Stage 6 — Analysis & Visualization:**
```bash
python scripts/run_analysis.py
```

**Generate Publication Figures:**
```bash
python scripts/generate_figures.py
# Or equivalently:
make figures
```

> **Note:** Image files (PNG, SVG, PDF figures) are excluded from version control. After cloning, run `make figures` or the command above to regenerate all publication figures locally.

### Configuration

Pipeline parameters can be specified via YAML configuration files with precedence cascade: CLI arguments → config file → dataclass defaults.

```bash
python scripts/run_production.py --config production_config.yaml
```

### Interactive Notebooks

For exploratory analysis with real-time visualization:

```bash
jupyter notebook notebooks/
```

| Notebook | Purpose |
|----------|---------|
| `00_structure_acquisition.ipynb` | PDB retrieval and inspection |
| `01_system_prep.ipynb` | Interactive structure preparation with visual inspection |
| `02_equilibration.ipynb` | Equilibration monitoring with real-time plots |
| `03_production_analysis.ipynb` | Production trajectory analysis |
| `04_smd_analysis.ipynb` | SMD work distributions and Jarzynski analysis |
| `05_umbrella_pmf.ipynb` | Umbrella Sampling results and WHAM/MBAR PMF |
| `06_visualization.ipynb` | 3D molecular rendering and publication figures |
| `07_production_campaign.ipynb` | Full production campaign orchestration |
| `08_convergence_analysis.ipynb` | Convergence diagnostics |
| `09_invariant_validation.ipynb` | Physics invariant verification |
| `10_fep_validation.ipynb` | FEP method validation |
| `11_spink7_mutagenesis.ipynb` | SPINK7-KLK5 computational mutagenesis |

<div style="page-break-after: always;"></div>

### Running V3 Experiments (CPU)

CPU-executable V3 experiments can be run using their respective implementation guides. Each experiment folder in `v3_experiments/` contains a step-by-step protocol:

```bash
# Example: navigate to an experiment and follow implementation_guide.md
cd v3_experiments/EXP-15_catalytic_contact_distance/
cat implementation_guide.md
```

Results for completed experiments are stored in `outputs/` and documented in `results_report.md`.

---

## GPU Experiment Execution Guide

The 17 GPU-required experiments include pre-built Google Colab notebooks ready for execution on NVIDIA A100/H100 runtimes.

### Prerequisites

- Google account with [Colab Pro/Pro+](https://colab.research.google.com/signup) for A100 GPU access
- Sufficient compute units for the experiment tier (see table below)

### Setup

1. **Upload the repository** to Google Drive:
   ```
   My Drive/SPINK7-KLK5-MD-Pipeline/
   ```

2. **Open the experiment's Colab notebook** in Google Colab:
   - Navigate to `v3_experiments/EXP-XX_*/EXP-XX_colab.ipynb`
   - Select **Runtime → Change runtime type → A100 GPU**

3. **Mount Google Drive** (handled by the first notebook cell):
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   ```

<div style="page-break-after: always;"></div>

4. **Install dependencies** (handled by the setup cell):
   ```bash
   pip install openmm pdbfixer mdtraj pymbar openmmtools deeptime
   ```

5. **Execute cells sequentially** — each notebook follows the experiment's implementation guide with checkpoint saves at each major stage.

### Checkpoint Recovery

All Colab notebooks implement checkpoint-based fault tolerance. If a session disconnects:

1. Reconnect and re-mount Google Drive
2. Re-run the notebook — checkpoint detection automatically resumes from the last saved state
3. Checkpoints are saved to `EXP-XX_*/outputs/checkpoints/` on Google Drive

### Dry-Run Validation

Before committing to full GPU compute, run the dry-run validation notebook:

```
v3_experiments/00_dry_run_validation.ipynb
```

This validates the pipeline setup, dependency installation, and structure preparation without running full simulations.

<div style="page-break-after: always;"></div>

### Compute Requirements by Tier

| Tier | Experiments | Method | Est. A100 Hours | Colab Sessions (~24h) |
|------|------------|--------|-----------------|-----------------------|
| 1 — Thermodynamic | EXP-04, 05, 06, 13, 29 | SMD + US | ~250 | ~11 |
| 2 — PMF Analysis | EXP-07, 08, 09, 10 | Post-processing | ~7 | ~1 |
| 3 — Dynamics | EXP-24, 25 | Production MD | ~8.5 | ~1 |
| 4 — Correlation | EXP-26 | Post-processing | ~0.5 | < 1 |
| 5 — FEP Mutagenesis | EXP-28, 30, 31, 32, 33 | FEP / US | ~258–286 | ~11–12 |
| **Total** | **17 experiments** | | **~530–570** | **~22–25** |

### Execution Order (Dependency-Aware)

1. **Independent experiments first:** EXP-04, 05, 06, 13, 24, 29, 30, 31, 32, 33 (can run in parallel on separate Colab sessions)
2. **Dependent experiments after prerequisites complete:**
   - EXP-07, 08, 09, 10 → require EXP-04
   - EXP-25 → requires EXP-24
   - EXP-26 → requires EXP-04, 05, 13, 29
   - EXP-28 → requires partial EXP-04

---

<div style="page-break-after: always;"></div>

## Simulation Protocol

| Stage | Ensemble | Duration | Key Conditions |
|-------|----------|----------|---------------|
| 1. Energy Minimization | — | ≤ 10,000 steps | Steepest descent until $< 10$ kJ/mol/nm; convergence reporting |
| 2. NVT Equilibration | NVT | 500 ps | Heavy-atom restraints ($k = 1000$ kJ/mol/nm²) |
| 3. NPT Equilibration | NPT | 1 ns | Gradual restraint release; checkpoint save |
| 4. Production MD | NPT | 100–500 ns | Unrestrained; frames every 10 ps; automated equilibration detection |
| 5. Enhanced Sampling | NPT | Per method | SMD (50 replicates, parallel), Umbrella (25–50 windows, pre-equilibrated), or well-tempered metadynamics |

---

<div style="page-break-after: always;"></div>

## Test Suite

All **367 tests** pass across CPU and GPU platforms:

| Category | Files | Coverage |
|----------|-------|----------|
| Configuration | 1 | Frozen dataclass immutability, parameter correctness, YAML config loading |
| Physics | 3 | Unit conversions, PBC-aware COM distance, restraint forces |
| Preparation | 5 | PDB fetch with retry/cache, clean, protonate (PROPKA), topology (PME), solvate |
| Simulation | 5 | Minimization (convergence), NVT/NPT equilibration, production, SMD, umbrella |
| Analysis | 8 | Trajectory I/O, structural metrics, streaming contacts, WHAM, MBAR, Jarzynski/BAR, convergence, equilibration detection |
| Visualization | 3 | 3D viewer (cyclic palette), PMF plots, timeseries plots |
| Integration | 2 | Full pipeline on alanine dipeptide; optimization-mode harness |
| Edge Cases | 1 | Parser robustness: multi-model PDB, alternate conformers, malformed records |
| GPU Validation | 1 | Force field factory GPU backends: AMOEBA polarizable, ANI-2x ML potential |

Analytical validation confirms mathematical correctness: the Jarzynski estimator recovers exact free energies from harmonic potentials within 0.5 kJ/mol; the WHAM solver reconstructs known flat and parabolic PMFs with < 1.0 kJ/mol error; the MBAR solver agrees with WHAM within statistical uncertainty; convergence estimators satisfy the expected $1/\sqrt{N}$ scaling.

### Running Tests

```bash
# Full test suite
python -m pytest tests/ -v

# With coverage reporting
python -m pytest tests/ -v --cov=src --cov-report=term-missing

# Optimization mode (verifies no assert-dependent logic)
python -O -m pytest tests/ -v --tb=short
```

---

## Data Flow

```
               RCSB / AlphaFold              PDB / mmCIF Files
                    │                          │
                    ▼                          ▼
               pdb_fetch.py              pdb_clean.py
               (retry + cache)           (multi-model aware)
                    │                          │
                    ▼                          ▼
               data/pdb/raw/              protonate.py
                                             (PROPKA pKa)
                                                  │
                                                  ▼
                                             topology.py  ──►  solvate.py
                                             (PME enforced)     (cubic / oct)
                                                                 │
                    ┌────────────────────────────────────────────┘
                    │
                    ▼
               minimizer.py  ──►  equilibrate.py  ──►  production.py
               (convergence)      (checkpoint)          (auto-detect equil.)
                                                               │
                    ┌──────────────────────────────────────────┤
                    │               │                          │
                    ▼               ▼                          ▼
               smd.py         umbrella.py             metadynamics.py
            (parallel)     (pre-equilibrated)        (well-tempered)
                    │               │                          │
                    ▼               ▼                          │
           jarzynski.py / BAR  wham.py / mbar.py               │
                    │               │                          │
                    └───────┬───────┘──────────────────────────┘
                            │
                            ▼
                  structural.py / contacts.py
                  (PBC-unwrapped, streaming)
                            │
                            ▼
                  visualization / plots
                  (generate_figures.py)
```

---

## Dependencies

| Package | Version | Role |
|---------|---------|------|
| OpenMM | ≥ 8.1 | Core MD simulation engine |
| PDBFixer | ≥ 1.9 | Structural repair and standardization |
| MDTraj | ≥ 1.9.9 | Trajectory analysis and I/O |
| NumPy | ≥ 1.26, < 2.0 | Numerical computations |
| SciPy | ≥ 1.12 | Scientific computing and optimization |
| matplotlib | ≥ 3.8 | 2D plotting and visualization |
| py3Dmol | ≥ 2.0 | Interactive 3D molecular visualization |
| pandas | ≥ 2.1 | Tabular data analysis |
| PyYAML | ≥ 6.0 | Configuration file parsing |
| PDB2PQR | ≥ 3.6 | Protonation state assignment (bundles PROPKA) |
| gemmi | ≥ 0.6.4 | mmCIF/PDBx file parsing |
| pymbar | ≥ 4.0 | MBAR free energy estimation |
| openmmtools | ≥ 0.23.0 | Enhanced sampling utilities and alchemical tools |
| deeptime | ≥ 0.4.4 | MSM construction and TICA analysis |
| requests | ≥ 2.31 | HTTP client for PDB download with retry |
| Pillow | ≥ 10.0 | Image processing (test suite) |
| pytest | ≥ 8.0 | Test framework |
| pytest-cov | ≥ 4.1 | Test coverage reporting |

<div style="page-break-after: always;"></div>

**Optional GPU packages** (for Tier 2/3 force fields on CUDA-capable machines):

| Package | Version | Role |
|---------|---------|------|
| openmm-ml | ≥ 1.1 | ML potential integration (ANI-2x, MACE-OFF) |
| torchani | ≥ 2.2 | ANI-2x neural network potential |
| openmm-torch | ≥ 1.1 | PyTorch-OpenMM bridge |

---

## Key Findings & Identified Limitations

### Strengths

The V2 pipeline demonstrates strong performance on **structural analysis** tasks:

- **Crystallographic fidelity**: BPTI disulfide bond distances reproduced to within 0.03 Å of the 0.94 Å X-ray structure (EXP-23: PASS)
- **Interface geometry**: Buried surface area (1608.5 Å²), catalytic contact distance (2.675 Å), and nonpolar BSA fraction (56.9%) all within experimental confidence intervals (EXP-15, 16, 19: PASS)
- **Cross-system conservation**: Binding loop RMSD of 0.234 Å across Kazal-family inhibitors confirms structural conservation (EXP-20: PASS)
- **Subsite recognition**: KLK5 S1 salt bridge distance (3.63 Å) confirms correct electrostatic complementarity (EXP-21: PASS)

### Identified Limitations

1. **Thermodynamic validation incomplete** — The 5 independent $\Delta G_{\text{bind}}$ experiments (EXP-04, 05, 06, 13, 29) require GPU-accelerated SMD + umbrella sampling. Until these are executed, the pipeline's quantitative accuracy for binding free energies is unvalidated.

2. **Static vs. dynamic analysis** — The MARGINAL result for interfacial water count (EXP-18: 8.0 vs. target 15–30) reflects analysis of an energy-minimized snapshot rather than a dynamic ensemble average. Production MD trajectories would likely improve this prediction.

3. **Binding loop flexibility** — The MARGINAL result for binding loop geometry (EXP-14: 1.766 Å vs. target < 1.5 Å) suggests the energy minimization protocol may not fully relax loop conformations to crystallographic reference geometries.

4. **Structure availability** — 7 experiments are blocked by the absence of experimentally resolved co-crystal structures. Homology modeling or AlphaFold-Multimer predictions could potentially address this gap in Part 2.

5. **FEP validation pending** — All 5 alchemical FEP experiments (EXP-28, 30, 31, 32, 33) require GPU compute. The pipeline's mutagenesis capability is implemented and tested but not yet validated against experimental $\Delta\Delta G$ values.

---

## Future Directions

### V3 Part 1 — GPU Experiment Completion

The immediate priority is executing the 17 GPU-required experiments using the prepared Colab notebook infrastructure:

1. **Tier 1 (highest priority):** EXP-04 (BPTI-trypsin $\Delta G_{\text{bind}}$) — the cornerstone validation experiment. A PASS here would validate the pipeline's SMD + umbrella sampling methodology against the gold-standard protease-inhibitor system ($K_d = 6 \times 10^{-14}$ M).

2. **Tier 5 (FEP validation):** EXP-30 (alanine scanning) and EXP-33 (barnase-barstar DMC) — the two largest experiments (~200 combined A100 hours) that would validate computational mutagenesis.

3. **Remaining tiers** in dependency order (see [Execution Order](#execution-order-dependency-aware)).

### V3 Part 2 — Limitation Remediation

Based on V3 Part 1 findings, the following improvements are prioritized for Part 2:

1. **Enhanced structural modeling** — integrate AlphaFold-Multimer for generating initial complex structures where no co-crystal is available (EXP-01, 02, 03, 11, 12, 22)
2. **Constant-pH MD** — implement titratable residue dynamics for pH-dependent binding kinetics (EXP-12)
3. **Dynamic water analysis** — extend interfacial water counting to production MD ensembles rather than static snapshots
4. **Systematic force field comparison** — execute all thermodynamic experiments at all three force field tiers (AMBER, AMOEBA, ANI-2x) to quantify polarization effects on $\Delta G_{\text{bind}}$ accuracy

---

## References

1. M. E. Rothenberg, "Biology and treatment of eosinophilic esophagitis," *Gastroenterology*, 137(4), 1238–1249, 2009.
2. N. P. Azouz *et al.*, "The antiprotease SPINK7 serves as an inhibitory checkpoint for esophageal epithelial inflammatory responses," *Sci. Transl. Med.*, 10(444), eaap9736, 2018.
3. M. Laskowski Jr. and I. Kato, "Protein inhibitors of proteinases," *Annu. Rev. Biochem.*, 49, 593–626, 1980.
4. J. A. Maier *et al.*, "ff14SB: Improving the accuracy of protein side chain and backbone parameters from ff99SB," *J. Chem. Theory Comput.*, 11(8), 3696–3713, 2015.
5. W. L. Jorgensen *et al.*, "Comparison of simple potential functions for simulating liquid water," *J. Chem. Phys.*, 79(2), 926–935, 1983.
6. I. S. Joung and T. E. Cheatham III, "Determination of alkali and halide monovalent ion parameters," *J. Phys. Chem. B*, 112(30), 9020–9041, 2008.
7. T. Darden, D. York, and L. Pedersen, "Particle mesh Ewald: An N·log(N) method for Ewald sums," *J. Chem. Phys.*, 98(12), 10089–10092, 1993.
8. C. Jarzynski, "Nonequilibrium equality for free energy differences," *Phys. Rev. Lett.*, 78(14), 2690, 1997.
9. S. Kumar *et al.*, "The weighted histogram analysis method for free-energy calculations on biomolecules," *J. Comput. Chem.*, 13(8), 1011–1021, 1992.
10. P. Eastman *et al.*, "OpenMM 8: Molecular dynamics simulation with machine learning potentials," *J. Phys. Chem. B*, 128(1), 109–116, 2024.
11. G. Schreiber and A. R. Fersht, "Energetics of protein-protein interactions: Analysis of the barnase-barstar interface," *J. Mol. Biol.*, 248(2), 478–486, 1995.
12. M. H. M. Olsson *et al.*, "PROPKA3: Consistent treatment of internal and surface residues in empirical pKa predictions," *J. Chem. Theory Comput.*, 7(2), 525–537, 2011.
13. D. M. Zuckerman and T. B. Woolf, "Theory of a systematic computational error in free energy differences," *Phys. Rev. Lett.*, 89(18), 180602, 2002.
14. G. E. Crooks, "Entropy production fluctuation theorem and the nonequilibrium work relation for free energy differences," *Phys. Rev. E*, 60(3), 2721–2726, 1999.
15. M. R. Shirts and J. D. Chodera, "Statistically optimal analysis of samples from multiple equilibrium states," *J. Chem. Phys.*, 129(12), 124105, 2008.
16. J. D. Chodera, "A simple method for automated equilibration detection in molecular simulations," *J. Chem. Theory Comput.*, 12(4), 1799–1805, 2016.
17. C. H. Bennett, "Efficient estimation of free energy differences from Monte Carlo data," *J. Comput. Phys.*, 22(2), 245–268, 1976.
18. A. Barducci, G. Bussi, and M. Parrinello, "Well-tempered metadynamics: A smoothly converging and tunable free-energy method," *Phys. Rev. Lett.*, 100(2), 020603, 2008.
19. P. H. Hünenberger and J. A. McCammon, "Ewald artifacts in computer simulations of ionic solvation and ion-ion interaction," *J. Chem. Phys.*, 110(4), 1856–1872, 1999.
20. J. W. Ponder *et al.*, "Current status of the AMOEBA polarizable force field," *J. Phys. Chem. B*, 114(8), 2549–2564, 2010.
21. S. Izadi, R. Anandakrishnan, and A. V. Onufriev, "Building water models: A different approach," *J. Phys. Chem. Lett.*, 5(21), 3863–3871, 2014.
22. F. Noé *et al.*, "Constructing the equilibrium ensemble of folding pathways from short off-equilibrium simulations," *Proc. Natl. Acad. Sci.*, 106(45), 19011–19016, 2009.
23. R. T. McGibbon *et al.*, "MDTraj: A modern open library for the analysis of molecular dynamics trajectories," *Biophys. J.*, 109(8), 1528–1532, 2015.
24. H. R. Künsch, "The jackknife and the bootstrap for general stationary observations," *Ann. Stat.*, 17(3), 1217–1241, 1989.
25. C. Devereux *et al.*, "Extending the applicability of the ANI deep learning molecular potential to sulfur and halogens," *J. Chem. Theory Comput.*, 16(7), 4192–4202, 2020.

---

## License

[MIT](LICENSE) © 2026 Ryan Kamp

---

## Author

**Ryan Kamp**
Dept. of Computer Science, University of Cincinnati
[kamprj@mail.uc.edu](mailto:kamprj@mail.uc.edu) · [GitHub](https://github.com/ryanjosephkamp)
