# SPINK7-KLK5 MD Pipeline: Binding Free Energy via Enhanced Sampling

[![OpenMM](https://img.shields.io/badge/OpenMM-%E2%89%A58.1-blue.svg)](https://openmm.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Tests: 259 passed](https://img.shields.io/badge/tests-259%20passed-brightgreen.svg)]()
[![Force Field](https://img.shields.io/badge/force%20field-AMBER%20ff14SB-orange.svg)]()

An end-to-end molecular dynamics simulation pipeline for computing the binding free energy of the **SPINK7-KLK5** protease-antiprotease complex — a protein-protein interaction central to the pathogenesis of **Eosinophilic Esophagitis (EoE)**. The pipeline automates the full computational biophysics workflow: PDB retrieval, structure cleaning, protonation at physiological pH, AMBER ff14SB topology construction, explicit TIP3P solvation, energy minimization, NVT/NPT equilibration, and unrestrained production dynamics. Two complementary enhanced sampling strategies — **Steered Molecular Dynamics (SMD)** with the **Jarzynski equality** and **Umbrella Sampling** with **WHAM** reconstruction of the Potential of Mean Force — provide rigorous, cross-validated binding free energy estimates with quantified statistical uncertainties.

Ten physical validity invariants are enforced as runtime assertions throughout the pipeline. A comprehensive test suite of **259 unit, integration, and analytical tests** validates every module, all passing on a CPU-only platform.

---

<div style="page-break-after: always;"></div>

## Biological Motivation

Eosinophilic Esophagitis (EoE) is a chronic inflammatory disease of the esophagus driven by IL-13-mediated transcriptional silencing of *SPINK7* (Serine Peptidase Inhibitor, Kazal Type 7). Under homeostatic conditions, SPINK7 stoichiometrically inhibits KLK5 (Kallikrein-Related Peptidase 5), a trypsin-like serine protease. SPINK7 deficiency unleashes KLK5 proteolytic activity, which degrades Desmoglein-1 (DSG1) and compromises the epithelial barrier, permitting allergen penetration.

The SPINK7-KLK5 interaction follows the canonical **Laskowski mechanism** for Kazal-type inhibitor–serine protease binding: the reactive site loop (RSL) of SPINK7 inserts into the KLK5 active-site cleft in a substrate-like orientation, and the conformational rigidity imposed by three disulfide bonds renders the acyl-enzyme intermediate thermodynamically trapped with an extremely slow $k_{\text{off}}$. The binding interface buries approximately 800–1200 Å² of solvent-accessible surface area, stabilized by backbone hydrogen bonds, electrostatic complementarity (P1-Arg/Lys ↔ Asp189), and hydrophobic contacts at the P2/P2' sub-sites.

Understanding this interaction at the atomic level through molecular dynamics simulation is essential for developing targeted therapeutic strategies to restore epithelial barrier function in EoE patients.

---

## Mathematical Framework

### Force Field Potential Energy

The total potential energy is decomposed into bonded and nonbonded contributions:

$$V_{\text{total}}(\mathbf{r}) = V_{\text{bonded}}(\mathbf{r}) + V_{\text{nonbonded}}(\mathbf{r})$$

**Bonded terms** (covalent interactions within the molecular topology):

$$V_{\text{bonded}} = \sum_{\text{bonds}} K_b (b - b_0)^2 + \sum_{\text{angles}} K_\theta (\theta - \theta_0)^2 + \sum_{\text{dihedrals}} \frac{V_n}{2} [1 + \cos(n\phi - \gamma)] + \sum_{\text{impropers}} K_\omega (\omega - \omega_0)^2$$

where $K_b$, $K_\theta$, $V_n$, and $K_\omega$ are force constants; $b_0$, $\theta_0$, $\gamma$, and $\omega_0$ are equilibrium values; and $n$ is the dihedral periodicity.

**Nonbonded terms** (van der Waals + electrostatics):

$$V_{\text{nonbonded}} = \sum_{i < j} \left[ 4\epsilon_{ij} \left( \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6} \right) + \frac{q_i q_j}{4\pi\epsilon_0 r_{ij}} \right]$$

The **AMBER ff14SB** force field provides optimized backbone torsion parameters. Water is modeled using **TIP3P**, and ions follow the **Joung-Cheatham** monovalent parameters.

### Long-Range Electrostatics (Particle Mesh Ewald)

$$E_{\text{elec}} = E_{\text{direct}} + E_{\text{reciprocal}} + E_{\text{self-correction}}$$

PME reduces the computational cost from $O(N^2)$ to $O(N \log N)$ using B-spline interpolation (order 5) with grid spacing $\leq 1.0$ Å and direct-space cutoff $r_c = 10$ Å.

### Langevin Equation of Motion

$$m_i \ddot{\mathbf{r}}_i = -\nabla_i V(\mathbf{r}) - \gamma m_i \dot{\mathbf{r}}_i + \sqrt{2 \gamma m_i k_B T} \, \boldsymbol{\eta}_i(t)$$

where $\gamma = 1.0 \text{ ps}^{-1}$, $T = 310$ K, and $\boldsymbol{\eta}_i(t)$ is Gaussian white noise. Integration uses the **Langevin middle integrator** with $\Delta t = 2$ fs, enabled by **SHAKE** constraints on hydrogen bonds.

### Steered Molecular Dynamics

SMD applies a time-dependent harmonic bias to the center-of-mass (COM) distance $\xi$:

$$U_{\text{SMD}}(\xi, t) = \frac{k}{2} \left[ \xi(t) - \xi_0 - v \cdot t \right]^2$$

where $k = 1000$ kJ/mol/nm², $v = 0.001$ nm/ps, and the COM distance is:

$$\xi(\mathbf{r}) = \left\| \mathbf{R}_{\text{COM}}^{\text{SPINK7}} - \mathbf{R}_{\text{COM}}^{\text{KLK5}} \right\|_2$$

The **Jarzynski equality** connects non-equilibrium work to equilibrium free energy:

$$\Delta G = -k_B T \ln \left[ \frac{1}{N_{\text{traj}}} \sum_{j=1}^{N_{\text{traj}}} e^{-\beta W_j} \right]$$

with the second-order cumulant expansion for near-Gaussian work distributions:

$$\Delta G \approx \langle W \rangle - \frac{\beta}{2} \sigma_W^2$$

<div style="page-break-after: always;"></div>

### Umbrella Sampling & WHAM

Each of $M$ discrete windows is biased by a harmonic potential:

$$U_i^{\text{bias}}(\xi) = \frac{k_i}{2} \left( \xi - \xi_i^{\text{ref}} \right)^2$$

The **Weighted Histogram Analysis Method** recovers the unbiased PMF via self-consistent equations:

$$P^{\text{unbiased}}(\xi) = \frac{\sum_{i=1}^{M} n_i \, h_i(\xi)}{\sum_{i=1}^{M} n_i \, \exp\left[ \beta \left( f_i - U_i^{\text{bias}}(\xi) \right) \right]}$$

$$e^{-\beta f_i} = \int P^{\text{unbiased}}(\xi) \, \exp\left[ -\beta \, U_i^{\text{bias}}(\xi) \right] d\xi$$

iterated until convergence: $\max_i |f_i^{(n+1)} - f_i^{(n)}| < 10^{-6}$ kJ/mol. The PMF is then:

$$G(\xi) = -k_B T \ln P^{\text{unbiased}}(\xi) + C$$

### Binding Free Energy Extraction

$$\Delta G_{\text{bind}}^{\circ} = -k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{site}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right] + k_B T \ln \left[ \frac{C^{\circ}}{4\pi} \int_{\text{bulk}} e^{-\beta G(\xi)} \xi^2 \, d\xi \right]$$

where $C^{\circ} = 1/1660$ Å$^{-3}$ corresponds to the standard concentration of 1 M.

---

## Results

### Validation System: Barnase-Barstar Complex (PDB: 1BRS)

The pipeline was validated on the barnase-barstar complex, which has an experimentally measured binding free energy of $\Delta G_{\text{bind}} \approx -19$ kcal/mol. This system provides a stringent benchmark with a large binding interface (~1600 Å² buried SASA) dominated by electrostatic complementarity.

<p align="center">
  <img src="figures/fig3_protein_complex.png" alt="Barnase-Barstar Complex" width="600"/>
</p>

**Figure 1.** Three-dimensional Cα backbone rendering of the barnase-barstar complex (PDB: 1BRS). Chain A (barnase) is shown in blue, Chain D (barstar) in red. Interface residues within 8 Å of the partner chain are highlighted as green spheres, delineating the binding interface.

### PMF Profile and Binding Free Energy

<p align="center">
  <img src="figures/fig4_pmf_profile.png" alt="Potential of Mean Force" width="700"/>
</p>

**Figure 2.** Potential of Mean Force (PMF) profile as a function of COM distance $\xi$. The deep minimum near $\xi \approx 2.0$ nm corresponds to the bound state; the plateau at large distances represents the dissociated reference state. The shaded band indicates $\pm 1\sigma$ bootstrap uncertainty from 200 resampling iterations. The orange marker annotates the binding free energy $\Delta G_{\text{bind}}$ at the PMF minimum.

### Simulation Quality Diagnostics

<p align="center">
  <img src="figures/fig5_simulation_timeseries.png" alt="Simulation Timeseries" width="600"/>
</p>

**Figure 3.** Simulation timeseries diagnostics. **(a)** Potential energy (navy) and kinetic energy (orange) demonstrating energy conservation. **(b)** Temperature stability around 310 K, confirming Langevin thermostat function (invariant IV-2). **(c)** Backbone RMSD evolution showing equilibration plateau at ~0.15 nm, indicating structural stability (invariant IV-4).

---

<div style="page-break-after: always;"></div>

## Physical Validity Invariants

Ten invariants are enforced as runtime assertions — violation halts execution immediately:

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

---

## Repository Structure

```
medium_project_2/
├── src/
│   ├── config.py                      # Central configuration (single source of truth)
│   ├── prep/                          # Structure preparation pipeline
│   │   ├── pdb_fetch.py               #   RCSB download & AlphaFold retrieval
│   │   ├── pdb_clean.py               #   Crystallographic artifact removal
│   │   ├── protonate.py               #   Protonation state assignment (pH 7.4)
│   │   ├── topology.py                #   OpenMM Topology & System construction
│   │   └── solvate.py                 #   Solvation box & ion placement
│   ├── simulate/                      # Molecular dynamics engines
│   │   ├── minimizer.py               #   Energy minimization (steepest descent)
│   │   ├── equilibrate.py             #   NVT → NPT equilibration
│   │   ├── production.py              #   Unrestrained production MD
│   │   ├── smd.py                     #   Steered Molecular Dynamics engine
│   │   └── umbrella.py                #   Umbrella Sampling window manager
│   ├── analyze/                       # Post-processing & thermodynamic analysis
│   │   ├── trajectory.py              #   Trajectory I/O with streaming support
│   │   ├── structural.py              #   RMSD, RMSF, Rg, SASA computation
│   │   ├── contacts.py                #   Interface contact maps & H-bond analysis
│   │   ├── wham.py                    #   WHAM solver for PMF extraction
│   │   ├── jarzynski.py               #   Jarzynski free energy estimator
│   │   └── convergence.py             #   Block averaging & bootstrap resampling
│   ├── physics/                       # Physical models & collective variables
│   │   ├── collective_variables.py    #   COM distance calculations
│   │   ├── restraints.py              #   Positional & distance restraints
│   │   └── units.py                   #   Unit conversion utilities
│   └── visualization/                 # Rendering & plotting
│       ├── viewer_3d.py               #   py3Dmol interactive 3D widgets
│       ├── plot_pmf.py                #   PMF profile plotting
│       └── plot_timeseries.py         #   Energy, temperature, RMSD timeseries
├── scripts/                           # CLI entry points for each pipeline stage
│   ├── run_prep.py                    #   Structure preparation
│   ├── run_equilibration.py           #   Minimization + NVT/NPT equilibration
│   ├── run_production.py              #   Unrestrained production dynamics
│   ├── run_smd.py                     #   SMD campaign (N replicates)
│   ├── run_umbrella.py                #   Umbrella Sampling campaign
│   ├── run_analysis.py                #   WHAM, Jarzynski, structural analysis
│   └── generate_figures.py            #   Publication-quality figure generation
├── notebooks/                         # Interactive Jupyter workflows
│   ├── 01_system_prep.ipynb
│   ├── 02_equilibration.ipynb
│   ├── 03_production_analysis.ipynb
│   ├── 04_smd_analysis.ipynb
│   ├── 05_umbrella_pmf.ipynb
│   └── 06_visualization.ipynb
├── tests/                             # 259 tests across 24 test files
├── latex/                             # IEEE final report (LaTeX source)
├── figures/                           # Publication-quality figures
├── data/                              # Raw and prepared structures
└── reports/                           # Project documentation
```

---

## Simulation Protocol

| Stage | Ensemble | Duration | Key Conditions |
|-------|----------|----------|---------------|
| 1. Energy Minimization | — | ≤ 10,000 steps | Steepest descent until $< 10$ kJ/mol/nm |
| 2. NVT Equilibration | NVT | 500 ps | Heavy-atom restraints ($k = 1000$ kJ/mol/nm²) |
| 3. NPT Equilibration | NPT | 1 ns | Gradual restraint release |
| 4. Production MD | NPT | 100–500 ns | Unrestrained; frames every 10 ps |
| 5. Enhanced Sampling | NPT | Per method | SMD (50 replicates) or Umbrella Sampling (25–50 windows) |

---

## Test Suite

All **259 tests** pass on a CPU-only platform:

| Category | Files | Coverage |
|----------|-------|----------|
| Configuration | 1 | Frozen dataclass immutability, parameter correctness |
| Physics | 3 | Unit conversions, COM distance, restraint forces |
| Preparation | 5 | PDB fetch, clean, protonate, topology, solvate |
| Simulation | 5 | Minimization, NVT/NPT equilibration, production, SMD, umbrella |
| Analysis | 6 | Trajectory I/O, structural metrics, contacts, WHAM, Jarzynski, convergence |
| Visualization | 3 | 3D viewer, PMF plots, timeseries plots |
| Integration | 1 | Full pipeline on alanine dipeptide |

<div style="page-break-after: always;"></div>

Analytical validation confirms mathematical correctness: the Jarzynski estimator recovers exact free energies from harmonic potentials within 0.5 kJ/mol; the WHAM solver reconstructs known flat and parabolic PMFs with < 1.0 kJ/mol error; convergence estimators satisfy the expected $1/\sqrt{N}$ scaling.

### Stochastic Test Robustness for Small Systems

Enforcing the IV-2 temperature invariant ($|T_{\text{avg}} - 310| < 5$ K) on the miniature test systems used for unit and integration testing presents a finite-size statistical challenge distinct from production-scale simulations. For a solvated alanine dipeptide with $N_{\text{dof}} \approx 1000$ degrees of freedom, the equipartition theorem predicts instantaneous temperature fluctuations of $\sigma_{T_{\text{inst}}} = T\sqrt{2/N_{\text{dof}}} \approx 13.9$ K — individual frame temperatures fluctuate with a standard deviation nearly three times the 5 K tolerance. With the original test configuration (`nvt_duration_ps = 10.0`, `friction_per_ps = 5.0`), only $n = 10$ equilibrated temperature samples were available for averaging, and temporal autocorrelation imposed by the Langevin thermostat ($\tau_{\text{relax}} = 1/\gamma = 0.2$ ps) reduced the effective sample size to $n_{\text{eff}} \approx 7$. The corrected standard error:

$$\text{SE}_{\text{corr}}(\bar{T}) \approx \frac{\sigma_{T_{\text{inst}}}}{\sqrt{n_{\text{eff}}}} \approx \frac{13.9}{\sqrt{7}} \approx 5.3 \text{ K}$$

exceeded the 5 K tolerance, yielding a theoretical failure probability of $P(|\bar{T} - T_{\text{target}}| \geq 5 \text{ K}) \approx 2\Phi(-0.95) \approx 34\%$ per invocation — rendering the test unreliable for continuous integration despite correct underlying physics. Nondeterministic water placement during solvation (via OpenMM's `modeller.addSolvent`) further prevents seed-based reproducibility from stabilizing outcomes across runs.

**Resolution.** The equilibration parameters in all affected test configurations (`test_production.py`, `test_equilibrate.py`, `test_integration.py`) were adjusted:

| Parameter | Original | Updated | Rationale |
|-----------|----------|---------|----------|
| `nvt_duration_ps` | 10.0 | 100.0 | 10× more temperature samples in the equilibrated segment |
| `npt_duration_ps` | 40.0 | 100.0 | Consistent statistical margin for the NPT density invariant (IV-3) |
| `friction_per_ps` | 5.0 | 10.0 | Halves $\tau_{\text{relax}}$ to 0.1 ps, decorrelating 0.5 ps-spaced samples |

The increased friction renders consecutive temperature samples at 0.5 ps intervals nearly statistically independent ($\Delta t_{\text{save}} / \tau_{\text{relax}} = 5.0$). With approximately 100 frames in the equilibrated segment and near-independence between samples, the standard error drops below 2 K, and the predicted failure probability falls below 1%.

Empirical verification confirmed 16 consecutive passes of the previously intermittent test — an outcome with probability $P = 0.7^{16} \approx 0.3\%$ under the original failure rate — and a full regression of the complete test suite passed with zero failures. This fix is confined entirely to test-level equilibration configurations; the production IV-2 tolerance of $\pm 5$ K remains unchanged, as it is appropriate for production-scale systems with orders of magnitude more degrees of freedom.

---

## Installation & Usage

### Prerequisites

- Python 3.10+
- OpenMM ≥ 8.1 (may require [conda installation](http://docs.openmm.org/latest/userguide/application/01_getting_started.html) on some platforms)

### Setup

```bash
git clone https://github.com/ryanjosephkamp/SPINK7-KLK5-MD-Pipeline.git
cd SPINK7-KLK5-MD-Pipeline
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Verify Installation

```bash
python -m pytest tests/ -v
```

All 259 tests should pass with no failures.

### Running the Pipeline

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

<div style="page-break-after: always;"></div>

**Stage 5 — Umbrella Sampling:**
```bash
python scripts/run_umbrella.py
```

**Stage 6 — Analysis & Visualization:**
```bash
python scripts/run_analysis.py
```

### Interactive Notebooks

For exploratory analysis with real-time visualization:

```bash
jupyter notebook notebooks/
```

| Notebook | Purpose |
|----------|---------|
| `01_system_prep.ipynb` | Interactive structure preparation with visual inspection |
| `02_equilibration.ipynb` | Equilibration monitoring with real-time plots |
| `03_production_analysis.ipynb` | Production trajectory analysis |
| `04_smd_analysis.ipynb` | SMD work distributions and Jarzynski analysis |
| `05_umbrella_pmf.ipynb` | Umbrella Sampling results and WHAM PMF |
| `06_visualization.ipynb` | 3D molecular rendering and publication figures |

---

## Data Flow

```
RCSB / AlphaFold              PDB Files (.pdb)
        │                          │
        ▼                          ▼
   pdb_fetch.py              pdb_clean.py
        │                          │
        ▼                          ▼
   data/pdb/raw/              protonate.py
                                   │
                                   ▼
                              topology.py  ──►  solvate.py
                                                     │
        ┌──────────────────────────────────────────┘
        │
        ▼
   minimizer.py  ──►  equilibrate.py  ──►  production.py
                                                │
                    ┌───────────────────────────┤
                    │                           │
                    ▼                           ▼
               smd.py                     umbrella.py
                    │                           │
                    ▼                           ▼
             jarzynski.py                   wham.py
                    │                           │
                    └───────────┬───────────────┘
                                │
                                ▼
                    structural.py / contacts.py
                                │
                                ▼
                    visualization / plots
```

---

<div style="page-break-after: always;"></div>

## Dependencies

| Package | Version | Role |
|---------|---------|------|
| OpenMM | ≥ 8.1 | Core MD simulation engine |
| PDBFixer | ≥ 1.9 | Structural repair and standardization |
| MDTraj | ≥ 1.9.9 | Trajectory analysis and I/O |
| NumPy | ≥ 1.26, < 2.0 | Numerical computations |
| SciPy | ≥ 1.12 | Scientific computing |
| matplotlib | ≥ 3.8 | 2D plotting and visualization |
| py3Dmol | ≥ 2.0 | Interactive 3D molecular visualization |
| Biopython | ≥ 1.83 | PDB parsing and biological sequence utilities |
| PDB2PQR | ≥ 3.6 | Protonation state assignment |

---

<div style="page-break-after: always;"></div>

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

---

## Author

**Ryan Kamp**
Dept. of Computer Science, University of Cincinnati
[kamprj@mail.uc.edu](mailto:kamprj@mail.uc.edu) · [GitHub](https://github.com/ryanjosephkamp)

## License

[MIT](LICENSE) © 2026 Ryan Kamp
