# Source Report: S-27

## Bibliographic Information
- **ID:** S-27
- **Authors:** Park S, Khalili-Araghi F, Tajkhorshid E, Schulten K
- **Title:** Free Energy Calculation from Steered Molecular Dynamics Simulations Using Jarzynski's Equality
- **Journal:** Journal of Chemical Physics
- **Year:** 2003
- **Volume/Pages:** 119(6):3559–3566
- **DOI:** 10.1063/1.1590311
- **PDF Filename:** S-27.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VIII: Computational Benchmarking

## Summary
Methodological paper demonstrating the extraction of free energy profiles from steered molecular dynamics (SMD) simulations using Jarzynski's equality. Uses deca-alanine in vacuum as the test system, comparing irreversible work from SMD at different pulling speeds (10 and 100 Å/ns) with reversible reference PMF (0.1 Å/ns). Shows that the 2nd-order cumulant expansion with ~10 trajectories provides the most accurate estimates. Provides practical guidance: at 10 Å/ns, mean error is 1.6 kcal/mol (7.6%) from 10 trajectories; at 100 Å/ns, error increases to 6.7 kcal/mol (31%).

## Extracted Features

### Feature 1: Reversible PMF Reference
- **Type:** Thermodynamic
- **Value:** ΔF_total = 21.4 kcal/mol (helix → extended coil); PMF minimum at ξ ≈ 15.2 Å; reversible work standard deviation < 0.5 k_BT
- **Conditions:** Deca-alanine vacuum, 300 K, Langevin dynamics, v = 0.1 Å/ns (200 ns), CHARMM22, NAMD
- **Confidence:** High — converged reversible simulation
- **Benchmarkable:** Yes — internal methodology benchmark

### Feature 2: SMD Accuracy at 10 Å/ns
- **Type:** Methodological
- **Value:** Irreversible work up to 2.7 kcal/mol (4.5 k_BT); work standard deviation 1.9 kcal/mol (3.1 k_BT); mean error from 10 trajectories = 1.6 kcal/mol (7.6%)
- **Conditions:** Deca-alanine, k = 500 pN/Å, 10 trajectories
- **Confidence:** High
- **Benchmarkable:** Methodology reference — establishes accuracy vs. speed tradeoff

### Feature 3: SMD Accuracy at 100 Å/ns
- **Type:** Methodological
- **Value:** Irreversible work up to 18.8 kcal/mol (31.3 k_BT); work standard deviation 4.3 kcal/mol (7.1 k_BT); mean error 6.7 kcal/mol (31%); all total work values in 35–50 kcal/mol range
- **Conditions:** Deca-alanine, k = 500 pN/Å, 10 trajectories
- **Confidence:** High
- **Benchmarkable:** Methodology — shows limits of fast pulling

### Feature 4: Optimal Estimator Selection
- **Type:** Methodological
- **Value:** 2nd-order cumulant expansion yields most accurate estimates with ~10 trajectories; unbiased variance estimator (Eq. 19) improves small-N results; stiff-spring approximation correction < 0.5 kcal/mol for k = 500 pN/Å
- **Conditions:** Comparison of estimators
- **Confidence:** High
- **Benchmarkable:** Framework — guides analysis of SMD-based ΔG

### Feature 5: Umbrella Sampling Comparison
- **Type:** Methodological
- **Value:** 10 harmonic biasing potentials, A = 70 pN/Å, positions spanning 13.4–33.0 Å; WHAM combination; efficiency comparable to SMD/Jarzynski
- **Conditions:** Same system
- **Confidence:** High
- **Benchmarkable:** Methodology — validates SMD against US

## PDB Codes Referenced
None

## Usefulness Assessment
- **Overall Rating:** Medium
- **Rationale:** Provides the theoretical and practical foundation for SMD-based free energy calculations using Jarzynski's equality — applicable to protease-inhibitor unbinding trajectories. The accuracy-vs-speed tradeoffs and optimal estimator selection guide experimental design. However, the test system (deca-alanine in vacuum, 104 atoms) is far simpler than SPINK7-KLK5 in explicit solvent, and the much larger dissipation expected in PPI systems will require substantially more trajectories.
- **Key Limitation:** Vacuum test system with minimal solvent friction. Direct application to SPINK7-KLK5 (>30,000 atoms, explicit water, large BSA) will involve much larger irreversible work and poorer convergence. The methodology is sound but scaling to PPI systems requires careful adaptation.
