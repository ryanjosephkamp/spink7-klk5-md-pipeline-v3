# Source Report: S-28

## Bibliographic Information
- **ID:** S-28
- **Authors:** Gumbart JC, Roux B, Chipot C
- **Title:** Standard Binding Free Energies from Computer Simulations: What Is the Best Strategy?
- **Journal:** Journal of Chemical Theory and Computation
- **Year:** 2013
- **Volume/Pages:** 9(1):794–802
- **DOI:** 10.1021/ct3008099
- **PDF Filename:** S-28.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VIII: Computational Benchmarking

## Summary
Gold-standard methodological paper comparing alchemical and PMF-based routes for computing absolute binding free energies. Uses the SH3-p41 peptide system (PDB 1bbz) as the test case, achieving remarkable convergence: PMF = −7.8 kcal/mol, alchemical = −7.7 kcal/mol, experimental = −7.99 kcal/mol. Both approaches agree within ~0.1 kcal/mol of each other. Provides complete thermodynamic decomposition (Table 1) showing contributions from conformational restraints, orientational restraints, and separation/coupling. The PMF-based separation approach using REMD-US is the most computationally efficient (~115 ns) and directly applicable to protein-protein systems.

## Extracted Features

### Feature 1: Experimental Binding Free Energy
- **Type:** Thermodynamic
- **Value:** ΔG_bind = −7.99 kcal/mol (SH3-p41 peptide)
- **Conditions:** From Pisabarro & Serrano 1996
- **Confidence:** High
- **Benchmarkable:** Yes — validation target

### Feature 2: Computed Binding Free Energies
- **Type:** Thermodynamic
- **Value:** PMF (REMD-US): −7.8 ± 0.9 kcal/mol; PMF (ABF): −7.7 ± 0.9 kcal/mol; Alchemical: −7.7 ± 1.0 kcal/mol; MM/PBSA: −2.6 kcal/mol (poor accuracy)
- **Conditions:** CHARMM22/27, TIP3P, NAMD 2.9, NPT 300 K 1 atm, PME, 2/4 fs timestep
- **Confidence:** High — three independent methods agree
- **Benchmarkable:** Yes — method accuracy comparison

### Feature 3: Separation PMF Profile
- **Type:** Thermodynamic
- **Value:** Deep well at r = 20.3 Å of −17.5 kcal/mol; S*I* product varies by 0.3 kcal/mol across r* values; reference point r* = 35 Å
- **Conditions:** REMD-US, 37 windows (17–43 Å), k = 2.5 kcal/mol·Å²
- **Confidence:** High
- **Benchmarkable:** Yes — PMF shape/depth is a direct output of protease-inhibitor MD

### Feature 4: Thermodynamic Decomposition (Table 1)
- **Type:** Thermodynamic
- **Value:** Key contributions (PMF route): ΔG_c,site (RMSD) = −3.6, ΔG_separation = −14.5, ΔG_c,bulk = +5.8 kcal/mol. Key contributions (alchemical): ΔG_decouple,site = +35.9, ΔG_couple,bulk = −53.3, ΔG_c,bulk = +6.1 kcal/mol.
- **Conditions:** Full decomposition into conformational, orientational, and separation/coupling terms
- **Confidence:** High — rigorous decomposition
- **Benchmarkable:** Yes — component analysis framework

### Feature 5: Computational Cost Comparison
- **Type:** Methodological
- **Value:** PMF route (REMD-US): 115 ns total; PMF route (ABF): 155 ns; Alchemical route: 280 ns. REMD-US: 37 windows × 0.5 ns = 18.5 ns for separation PMF alone.
- **Conditions:** NAMD 2.9, ~32,000 atom system
- **Confidence:** High
- **Benchmarkable:** Framework — guides computational budget

### Feature 6: Restraint Framework (6 Degrees of Freedom)
- **Type:** Methodological
- **Value:** Six geometric restraints: r (separation), θ (polar angle), ϕ (azimuthal angle), Θ (Euler angle 1), Φ (Euler angle 2), Ψ (Euler angle 3). Plus RMSD conformational restraint. Spring constants: k = 6–15 kcal/mol·deg² (angular), 250–1000 kcal/mol·Å² (RMSD), 60–250 kcal/mol·Å² (separation).
- **Conditions:** Colvars module in NAMD
- **Confidence:** High
- **Benchmarkable:** Framework — template for SPINK7-KLK5 collective variable definitions

### Feature 7: REMD-US Parameters
- **Type:** Methodological
- **Value:** 37 windows, 17–43 Å, spacing 0.5 Å (20–30 Å) and 1 Å elsewhere; harmonic restraint k = 2.5 kcal/mol·Å²; exchanges every 1000 timesteps; WHAM reconstruction; 0.5 ns per window
- **Conditions:** SH3-p41 system
- **Confidence:** High
- **Benchmarkable:** Framework — template for SPINK7-KLK5 US parameters

### Feature 8: MM/PBSA Accuracy Comparison
- **Type:** Methodological
- **Value:** MM/PBSA: ΔG = −2.6 kcal/mol vs. experimental −7.99 kcal/mol (error > 5 kcal/mol). Rigorous FE methods: error < 0.3 kcal/mol.
- **Conditions:** Same system, same trajectories
- **Confidence:** High
- **Benchmarkable:** Yes — demonstrates inadequacy of MM/PBSA for absolute ΔG

### Feature 9: ABF Sensitivity Warning
- **Type:** Methodological
- **Value:** ABF very sensitive to initial conditions; ABF bin widths: 1° (angular), 0.05 Å (RMSD), 0.1 Å (separation); soft positional restraints: k = 0.03 kcal/mol·deg²
- **Conditions:** ABF exploration
- **Confidence:** High — documented pitfall
- **Benchmarkable:** Warning — ABF requires careful setup

## PDB Codes Referenced
**1bbz** (SH3 domain of Abl kinase complexed with p41 peptide)

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** The most directly applicable methodological reference for computing SPINK7-KLK5 binding free energy. The PMF-based separation approach with REMD-US is the recommended strategy for protein-protein/protein-peptide systems (interfacial binding mode, physically meaningful unbinding pathway). Provides: (1) complete protocol with simulation parameters; (2) known accuracy (~1 kcal/mol uncertainty); (3) computational cost estimates (~115 ns for PMF route); (4) comparison showing MM/PBSA inadequacy; (5) restraint framework for defining collective variables; (6) ABF sensitivity warning.
- **Key Limitation:** Test system (SH3-p41, 10-residue peptide, ΔG = −8 kcal/mol) is smaller and less tightly bound than SPINK7-KLK5 (~56 residue inhibitor, ΔG ≈ −9.4 kcal/mol). Scaling to the larger SPINK7-KLK5 interface will require more US windows, longer sampling per window, and careful treatment of the larger conformational space.
