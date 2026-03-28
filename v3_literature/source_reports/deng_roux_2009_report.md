# Source Report: S-25

## Bibliographic Information
- **ID:** S-25
- **Authors:** Deng Y, Roux B
- **Title:** Computations of Standard Binding Free Energies with Molecular Dynamics Simulations
- **Journal:** Journal of Physical Chemistry B
- **Year:** 2009
- **Volume/Pages:** 113:2234–2246
- **DOI:** 10.1021/jp807701h
- **PDF Filename:** S-25.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VIII: Computational Benchmarking

## Summary
Centennial Feature Article describing the rigorous statistical mechanical framework for computing standard binding free energies from molecular dynamics simulations. Presents two approaches: the Double Decoupling Method (DDM) with alchemical transformations and restraints, and the Potential of Mean Force (PMF) method with separation along a reaction coordinate. Demonstrates accuracy on multiple test systems: T4 lysozyme mutants (~1.9 kcal/mol RMS error), FKBP12-ligand (computed −10.2 vs. experimental −10.9 kcal/mol), SH2-peptide, and cytochrome P450-camphor. Expected accuracy: ~1–2.5 kcal/mol RMS error for well-behaved systems.

## Extracted Features

### Feature 1: Double Decoupling Method (DDM) Framework
- **Type:** Methodological
- **Value:** ΔG°_b = ΔΔG_int + ΔΔG°_t + ΔΔG_r + ΔΔG_c; force decomposition: ΔΔG_int = ΔΔG_rep + ΔΔG_dis + ΔΔG_elec
- **Conditions:** Alchemical transformation with restraints
- **Confidence:** High — rigorous statistical mechanics
- **Benchmarkable:** Framework — defines the computational approach

### Feature 2: T4 Lysozyme Binding Free Energies (Table 1)
- **Type:** Thermodynamic
- **Value:** L99A/benzene: computed −5.96 vs. experimental −5.19 kcal/mol; M102Q-L99A/phenol: computed −5.64 vs. experimental −5.55 kcal/mol; 13-ligand retrospective RMS error ~1.9 kcal/mol
- **Conditions:** CHARMM force field; PDB 183L
- **Confidence:** High
- **Benchmarkable:** Yes — establishes expected accuracy for small-molecule binding

### Feature 3: FKBP12-Ligand 8 Binding (Table 2)
- **Type:** Thermodynamic
- **Value:** ΔG°_bind computed = −10.2 kcal/mol vs. experimental = −10.9 kcal/mol; ΔΔG_rep = −1.1, ΔΔG_dis = −21.1, ΔΔG_elec = −3.7, ΔΔG_c = 6.9 kcal/mol
- **Conditions:** CHARMM + GAFF/AM1-BCC charges
- **Confidence:** High
- **Benchmarkable:** Yes — medium-affinity system benchmark

### Feature 4: Protein-Peptide Binding
- **Type:** Thermodynamic
- **Value:** SH2-pYEEI: computed −8.8 vs. experimental −7.1 kcal/mol; cytochrome P450-camphor: computed −8.25 vs. experimental −7.75 kcal/mol
- **Conditions:** PMF + DDM approaches
- **Confidence:** Moderate
- **Benchmarkable:** Yes — validates approach for protein-peptide systems

### Feature 5: Expected Accuracy Benchmarks
- **Type:** Methodological
- **Value:** RMS error 1–2.5 kcal/mol for well-behaved systems; ~60% of favorable interaction energy opposed by unfavorable motional freedom loss; translational freedom loss ~4.4 kcal/mol (standard state)
- **Conditions:** Multiple test systems
- **Confidence:** High — consistent across systems
- **Benchmarkable:** Framework — defines accuracy expectations for our pipeline

### Feature 6: Force Field Recommendations
- **Type:** Methodological
- **Value:** CHARMM, AMBER, OPLS recommended; GAFF with AM1-BCC charges for small molecules; long-range electrostatics via PME or reaction field; mixed explicit/implicit solvent (SSBP, GSBP) possible
- **Conditions:** Literature review and benchmarking
- **Confidence:** High
- **Benchmarkable:** Framework — informs force field selection for SPINK7-KLK5

## PDB Codes Referenced
183L (T4 lysozyme L99A-indene complex)

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** Essential methodological reference for the computational approach to binding free energy calculations. The DDM and PMF frameworks are directly applicable to SPINK7-KLK5. The accuracy benchmarks (1–2.5 kcal/mol RMS error) establish realistic expectations for our pipeline's computed ΔG values relative to the experimental Ki = 132 nM (ΔG ≈ −9.4 kcal/mol). The force decomposition framework enables understanding of individual contributions. Force field recommendations inform setup.
- **Key Limitation:** All test systems are protein-small molecule or protein-peptide; no protein-protein binding free energy calculations demonstrated. Extrapolation to the larger SPINK7-KLK5 interface involves additional challenges (larger conformational changes, more complex solvation).
