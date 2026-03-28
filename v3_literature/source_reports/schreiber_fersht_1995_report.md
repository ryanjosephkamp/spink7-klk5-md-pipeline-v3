# Source Report: S-22

## Bibliographic Information
- **ID:** S-22
- **Authors:** Schreiber G, Fersht AR
- **Title:** Energetics of Protein-Protein Interactions: Analysis of the Barnase-Barstar Interface by Single Mutations and Double Mutant Cycles
- **Journal:** Journal of Molecular Biology
- **Year:** 1995
- **Volume/Pages:** 248:478–486
- **DOI:** N/A
- **PDF Filename:** S-22.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VII: PPI Energetics

## Summary
Landmark study establishing the double mutant cycle (DMC) methodology for analyzing protein-protein interaction energetics at the single-residue level, using the barnase-barstar system. Wild-type Ka = 10¹⁴ M⁻¹ (Kd = 0.01 pM). Provides coupling energies (ΔΔGint) for ~45 single and double mutant pairs, demonstrating that coupling energy decreases with inter-residue distance and that residues within ~7 Å interact cooperatively. Identifies electrostatic hot spots (Asp39b*–Arg87bn: 7.0 kcal/mol coupling) at the interface.

## Extracted Features

### Feature 1: Wild-Type Barnase-Barstar Affinity
- **Type:** Thermodynamic
- **Value:** Ka = 10¹⁴ M⁻¹ (Kd = 10⁻¹⁴ M = 0.01 pM); ΔG = 19.0 kcal/mol
- **Conditions:** 50 mM Tris-HCl pH 8.0, 25°C
- **Confidence:** High
- **Benchmarkable:** Yes — ultra-tight binding reference system

### Feature 2: Association/Dissociation Kinetics
- **Type:** Kinetic
- **Value:** kon = 3.7 × 10⁸ M⁻¹s⁻¹; koff = 3.7 × 10⁻⁶ s⁻¹
- **Conditions:** pH 8.0, 25°C
- **Confidence:** High
- **Benchmarkable:** Yes — near diffusion-limited association

### Feature 3: Double Mutant Cycle Coupling Energies (Table 1)
- **Type:** Thermodynamic
- **Value:** ~45 mutant pairs analyzed. Strongest coupling: Asp39b*–Arg87bn = 7.0 kcal/mol; Asp39b*–His102bn = 5.5 kcal/mol; Asp39b*–Arg83bn = 5.5 kcal/mol; Asp39b*–Lys27bn = 5.0 kcal/mol. Medium: Asp35b*–Arg59bn = 3.0–4.0 kcal/mol; Tyr29b*–His102bn = 3.2 kcal/mol. Weak: Gly31b*–Arg83bn = 1.3 kcal/mol; Asn33b*–Lys27bn = 0.7 kcal/mol.
- **Conditions:** pH 8.0, 25°C; stopped-flow kinetics
- **Confidence:** High (±0.2 kcal/mol, 2σ)
- **Benchmarkable:** Yes — per-residue coupling energies for FEP/DMC validation

### Feature 4: Distance-Coupling Correlation
- **Type:** Structural/Thermodynamic
- **Value:** Coupling energy decreases with inter-residue distance; residues within ~7 Å show cooperative interactions
- **Conditions:** Correlation analysis
- **Confidence:** High
- **Benchmarkable:** Yes — testable by computational DMC simulations

### Feature 5: His102bn pKa Shift
- **Type:** Thermodynamic
- **Value:** pKa shifts from 6.3 (free barnase) to ~5 (complexed)
- **Conditions:** pH titration + kinetics
- **Confidence:** Moderate
- **Benchmarkable:** Partially — constant-pH MD could reproduce this

## PDB Codes Referenced
- **1brs** (implicit via Buckle et al. 1994 — barnase-barstar complex at 2.0 Å)

## Usefulness Assessment
- **Overall Rating:** Medium
- **Rationale:** Establishes the double mutant cycle methodology for decomposing PPI binding energetics — directly applicable to SPINK7-KLK5 interface analysis. The comprehensive dataset of ~45 coupling energies provides a validation framework for FEP and computational mutagenesis approaches. However, barnase-barstar is not a protease-inhibitor system, so data cannot be directly transferred. The methodology and analytical framework are the primary contributions.
- **Key Limitation:** Different protein system (barnase-barstar, not protease-inhibitor). Provides methodology rather than directly applicable benchmark data for SPINK-KLK.
