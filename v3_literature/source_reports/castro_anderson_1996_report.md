# Source Report: S-17

## Bibliographic Information
- **ID:** S-17
- **Authors:** Castro MJM, Anderson S
- **Title:** Alanine Point-Mutations in the Reactive Region of Bovine Pancreatic Trypsin Inhibitor: Effects on the Kinetics and Thermodynamics of Binding to β-Trypsin and α-Chymotrypsin
- **Journal:** Biochemistry
- **Year:** 1996
- **Volume/Pages:** 35(35):11435–11446
- **DOI:** S0006-2960(96)00515-6
- **PDF Filename:** S-17.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category V: BPTI-Trypsin Benchmark

## Summary
Comprehensive alanine-scanning mutagenesis of the BPTI reactive region (15 single-point mutants at positions T11A through G36A) with full kinetic (kon, koff, Ki) and thermodynamic (ΔG, ΔH, TΔS) characterization against both trypsin and chymotrypsin. Confirms K_d = 5 × 10⁻¹⁴ M for wild-type BPTI-trypsin. The P1 mutation K15A causes the most dramatic effect: ΔΔG ≈ 10 kcal/mol for trypsin (K_d increases to 1.4 × 10⁻⁶ M) with >200-fold reduction in kon (electrostatic steering). Demonstrates that binding affinity is predominantly determined by koff, and reveals enthalpy-entropy compensation among mutants.

## Extracted Features

### Feature 1: Wild-Type BPTI-Trypsin Kinetics
- **Type:** Kinetic/Thermodynamic
- **Value:** Ki = 5 × 10⁻¹⁴ M; kon = (9.9 ± 2.5) × 10⁵ M⁻¹s⁻¹; koff = 5 × 10⁻⁸ s⁻¹
- **Conditions:** 50 mM Tris-HCl, 20 mM CaCl₂, 0.005% Triton X-100, pH 8.2, 22°C
- **Confidence:** High
- **Benchmarkable:** Yes — confirms S-16 benchmark value

### Feature 2: P1 Mutant (K15A) BPTI-Trypsin
- **Type:** Kinetic/Thermodynamic
- **Value:** Ki = (1.4 ± 0.1) × 10⁻⁶ M; kon = (4.3 ± 0.6) × 10³ M⁻¹s⁻¹; koff = 4.2 × 10⁻⁵ s⁻¹; ΔΔG ≈ 10 kcal/mol
- **Conditions:** pH 8.2, 22°C
- **Confidence:** High
- **Benchmarkable:** Yes — key FEP/computational alanine scanning benchmark

### Feature 3: Complete Alanine Scanning for Trypsin (Table 1)
- **Type:** Kinetic
- **Value:** 10 mutants measured: T11A (<5×10⁻¹¹), G12A (8.3×10⁻¹¹), P13A (<5×10⁻¹¹), K15A (1.4×10⁻⁶), R17A (<5×10⁻¹¹), I18A (2.4×10⁻¹⁰), I19A (<5×10⁻¹¹), R20A (<5×10⁻¹¹), G36A (2.1×10⁻¹⁰). All with kon and koff values.
- **Conditions:** pH 8.2, 22°C
- **Confidence:** High (some are upper bounds)
- **Benchmarkable:** Yes — per-residue ΔΔG values for computational alanine scanning validation

### Feature 4: Complete Alanine Scanning for Chymotrypsin (Table 1)
- **Type:** Kinetic
- **Value:** wt Ki = (1.1 ± 0.1) × 10⁻⁸ M; K15A = (3.3 ± 0.1) × 10⁻⁷ M; I18A = (1.2 ± 0.2) × 10⁻⁷ M; plus 12 additional mutants with kon/koff/Ki
- **Conditions:** pH 8.2, 22°C
- **Confidence:** High
- **Benchmarkable:** Yes — cross-protease selectivity benchmark

### Feature 5: Thermodynamic Decomposition for Chymotrypsin (Table 2)
- **Type:** Thermodynamic
- **Value:** wt ΔG = −10.7 kcal/mol, ΔH = +2.5 kcal/mol, −TΔS = −13.2 kcal/mol; G12A ΔG = −10.1, ΔH = −0.6, −TΔS = −9.5; K15A ΔG = −8.7, ΔH = +3.4, −TΔS = −12.1; V34A ΔG = −10.7, ΔH = −2.7, −TΔS = −8.0; plus 14 total mutants
- **Conditions:** ITC, pH 8.2, 22°C
- **Confidence:** High
- **Benchmarkable:** Yes — ΔH/ΔS decomposition for force field validation

### Feature 6: Entropy-Driven Binding and Enthalpy-Entropy Compensation
- **Type:** Thermodynamic
- **Value:** ΔH for wt BPTI-chymotrypsin is +2.5 kcal/mol (endothermic); binding is entirely entropy-driven (−TΔS = −13.2 kcal/mol). Compensation observed: ΔΔH and ΔΔ(−TΔS) anti-correlated across mutants.
- **Conditions:** ITC data
- **Confidence:** High
- **Benchmarkable:** Yes — critical for validating force field energetics

### Feature 7: Thermal Stability of Mutants (Table 3)
- **Type:** Thermodynamic
- **Value:** wt Tm = 94.2°C, ΔHm = 84.3 kcal/mol; K15A ΔTm = −4.2°C; F33A ΔTm = −16.5°C; G12A ΔTm = −16.1°C
- **Conditions:** DSC, pH 3.0
- **Confidence:** High
- **Benchmarkable:** Partially — protein stability, not binding

### Feature 8: kon Dominance by Electrostatic Steering
- **Type:** Kinetic
- **Value:** K15A reduces kon by >200-fold for trypsin (10⁶ → 10³ M⁻¹s⁻¹) but only ~3-fold for chymotrypsin, demonstrating electrostatic steering
- **Conditions:** pH 8.2, 22°C
- **Confidence:** High
- **Benchmarkable:** Yes — electrostatic steering effect testable by Brownian dynamics or long MD

## PDB Codes Referenced
- **5PTI** (free BPTI)
- **2PTC** (trypsin-BPTI complex)

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** The richest per-residue mutagenesis dataset available for any canonical protease-inhibitor system. The 15-mutant alanine scan with full kinetics (kon, koff, Ki) against two proteases plus thermodynamic decomposition (ΔG, ΔH, TΔS) for chymotrypsin provides the most comprehensive validation dataset for computational alanine scanning and FEP calculations. The entropy-driven binding result is important for force field validation. PDB structures 5PTI and 2PTC are readily available.
- **Key Limitation:** BPTI is Kunitz-type, not Kazal-type. Direct transfer to SPINK7-KLK5 requires structural analogy. However, the methodology and validation framework are directly applicable.
