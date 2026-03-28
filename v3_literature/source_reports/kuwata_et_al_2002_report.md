# Source Report: S-20

## Bibliographic Information
- **ID:** S-20
- **Authors:** Kuwata K, Hirota M, Shimizu H, Nakae M, Nishihara S, Takimoto A, Mitsushima K, Kikuchi N, Endo K, Inoue M, Ogawa M
- **Title:** Functional Analysis of Recombinant Pancreatic Secretory Trypsin Inhibitor Protein with Amino-acid Substitution
- **Journal:** Journal of Gastroenterology
- **Year:** 2002
- **Volume/Pages:** 37:928–934
- **DOI:** N/A
- **PDF Filename:** S-20.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VI: SPINK1-Trypsin

## Summary
Functional analysis of recombinant wild-type and N34S mutant PSTI/SPINK1 proteins. Establishes Ki ≈ 7.6–8.7 × 10⁻⁹ M for SPINK1-trypsin interaction (Green & Work method). Demonstrates that the clinically important N34S mutation (found in 9.1% of pancreatitis patients vs. 0.7% controls) does NOT alter trypsin inhibitor function under any tested condition (pH 5–9, ±CaCl₂, ±EDTA, dissociation kinetics). This is a critical negative result suggesting N34S pathogenicity operates through a non-functional mechanism.

## Extracted Features

### Feature 1: Wild-Type SPINK1-Trypsin Ki
- **Type:** Thermodynamic
- **Value:** Recombinant wt Ki ≈ 7.6 × 10⁻⁹ M; natural PSTI Ki ≈ 8.7 × 10⁻⁹ M
- **Conditions:** 0.1 M Tris-HCl pH 8.0, 0.02 M CaCl₂, 0.01% Triton X-100, 37°C; Green & Work method
- **Confidence:** Moderate — single method, no full kinetic characterization
- **Benchmarkable:** Yes — ΔG ≈ −11.1 kcal/mol, Kazal-family benchmark

### Feature 2: N34S Mutation Functional Equivalence
- **Type:** Kinetic/Functional
- **Value:** N34S shows no functional difference from wild-type under all conditions tested (pH 5–9, ±CaCl₂, ±EDTA, dissociation kinetics up to 24h)
- **Conditions:** Multiple buffer/pH conditions
- **Confidence:** High — comprehensive negative result
- **Benchmarkable:** Yes — mutation should yield ΔΔG ≈ 0 in computational mutagenesis

### Feature 3: pH Independence of Inhibition
- **Type:** Kinetic
- **Value:** Both wt and N34S retain inhibitory activity across pH 5–9 (MES, Tris, CHES buffers)
- **Conditions:** Multiple pH values, 37°C
- **Confidence:** High
- **Benchmarkable:** Context — validates use of standard simulation pH conditions

### Feature 4: Calcium Independence
- **Type:** Kinetic
- **Value:** Inhibitory activity unchanged with 0.02 M CaCl₂ or 0.01 M EDTA
- **Conditions:** ±CaCl₂, ±EDTA
- **Confidence:** High
- **Benchmarkable:** Context — SPINK1 does not require Ca²⁺ for function (unlike some serine proteases)

### Feature 5: N34S Clinical Prevalence
- **Type:** Epidemiological
- **Value:** N34S in pancreatitis patients: 83/908 (9.1%); controls: 11/1474 (0.7%)
- **Conditions:** Japanese population study
- **Confidence:** High — large sample
- **Benchmarkable:** No — clinical data

## PDB Codes Referenced
None

## Usefulness Assessment
- **Overall Rating:** Medium
- **Rationale:** Provides a Kazal-family (SPINK1) Ki benchmark value for trypsin (7.6–8.7 nM). The N34S negative result is a useful computational control — FEP of this mutation should predict ΔΔG ≈ 0. The broad pH and CaCl₂ independence validates standard simulation conditions. However, no structural data, no kon/koff kinetics, and no thermodynamic decomposition limit the quantitative benchmarking value.
- **Key Limitation:** Single Ki measurement method (Green & Work). No structural coordinates, no detailed kinetics, no thermodynamic parameters (ΔH, ΔS). Limited value for quantitative MD validation beyond the Ki value itself.
