# Source Report: S-26

## Bibliographic Information
- **ID:** S-26
- **Authors:** Mobley DL, Gilson MK
- **Title:** Predicting Binding Free Energies: Frontiers and Benchmarks
- **Journal:** Annual Review of Biophysics
- **Year:** 2017
- **Volume/Pages:** 46:531–558
- **DOI:** 10.1146/annurev-biophys-070816-033654
- **PDF Filename:** S-26.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VIII: Computational Benchmarking

## Summary
Comprehensive review defining the framework for benchmarking binding free energy calculations. Distinguishes between "hard" benchmarks (well-controlled systems with precise experimental data and convergent calculations) and "soft" benchmarks (complex systems where convergence is uncertain). Presents curated benchmark datasets for cucurbit[7]uril (CB7) host-guest systems, GDCC (Octa Acid) systems, and T4 lysozyme mutant cavities. Establishes that current methods achieve 1–2 kcal/mol RMS error for well-behaved systems, and that 2 kcal/mol accuracy is sufficient for practical prioritization.

## Extracted Features

### Feature 1: Current Accuracy Benchmarks
- **Type:** Methodological
- **Value:** Free energy calculations: 1–2 kcal/mol RMS error currently achievable; 2 kcal/mol sufficient for ~3-fold reduction in compounds needing synthesis; target precision for hard benchmarks: ~0.1 kcal/mol uncertainty
- **Conditions:** Multiple test systems
- **Confidence:** High — consensus from multiple groups
- **Benchmarkable:** Framework — defines accuracy expectations

### Feature 2: CB7 Host-Guest Benchmark Data (Tables 1–2)
- **Type:** Thermodynamic
- **Value:** Adamantane series ΔG: −5.99 to −18.22 kcal/mol; aromatic series ΔG: −6.31 to −12.62 kcal/mol. SAMPL4 results: R² = 0.1–0.8, RMS errors 1.9–4.9 kcal/mol (best R² = 0.92).
- **Conditions:** 50 mM NaOAc, D₂O, pH 4.74, 298 K
- **Confidence:** High — ITC measurements
- **Benchmarkable:** Yes — host-guest benchmark

### Feature 3: T4 Lysozyme L99A Benchmark Data (Table 5)
- **Type:** Thermodynamic
- **Value:** Benzene ΔG = −5.19 ± 0.16; toluene −5.52 ± 0.04; ethylbenzene −5.76 ± 0.07; propylbenzene −6.55 ± 0.02; butylbenzene −6.70 ± 0.02; p-xylene −4.67 ± 0.06; benzofuran −5.46 ± 0.03 kcal/mol
- **Conditions:** 302 K, 50 mM NaOAc, pH 5.5
- **Confidence:** High — well-characterized system
- **Benchmarkable:** Yes — model cavity benchmark with PDB codes

### Feature 4: Sampling Timescale Requirements
- **Type:** Methodological
- **Value:** CB7 guest binding-mode flips: 0.07 flips/ns; CB7 wetting/dewetting: ~50 ns; GDCC water motions: ~5 ns; T4L in-plane transitions: 1–10 ns; T4L out-of-plane: hundreds of ns; T4L helix F motion: ~50 ns
- **Conditions:** Various systems
- **Confidence:** High
- **Benchmarkable:** Framework — guides simulation length requirements

### Feature 5: Salt/Buffer Sensitivity
- **Type:** Methodological
- **Value:** CB7 salt-dependent ΔG variations: 2.5–2.8 kcal/mol (50 mM NaOAc vs 100 mM Na₃PO₄); GDCC salt shifts: 100 mM NaClO₄/NaClO₃/NaSCN → ΔH shift ~10 kcal/mol, ΔG shift ~2 kcal/mol
- **Conditions:** Various buffers/salts
- **Confidence:** High
- **Benchmarkable:** Warning — salt/buffer conditions must be carefully matched

## PDB Codes Referenced
181L, 4W52, 4W53, 1NHB, 4W54, 4W55, 186L, 4W57, 4W59, 187L, 182L, 1LI2, 1XEP, 3HU8, 3HUK, 3HUA, 2RBN, 1LI3

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** Defines the conceptual and practical framework for benchmarking binding free energy calculations. The hard/soft benchmark distinction is directly applicable — SPINK7-KLK5 is inherently a "soft" benchmark due to the large PPI interface and conformational complexity. The accuracy expectations (1–2 kcal/mol achievable, 2 kcal/mol sufficient) calibrate expectations for our pipeline. Salt/buffer sensitivity warnings are critical for ensuring simulation conditions match the experimental Ki = 132 nM measurement conditions.
- **Key Limitation:** All benchmark systems are host-guest or protein-small molecule; no protein-protein benchmarks. Directly applicable framework, but the complexity gap to SPINK7-KLK5 is large.
