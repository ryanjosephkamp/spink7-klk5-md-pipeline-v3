# Source Report: S-24

## Bibliographic Information
- **ID:** S-24
- **Authors:** Lo Conte L, Chothia C, Janin J
- **Title:** The Atomic Structure of Protein-Protein Recognition Sites
- **Journal:** Journal of Molecular Biology
- **Year:** 1999
- **Volume/Pages:** 285:2177–2198
- **DOI:** N/A
- **PDF Filename:** S-24.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VII: PPI Energetics

## Summary
Landmark structural analysis of 75 protein-protein complexes defining the quantitative anatomy of PPI interfaces. Establishes that standard-size interfaces bury 1600 ± 400 Å² and that protease-inhibitor interfaces average 1530 Å² with ~10 H-bonds and ~15 interface water molecules. Critically, includes 1tgs (trypsinogen-PSTI/SPINK1) and 2kai (kallikrein A-PTI) — the two most directly relevant structures for SPINK-KLK benchmarking. Shows protease-inhibitor interfaces are 61% non-polar with perfect close-packing (V/Vo = 1.00 ± 0.01).

## Extracted Features

### Feature 1: Protease-Inhibitor Interface Area
- **Type:** Structural
- **Value:** Mean BSA = 1530 Å² (σ = 170) for protease-inhibitor complexes. Specific entries: 1tgs = 1730 Å², 2kai = 1440 Å², 2ptc = 1570 Å², 1cho = 1400 Å², 1acb = 1460 Å², 1ppf = 1500 Å²
- **Conditions:** Probe radius 1.4 Å, compiled from crystal structures
- **Confidence:** High
- **Benchmarkable:** Yes — direct MD validation target for SPINK7-KLK5 interface area

### Feature 2: Interface Hydrogen Bonds
- **Type:** Structural
- **Value:** Protease-inhibitor mean: 9.5 H-bonds. 1tgs = 10, 2kai = 10, 2ptc = 13, 1cho = 10. One H-bond per ~170 Å² of interface.
- **Conditions:** Standard H-bond criteria
- **Confidence:** High
- **Benchmarkable:** Yes — H-bond count from MD must match

### Feature 3: Interface Water Molecules
- **Type:** Structural
- **Value:** Protease-inhibitor mean: 14.6 waters. 1tgs = 18, 2ptc = 11, 1cho = 16, 1acb = 8. One water per ~100 Å² of interface.
- **Conditions:** Crystal structure ordered waters
- **Confidence:** High
- **Benchmarkable:** Yes — solvation structure from MD

### Feature 4: Interface Chemical Composition
- **Type:** Structural
- **Value:** Protease-inhibitor: 61% non-polar. All interfaces: Arg contributes most (10%), followed by Tyr, Asp, Leu, Ile. Contact atoms 51%, buried 35%, rim 14%.
- **Conditions:** Classification by SASA change
- **Confidence:** High
- **Benchmarkable:** Yes — composition check for MD models

### Feature 5: Interface Packing Density
- **Type:** Structural
- **Value:** Protease-inhibitor V/Vo = 1.00 ± 0.01 (perfect close-packing, same as protein interior). All complexes: V/Vo = 1.01 ± 0.03.
- **Conditions:** Voronoi analysis of crystal structures
- **Confidence:** High
- **Benchmarkable:** Yes — packing density from MD trajectories

### Feature 6: Polar Interaction Types
- **Type:** Structural
- **Value:** 24% main chain–main chain, 40% main chain–side chain, 36% side chain–side chain. 30% charged; 13% salt bridges
- **Conditions:** Compiled across 75 complexes
- **Confidence:** High
- **Benchmarkable:** Yes — interaction type distribution from MD

### Feature 7: Standard Interface Size Distribution
- **Type:** Structural
- **Value:** Standard: 1600 ± 400 Å² (52/75 = 70%); large: 2000–4660 Å² (23/75 = 30%). Average: 211 atoms, 52 residues per interface.
- **Conditions:** 75-complex dataset
- **Confidence:** High
- **Benchmarkable:** Framework — expected range for SPINK7-KLK5

### Feature 8: 1tgs (Trypsinogen-PSTI) Structural Data
- **Type:** Structural
- **Value:** BSA = 1730 Å², 10 H-bonds, 18 waters, 1.8 Å resolution
- **Conditions:** Crystal structure analysis
- **Confidence:** High
- **Benchmarkable:** Yes — SPINK1-trypsinogen is closest analog to SPINK7-KLK5

### Feature 9: 2kai (Kallikrein-PTI) Structural Data
- **Type:** Structural
- **Value:** BSA = 1440 Å², 10 H-bonds, 2.5 Å resolution
- **Conditions:** Crystal structure analysis
- **Confidence:** High
- **Benchmarkable:** Yes — kallikrein-family complex directly relevant to KLK5

## PDB Codes Referenced
75 complexes including: **1tgs**, **2kai**, **2ptc**, 1avw, 1mct, 3tpi, **1cho**, **1acb**, 1cbw, **1ppf**, **1fle**, 1hia, **3sgb**, 1mkw, **1cse**, **2sic**, **2sni**, 1stf, 4cpa, 1bth, 4htc, 1tbq, 1toc, 1dan, 1vfb, 1mlc, 1jhl, 3hfl, 3hfm, 1fbi, 1mel, 1jel, 1nsn, 1osp, 1nca, 1nmb, 1dvf, 1iai, 1nfd, 1kb5, 1ao7, **1brs**, 1dfj, 1dhk, 1fss, 1gla, 1udi, 1ydr, 2pcc, 1tx4, 1gua, 1a2k, 1efu, 1aip, 1gg2, 1got, 2trc, 1agr, 1fin, 1a0b, 1fc2, 1igc, 1ak4, 1efn, 1atn, 2btf, 1dkg, 1ebp, 1hwg, 1seb, 1tco, 1ycs

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** Critical reference for SPINK7-KLK5 benchmarking. Provides quantitative structural targets for MD validation: protease-inhibitor interfaces should have ~1530 Å² BSA, ~10 H-bonds, ~15 interface waters, 61% non-polar character, and V/Vo ≈ 1.00. The inclusion of 1tgs (SPINK1-trypsinogen) and 2kai (kallikrein-PTI) makes this directly applicable to the SPINK7-KLK5 system. These structural parameters are the primary criteria for evaluating MD equilibration quality.
- **Key Limitation:** Structural analysis only — no binding free energies, kinetics, or thermodynamic data beyond interface geometry.
