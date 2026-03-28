# Source Report: S-18

## Bibliographic Information
- **ID:** S-18
- **Authors:** Wlodawer A, Walter J, Huber R, Sjölin L
- **Title:** Structure of Bovine Pancreatic Trypsin Inhibitor: Results of Joint Neutron and X-ray Refinement of Crystal Form II
- **Journal:** Journal of Molecular Biology
- **Year:** 1984
- **Volume/Pages:** 180:301–329
- **DOI:** 0022-2836/84/340301-29
- **PDF Filename:** S-18.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category V: BPTI-Trypsin Benchmark

## Summary
Ultra-high resolution joint X-ray (0.94 Å) and neutron (1.8 Å) refinement of BPTI crystal form II. Provides hydrogen atom positions, methyl rotor orientations, and detailed ordered water structure (63 water molecules, 4 internal). Hydrogen exchange data from neutron diffraction identifies 11 protected amide hydrogens after 3 months in D₂O, providing dynamics benchmarks. Comparison of form I and form II shows 0.40 Å r.m.s. deviation (main chain), establishing real conformational variability.

## Extracted Features

### Feature 1: Ultra-High Resolution Crystal Structure
- **Type:** Structural
- **Value:** X-ray 0.94 Å resolution (R = 0.200); neutron 1.8 Å (R = 0.197); space group P2₁2₁2₁; unit cell a = 74.1, b = 23.4, c = 28.9 Å
- **Conditions:** Crystal form II; 63 ordered waters + 1 phosphate
- **Confidence:** High
- **Benchmarkable:** Yes — highest-resolution BPTI reference structure for MD validation

### Feature 2: Hydrogen Atom Positions
- **Type:** Structural
- **Value:** All hydrogen atoms refined from neutron data. R.m.s. bond length deviation: 0.020 Å (non-H), 0.018 Å (H)
- **Conditions:** Joint X-ray + neutron refinement
- **Confidence:** High — neutron data directly locates hydrogens
- **Benchmarkable:** Yes — hydrogen positions for force field validation

### Feature 3: Hydrogen Exchange Protection Factors
- **Type:** Dynamic
- **Value:** 11 amide hydrogens protected from exchange after 3 months D₂O at pH 8.2; 43 deuteriums exchanged (occupancy 0.94 ± 0.12); maximum exchange rate ≤ 0.6 × 10⁻⁶ min⁻¹ for protected sites
- **Conditions:** D₂O soak, pH 8.2, crystal
- **Confidence:** High — agrees with 2D NMR assignments (Wagner & Wüthrich 1982)
- **Benchmarkable:** Yes — MD-predicted solvent accessibility and dynamics can be compared to H/D exchange data

### Feature 4: Ordered Water Structure
- **Type:** Structural
- **Value:** 63 ordered water molecules; 20/63 conserved between crystal forms (deviation < 1 Å); 4 internal (buried) waters in both forms; average B-factor water = 40.9 Å²
- **Conditions:** X-ray + neutron refinement
- **Confidence:** High
- **Benchmarkable:** Yes — solvation structure directly comparable to MD water density maps

### Feature 5: Crystal Form Conformational Variability
- **Type:** Structural
- **Value:** R.m.s. deviation form I vs form II (main chain) = 0.40 Å; solvent content 38% (II) vs 36% (I)
- **Conditions:** Two independent crystal forms
- **Confidence:** High
- **Benchmarkable:** Yes — establishes expected range of conformational variation

### Feature 6: Secondary Structure Parameters
- **Type:** Structural
- **Value:** α-helix Ser47–Gly56 (6 H-bonds, average φ/ψ = −63°/−41.5°, average H-bond length 2.91 Å); 3₁₀ helix near N-terminus (1.5 turns, 3 H-bonds); antiparallel β-sheet residues 16–36
- **Conditions:** Crystal structure
- **Confidence:** High
- **Benchmarkable:** Yes — secondary structure geometry for MD validation

### Feature 7: B-factor Distribution
- **Type:** Structural/Dynamic
- **Value:** Average B-factor protein = 15.3 Å²; average B-factor solvent = 40.9 Å²
- **Conditions:** Crystal structure
- **Confidence:** High
- **Benchmarkable:** Yes — B-factors correlate with MD root-mean-square fluctuations

## PDB Codes Referenced
- This work produced the structure deposited as **4PTI** (form II); references form I structure **5PTI** (Deisenhofer & Steigemann 1975)

## Usefulness Assessment
- **Overall Rating:** Medium-High
- **Rationale:** Provides the highest-resolution structural reference for any canonical protease inhibitor. Hydrogen positions, ordered water structure, hydrogen exchange protection factors, and B-factor distributions are all directly comparable to MD simulation output. The form I/form II comparison establishes a baseline for expected structural variability. This structure (4PTI/5PTI) is the most widely used starting structure for protease-inhibitor MD simulations.
- **Key Limitation:** Structure of free BPTI only — no complex structure in this paper. Kunitz-type inhibitor, not Kazal-type. Structural data, not binding data.
