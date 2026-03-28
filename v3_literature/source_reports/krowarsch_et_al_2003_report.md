# Source Report: S-14

## Bibliographic Information
- **ID:** S-14
- **Authors:** Krowarsch D, Cierpicki T, Jelen F, Otlewski J
- **Title:** Canonical Protein Inhibitors of Serine Proteases
- **Journal:** Cellular and Molecular Life Sciences (CMLS)
- **Year:** 2003
- **Volume/Pages:** 60:2427–2444
- **DOI:** 10.1007/s00018-003-3120-x
- **PDF Filename:** S-14.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category IV: Kazal-Type Inhibitor Biology

## Summary
The most comprehensive single reference for canonical protease inhibitor structural and energetic data. Provides explicit PDB codes for free inhibitors and protease-inhibitor complexes across all major families (Table 2), quantitative binding interface metrics (1400 Å² buried area, >100 van der Waals contacts, 8–15 hydrogen bonds, 10–18 inhibitor + 17–30 protease residues), energetic decomposition of binding (P1 residue = 70% of association energy, ~50% of interface contact area), and the Laskowski additivity framework for predicting mutation effects. Contains critical benchmark PDB codes for Kazal family: 2ovo (OMSVP3 free, 1.5 Å) and 1cho (OMTKY3:chymotrypsin complex, 1.8 Å).

## Extracted Features

### Feature 1: PDB Codes for Kazal Family Structures
- **Type:** Structural
- **Value:** Free inhibitor: 2ovo (OMSVP3, 1.5 Å resolution); Complex: 1cho (OMTKY3:chymotrypsin, 1.8 Å)
- **Conditions:** X-ray crystallography
- **Confidence:** High — established structures
- **Benchmarkable:** Yes — starting structures for Kazal-protease MD benchmarking

### Feature 2: Comprehensive PDB Code Table (Table 2)
- **Type:** Structural
- **Value:** 16 families × (free + complex) PDB codes. Key entries: BPTI 5pti/1f7z (1.0/1.5 Å), Kazal 2ovo/1cho, Potato I 2ci2/1cse (eglin c:subtilisin, 1.2 Å), Squash 1lu0/2btc, Ecotin 1ifg/1aaz, STI 1avu/1avw, BBI 1c2a/1df9, SFTI-1 1jbl/1sfi, Antistasin 1bx7/1eja, Ascaris 1ccv/1eai, Grasshopper 1pmc/1gl1, SSI 3ssi/2sic, Potato II 1tih/4sgb, Cereal 1bea, Chelonianin 2rel/1fle, Rapeseed 1jxc
- **Conditions:** X-ray and NMR structures
- **Confidence:** High
- **Benchmarkable:** Yes — comprehensive starting structure database for multi-system benchmarking

### Feature 3: Binding Interface Metrics
- **Type:** Structural
- **Value:** 10–18 inhibitor residues + 17–30 protease residues at interface; >100 van der Waals contacts; 8–15 hydrogen bonds; ~1400 Å² buried surface area
- **Conditions:** Compiled across canonical complexes
- **Confidence:** High — averaged from multiple crystal structures
- **Benchmarkable:** Yes — all directly measurable from MD trajectories

### Feature 4: P1 Residue Dominance in Binding
- **Type:** Thermodynamic
- **Value:** P1 residue contributes up to 70% of association energy and ~50% of interface contact area
- **Conditions:** Mutagenesis and structural analysis
- **Confidence:** High
- **Benchmarkable:** Yes — can be tested via computational alanine scanning or FEP

### Feature 5: Association Rate and Equilibrium Constants
- **Type:** Kinetic
- **Value:** k_on ≈ 10⁶ M⁻¹s⁻¹; Ka up to 10¹⁷ M⁻¹ achievable by engineering; Ki values typically 10⁶–10⁹-fold lower than Km
- **Conditions:** General for canonical inhibitors
- **Confidence:** High
- **Benchmarkable:** Ka/Ki ranges provide context for computed values

### Feature 6: Catalytic Contact Distance
- **Type:** Structural
- **Value:** P1 carbonyl carbon to catalytic serine Oγ ≈ 2.7 Å
- **Conditions:** Crystal structures
- **Confidence:** High
- **Benchmarkable:** Yes — direct MD validation target

### Feature 7: Hydrolysis Equilibrium Constants
- **Type:** Thermodynamic
- **Value:** K_hyd = 0.4–35 for natural ovomucoid variants; typical for normal proteins: 10⁻³ to 10⁻⁸
- **Conditions:** Kazal-family ovomucoid reactive-site variants
- **Confidence:** High
- **Benchmarkable:** Context — explains why canonical inhibitors are "temporary" (K_hyd ≈ 1)

### Feature 8: Hydrogen Bond Energetic Contributions
- **Type:** Thermodynamic
- **Value:** Single intermolecular H-bond (P1 amide): ~1.5 kcal/mol; conversion of P2–P1 amide to ester reduces ΔG_association by ~1.5 kcal/mol
- **Conditions:** BPTI and OMTKY3 mutagenesis studies
- **Confidence:** High
- **Benchmarkable:** Yes — directly comparable to FEP-computed ΔΔG values

### Feature 9: Inhibitor Thermal Stability
- **Type:** Thermodynamic
- **Value:** BPTI T_den ≈ 100°C, stable in 6M GdmCl; Kazal T_den up to 90°C; STI T_den 65°C
- **Conditions:** Denaturation experiments
- **Confidence:** High
- **Benchmarkable:** No — folding stability, not binding

### Feature 10: P1 Side Chain pK Shifts
- **Type:** Thermodynamic
- **Value:** Up to 5 pH unit shifts for charged P1 side chains in hydrophobic S1 pockets
- **Conditions:** pH-titration of reactive-site variants
- **Confidence:** Moderate
- **Benchmarkable:** Partially — can be probed by constant-pH MD if available

## PDB Codes Referenced
**Table 2 (explicit):** 5pti, 1f7z, 2ovo, 1cho, 2ci2, 1cse, 1lu0, 2btc, 1ifg, 1aaz, 1avu, 1avw, 1c2a, 1df9, 1jbl, 1sfi, 1bx7, 1eja, 1ccv, 1eai, 1pmc, 1gl1, 3ssi, 2sic, 1tih, 4sgb, 1bea, 2rel, 1fle, 1jxc

**Additional from figures:** 1ppe, 1toc, 1ezx, 1pi2, 2sta, 1ba7, 1bbi, 1bik, 1tbq, 1fyb, 1qh2, 1ecz

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** The single most data-rich paper in this batch for MD benchmarking purposes. Provides: (1) explicit PDB codes for Kazal-family free and complexed structures as starting structures; (2) quantitative interface metrics directly measurable from MD trajectories; (3) energetic decomposition data (P1 dominance, H-bond contributions) for validating FEP calculations; (4) the Laskowski additivity framework for predicting and validating mutation effects computationally. The Kazal entries (2ovo, 1cho) are the closest structural analogs to SPINK7-KLK5 for which crystal structures exist.
- **Key Limitation:** No SPINK7 or KLK5-specific data. All Kazal data is for ovomucoid domains (OMTKY3, OMSVP3). Transfer to SPINK7-KLK5 requires structural homology arguments.
