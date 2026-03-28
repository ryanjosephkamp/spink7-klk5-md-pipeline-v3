# Source Report: S-15

## Bibliographic Information
- **ID:** S-15
- **Authors:** Rawlings ND, Tolle DP, Barrett AJ
- **Title:** Evolutionary Families of Peptidase Inhibitors
- **Journal:** Biochemical Journal
- **Year:** 2004
- **Volume/Pages:** 378:705–716
- **DOI:** N/A
- **PDF Filename:** S-15.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category IV: Kazal-Type Inhibitor Biology

## Summary
MEROPS database-based taxonomic classification of peptidase inhibitors into 48 families and 26 clans based on protein fold. Establishes the nomenclature system used throughout the field. SPINK5 (LEKTI) is classified as family I1 (Kazal), identifier LI01-004, containing 15 Kazal domains — the largest known compound inhibitor in the Kazal family. Provides PDB fold-type examples for each clan. Mentions Netherton syndrome as a genetic disease resulting from SPINK5 mutations.

## Extracted Features

### Feature 1: SPINK5/LEKTI Classification
- **Type:** Classification
- **Value:** Family I1 (Kazal), identifier LI01-004, 15 Kazal domains; targets family S1 serine peptidases
- **Conditions:** MEROPS database classification
- **Confidence:** High — authoritative database
- **Benchmarkable:** No — nomenclature/classification

### Feature 2: Kazal Family Type Example
- **Type:** Structural
- **Value:** Type example: ovomucoid unit 3 (Meleagris gallopavo), SWISS-PROT P01004 (residues 135–185), Pfam PF00050
- **Conditions:** Family I1 definition
- **Confidence:** High
- **Benchmarkable:** Structural reference — connects to PDB 1OVO

### Feature 3: Clan IA Fold Type PDB
- **Type:** Structural
- **Value:** PDB 1OVO (clan IA, Kazal fold)
- **Conditions:** Clan fold type-example
- **Confidence:** High
- **Benchmarkable:** Yes — structural reference for Kazal fold

### Feature 4: Inhibitor Family Diversity
- **Type:** Classification
- **Value:** 48 families, 26 clans, 2500+ sequences; family I4 (serpins) largest with >500 sequences; I1 (Kazal), I2 (Kunitz), I25 (cystatin) each >200 sequences; 19 families use standard mechanism
- **Conditions:** MEROPS database (2004)
- **Confidence:** High
- **Benchmarkable:** No — population statistics

### Feature 5: Compound Inhibitor Architectures
- **Type:** Structural
- **Value:** Table 2 lists compound inhibitors with 1–15 inhibitor units. SPINK5/LEKTI has 15 Kazal domains (largest).
- **Conditions:** Domain architecture analysis
- **Confidence:** High
- **Benchmarkable:** No — domain count, not binding data

### Feature 6: Clan Fold-Type PDB Codes
- **Type:** Structural
- **Value:** 24 clan fold-type examples with PDB codes: 1OVO (IA/Kazal), 4PTI (IB/Kunitz), 1AVU (IC/STI), 1ATU (ID/squash), 1F2S (IE), 1PI2 (IF/potato I), 1EGP (IG), 1CEW (IH/cystatin), 1B1U (IJ), 1SMP (IK), 4HTC (IM), 1ECZ (IN/ecotin), 1FLE (IP/chelonianin), 1I3S (IQ/SFTI), 1F32 (IR), 1DTD (IS/Ascaris), 1UEA (IT/antistasin), 1BHU (IU), 1E31 (IV), 1GL1 (IW/grasshopper), 1ICF (IX), 3SSI (IY/SSI), 1AVG (IZ), 1DPJ (JA)
- **Conditions:** Structural classification
- **Confidence:** High
- **Benchmarkable:** 1OVO is relevant as Kazal clan representative

## PDB Codes Referenced
1OVO, 4PTI, 1AVU, 1ATU, 1F2S, 1PI2, 1EGP, 1CEW, 1B1U, 1SMP, 4HTC, 1ECZ, 1FLE, 1I3S, 1F32, 1DTD, 1UEA, 1BHU, 1E31, 1GL1, 1ICF, 3SSI, 1AVG, 1DPJ

## Usefulness Assessment
- **Overall Rating:** Medium
- **Rationale:** Essential reference for nomenclature and classification. Confirms SPINK5/LEKTI as a 15-domain Kazal family compound inhibitor (the largest known). The MEROPS classification and Netherton syndrome connection provide important biological framing. PDB 1OVO serves as a structural reference for the Kazal fold. However, no quantitative binding data, kinetics, or system-specific benchmarking data are provided.
- **Key Limitation:** Taxonomy/classification paper — no experimental binding data. PDB codes are fold-type examples rather than benchmark systems with associated thermodynamic data.
