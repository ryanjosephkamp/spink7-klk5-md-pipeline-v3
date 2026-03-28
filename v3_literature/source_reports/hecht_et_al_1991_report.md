# Source Report: S-19

## Bibliographic Information
- **ID:** S-19
- **Authors:** Hecht HJ, Szardenings M, Collins J, Schomburg D
- **Title:** Three-dimensional Structure of the Complexes between Bovine Chymotrypsinogen A and Two Recombinant Variants of Human Pancreatic Secretory Trypsin Inhibitor (Kazal-type)
- **Journal:** Journal of Molecular Biology
- **Year:** 1991
- **Volume/Pages:** 220:711–722
- **DOI:** 0022-2836/91/150711-12
- **PDF Filename:** S-19.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VI: SPINK1-Trypsin

## Summary
Crystal structures of two recombinant human PSTI (SPINK1) variants (PSTI3 and PSTI4) in complex with bovine chymotrypsinogen A at 2.3 Å resolution. This is the most directly relevant structural reference for the SPINK7-KLK5 system, as PSTI/SPINK1 belongs to the same Kazal family as SPINK7. Provides Ki values for PSTI variants (Ki = 1.6–2.4 × 10⁻¹¹ M for chymotrypsin), detailed interface contacts and H-bond distances, binding loop structural conservation (r.m.s. 0.25–0.35 Å across Kazal inhibitors), and evidence of inhibitor-induced zymogen refolding in the chymotrypsinogen active site.

## Extracted Features

### Feature 1: PSTI3/4-Chymotrypsin Ki Values
- **Type:** Thermodynamic
- **Value:** PSTI3-chymotrypsin Ki = 1.6 × 10⁻¹¹ M; PSTI4-chymotrypsin Ki = 2.4 × 10⁻¹¹ M
- **Conditions:** From Szardenings 1989 thesis
- **Confidence:** High
- **Benchmarkable:** Yes — ΔG ≈ −14.5 kcal/mol, directly comparable to FEP

### Feature 2: PSTI4-Elastase Ki
- **Type:** Thermodynamic
- **Value:** PSTI4-elastase (HLE) Ki = 3.7 × 10⁻¹¹ M
- **Conditions:** From Szardenings 1989
- **Confidence:** High
- **Benchmarkable:** Yes — cross-enzyme selectivity benchmark

### Feature 3: Crystallographic Resolution and Quality
- **Type:** Structural
- **Value:** 2.3 Å resolution; R = 19.5%; space group P4₁2₁2; a = 84.4, c = 86.7 Å; 50–56 water molecules; mean B-factor main chain 13.0–13.2 Å²
- **Conditions:** X-ray crystallography
- **Confidence:** High
- **Benchmarkable:** Yes — starting structure for MD

### Feature 4: Binding Loop Structural Conservation
- **Type:** Structural
- **Value:** PSTI3 vs PSTI4 backbone r.m.s. = 0.25 Å (I4–I56); PSTI3/4 vs porcine PSTI core (I20–I42) = 0.35 Å; binding loop (I15–I21) = 0.31 Å; comparison to other inhibitor complexes (I16–I20) ≈ 0.3 Å
- **Conditions:** Crystal structure superposition
- **Confidence:** High
- **Benchmarkable:** Yes — validates template-based modeling of SPINK7

### Feature 5: Catalytic Contact Geometry
- **Type:** Structural
- **Value:** Ser195 Oγ to PSTI P1 carbonyl carbon = 2.6 Å; to N = 3.0 Å; to O = 2.9 Å
- **Conditions:** Crystal structure (Table 7)
- **Confidence:** High
- **Benchmarkable:** Yes — key geometric validation target for MD

### Feature 6: Interface Hydrogen Bonds
- **Type:** Structural
- **Value:** Key H-bonds: Cys I16 N–Gly A216 O (2.9 Å), Leu I18 O–Gly A193 N (2.4 Å), Arg I21 Nη²–Cys A58 O (2.3 Å)
- **Conditions:** Crystal structure contacts
- **Confidence:** High
- **Benchmarkable:** Yes — H-bond distances directly measurable from MD

### Feature 7: Mutation Effects on Binding
- **Type:** Kinetic
- **Value:** Arg I21→Phe: 150-fold decrease; Glu I40→Ala: 25-fold decrease; Val I36→Asp in PSTI5: 10-fold increase; Val I36→Asp in PSTI4: 1000-fold decrease; P5/P6→Glu: 2-fold decrease; Leu I18→Tyr: 1.5-fold increase
- **Conditions:** Mutagenesis data from Szardenings 1989
- **Confidence:** Moderate — from thesis data
- **Benchmarkable:** Yes — per-residue ΔΔG for FEP validation

### Feature 8: Inhibitor-Induced Zymogen Refolding
- **Type:** Structural
- **Value:** Chymotrypsinogen active-site refolding upon PSTI binding: r.m.s. deviation from active chymotrypsin (5CHA) = 0.26 Å (A186–A195 + A217–A222)
- **Conditions:** Comparison: free chymotrypsinogen (2CGA) vs complexed
- **Confidence:** High
- **Benchmarkable:** Partially — represents induced-fit effects testable by MD

### Feature 9: PSTI Structural Features
- **Type:** Structural
- **Value:** Human PSTI: 56 residues, 3 disulfide bridges; 73% sequence homology with porcine PSTI
- **Conditions:** Sequence/structure analysis
- **Confidence:** High
- **Benchmarkable:** Structural context for SPINK7 homology model

## PDB Codes Referenced
- **1TGS** (porcine PSTI-trypsinogen complex — primary reference)
- **1CHG** (free chymotrypsinogen — molecular replacement model)
- **1CHO** (chymotrypsin-ovomucoid complex)
- **2CGA** (free chymotrypsinogen)
- **5CHA** (α-chymotrypsin)
- **4CHA** (α-chymotrypsin)
- **3SGB** (ovomucoid 3d–proteinase B)
- **2SNI** (CI-2–subtilisin novo)
- **2SEC** / **1CSE** (eglin c–subtilisin Carlsberg)
- **1OVO** (ovomucoid 3d free)
- **2CI2** (CI-2 free)
- PSTI3/4 structures: deposition was in progress at publication

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** The most directly relevant structural reference for SPINK7-KLK5 modeling. PSTI/SPINK1 is the same Kazal family as SPINK7. Provides: (1) a Kazal inhibitor-protease complex crystal structure at 2.3 Å for template-based modeling; (2) tight Ki values (10⁻¹¹ M) with mutation effects; (3) detailed interface geometry (H-bonds, catalytic contact) for MD validation; (4) evidence that the canonical binding loop is structurally conserved across Kazal inhibitors (r.m.s. ~0.3 Å); (5) induced-fit effects relevant to KLK5 activation state.
- **Key Limitation:** Complex is with chymotrypsinogen (a zymogen), not an active serine protease like KLK5. PSTI3/4 are engineered variants, not wild-type SPINK1. PDB codes for the actual PSTI3/4 structures may not be deposited.
