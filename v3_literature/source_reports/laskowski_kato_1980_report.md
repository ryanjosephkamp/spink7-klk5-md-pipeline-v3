# Source Report: S-12

## Bibliographic Information
- **ID:** S-12
- **Authors:** Laskowski M Jr, Kato I
- **Title:** Protein Inhibitors of Proteinases
- **Journal:** Annual Review of Biochemistry
- **Year:** 1980
- **Volume/Pages:** 49:593–626
- **DOI:** N/A (pre-DOI era)
- **PDF Filename:** S-12.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category IV: Kazal-Type Inhibitor Biology

## Summary
This is the foundational review article that defines the **Standard Mechanism** of serine proteinase inhibition: E + I ⇌ L ⇌ C ⇌ X ⇌ L* ⇌ E + I* (loose complex → tight complex → modified inhibitor complex → release of modified inhibitor). It establishes the thermodynamic and kinetic framework used by all subsequent work on canonical protease inhibitors, including the Kazal family to which SPINK5 and SPINK7 belong. The paper classifies 8+ families of serine proteinase inhibitors and defines key parameters: association equilibrium constants (Ka), hydrolysis equilibrium constants (K_hyd), and catalytic efficiency (k_cat/K_m).

## Extracted Features

### Feature 1: Standard Mechanism Kinetic Framework
- **Type:** Kinetic/Thermodynamic
- **Value:** E + I ⇌ L ⇌ C ⇌ X ⇌ L* ⇌ E + I*; Ka = 10⁷–10¹³ M⁻¹; k_cat/K_m = 10⁴–10⁶ M⁻¹s⁻¹ for reactive-site hydrolysis
- **Conditions:** General framework applicable to all canonical inhibitors
- **Confidence:** High — universally accepted
- **Benchmarkable:** Framework — defines how Ki/Ka values relate to MD-computed binding free energies

### Feature 2: Hydrolysis Equilibrium
- **Type:** Thermodynamic
- **Value:** K_hyd ≈ 1 at neutral pH (virgin and modified inhibitors near equimolar at equilibrium)
- **Conditions:** Neutral pH
- **Confidence:** High
- **Benchmarkable:** No — relates to covalent bond hydrolysis, not binding

### Feature 3: Catalytic Contact Distance
- **Type:** Structural
- **Value:** P1 carbonyl carbon to Oγ of catalytic Ser195 = 2.6 Å
- **Conditions:** Crystal structure measurements
- **Confidence:** High
- **Benchmarkable:** Yes — directly measurable from MD trajectories

### Feature 4: Kazal Family Classification
- **Type:** Classification
- **Value:** Kazal family includes pancreatic secretory trypsin inhibitors, ovomucoids (3 tandem domains), ovoinhibitors (6 tandem domains). Table 3 lists alternative P3–P3' residues for Kazal family.
- **Conditions:** N/A
- **Confidence:** High
- **Benchmarkable:** No — classification data

### Feature 5: Inhibitor Family Diversity
- **Type:** Classification
- **Value:** 8+ families: Kazal, Kunitz (BPTI), Bowman-Birk, STI, SSI, potato inhibitor I/II, squash seed, serpins
- **Conditions:** Table 4
- **Confidence:** High
- **Benchmarkable:** No — classification data

## PDB Codes Referenced
None (published 1980, pre-PDB era)

## Usefulness Assessment
- **Overall Rating:** Medium
- **Rationale:** Essential theoretical framework for interpreting all protease-inhibitor binding data. Defines the standard mechanism, Ka ranges, and K_hyd concepts. The catalytic contact distance of 2.6 Å is a direct structural benchmark. However, no directly usable PDB codes or system-specific quantitative data for SPINK-KLK interactions. This is a foundational reference, not a data source.
- **Key Limitation:** Pre-dates KLK and SPINK nomenclature entirely. No structural coordinates or PDB codes available.
