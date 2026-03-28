# Source Report: S-11

## Bibliographic Information
- **ID:** S-11
- **Authors:** Egelrud T, Brattsand M, Kreutzmann P, Walden K, Vitzithum U, Marx UC, Forssmann WG, Mägert HJ
- **Title:** hK5 and hK7, Two Serine Proteinases Abundant in Human Skin, Are Inhibited by LEKTI Domain 6
- **Journal:** British Journal of Dermatology
- **Year:** 2005
- **Volume/Pages:** 153:1200–1203
- **DOI:** 10.1111/j.1365-2133.2005.06834.x
- **PDF Filename:** S-11.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category IV: Kazal-Type Inhibitor Biology

## Summary
This short communication demonstrates that LEKTI domain 6 (LD6) is the primary inhibitor of KLK5 (hK5) and KLK7 (hK7) among the LEKTI domains tested. LD6 is a non-classical Kazal domain with only 2 disulfide bonds (vs. 3 in classical Kazal domains). The paper also tests recombinant LEKTI domain 15 (rLD15), a classical 3-disulfide Kazal domain, finding it inhibits plasmin but not KLK5/KLK7. This establishes domain-specific selectivity within the multi-domain LEKTI inhibitor.

## Extracted Features

### Feature 1: LEKTI Domain 6 IC₅₀ for KLK5
- **Type:** Kinetic
- **Value:** IC₅₀ ≈ 60 nmol/L (native hK5)
- **Conditions:** 80 mM Tris-HCl pH 8.0, 80 mM NaCl; hK5 at 2.5 µg/mL (~70 nmol/L); substrate 2.5 mmol/L
- **Confidence:** Approximate; single condition
- **Benchmarkable:** Partially — IC₅₀ is condition-dependent but can approximate Ki when [E] ≈ IC₅₀

### Feature 2: LEKTI Domain 6 IC₅₀ for KLK7
- **Type:** Kinetic
- **Value:** IC₅₀ ≈ 340 nmol/L
- **Conditions:** Same buffer; hK7 at ~90 nmol/L
- **Confidence:** Approximate
- **Benchmarkable:** Partially — selectivity ratio hK5/hK7 ≈ 6× is qualitatively benchmarkable

### Feature 3: LEKTI Domain 15 IC₅₀ for Plasmin
- **Type:** Kinetic
- **Value:** IC₅₀ ≈ 100 nmol/L
- **Conditions:** Same buffer; plasmin at 470 ng/mL (~6 nmol/L)
- **Confidence:** Approximate
- **Benchmarkable:** No — plasmin is not the target system

### Feature 4: Inhibition Specificity Panel
- **Type:** Selectivity
- **Value:** LD6 inhibits trypsin, hK5, hK7; does NOT inhibit chymotrypsin, thrombin, granzyme B, tryptase, urokinase, chymase. rLD15 inhibits trypsin, plasmin, weakly elastase; does NOT inhibit hK5, hK7, chymotrypsin, thrombin, granzyme B, tryptase, urokinase, chymase.
- **Conditions:** Table 1
- **Confidence:** Qualitative (+/−)
- **Benchmarkable:** Qualitatively — selectivity pattern can be checked against computed binding energies

### Feature 5: Non-classical Kazal Domain Architecture
- **Type:** Structural
- **Value:** LD6 has 2 disulfide bonds (non-classical); rLD15 has 3 disulfide bonds (classical Kazal)
- **Conditions:** Biochemical characterization
- **Confidence:** High
- **Benchmarkable:** Structural parameter verifiable in MD

## PDB Codes Referenced
None

## Usefulness Assessment
- **Overall Rating:** Medium
- **Rationale:** Provides IC₅₀ data for LD6-KLK5 interaction but not true Ki/Kd values. The IC₅₀ ≈ 60 nM for LD6-KLK5 is useful as a rough benchmark, and the ~6× selectivity for KLK5 over KLK7 is a qualitative target. However, no structural data or binding free energies are provided. The domain-specific selectivity findings are important for biological context.
- **Key Limitation:** IC₅₀ values are condition-dependent and enzyme concentration (~70 nM) is close to IC₅₀ (~60 nM), making tight-binding corrections necessary for true Ki estimation.
