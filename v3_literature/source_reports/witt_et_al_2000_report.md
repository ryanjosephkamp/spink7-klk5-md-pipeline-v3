# Source Report: S-21

## Bibliographic Information
- **ID:** S-21
- **Authors:** Witt H, Luck W, Hennies HC, Classen M, Kage A, Lass U, Landt O, Becker M
- **Title:** Mutations in the Gene Encoding the Serine Protease Inhibitor, Kazal Type 1, Are Associated with Chronic Pancreatitis
- **Journal:** Nature Genetics
- **Year:** 2000
- **Volume/Pages:** 25:213–216
- **DOI:** N/A
- **PDF Filename:** S-21.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VI: SPINK1-Trypsin

## Summary
Genetics/clinical paper establishing that SPINK1 mutations (particularly N34S) are associated with chronic pancreatitis. SPINK1 gene is ~7.5 kb with 4 exons on chromosome 5q32, encoding a 79-amino acid protein (23-aa signal peptide + 56-aa mature form) with reactive site Lys41–Ile42. N34S found in 23% of pediatric pancreatitis patients vs. 0.94% controls. Linkage disequilibrium maximum HRR lod score = 5.4. No biophysical or structural data.

## Extracted Features

### Feature 1: SPINK1 Gene Architecture
- **Type:** Genomic
- **Value:** ~7.5 kb, 4 exons, chromosome 5q32; protein: 79 aa (23 signal + 56 mature)
- **Conditions:** Genomic analysis
- **Confidence:** High
- **Benchmarkable:** No — genomic data

### Feature 2: SPINK1 Reactive Site
- **Type:** Structural
- **Value:** Lys41–Ile42 (positions 18–19 of mature peptide)
- **Conditions:** Sequence analysis
- **Confidence:** High
- **Benchmarkable:** Context — identifies P1 residue for MD

### Feature 3: N34S Mutation Association
- **Type:** Epidemiological
- **Value:** N34S in 23% of patients (18/79), 6 homozygous, 12 heterozygous; controls 0.94% (1/106 chromosomes); HRR lod score = 5.4
- **Conditions:** Pediatric pancreatitis cohort (96 patients)
- **Confidence:** High
- **Benchmarkable:** No — clinical association

### Feature 4: Trypsin Inhibition Capacity
- **Type:** Functional
- **Value:** SPINK1 inhibits up to 20% of total trypsin activity in pancreatic juice
- **Conditions:** In vivo estimate
- **Confidence:** Moderate
- **Benchmarkable:** No — physiological context

## PDB Codes Referenced
None

## Usefulness Assessment
- **Overall Rating:** Low
- **Rationale:** Purely genetics/clinical paper with no biophysical, structural, or binding data. Provides important biological context for the SPINK family (reactive site location, disease association) but no benchmarkable data for MD simulations. The N34S mutation association is relevant context for why SPINK-protease interactions matter clinically.
- **Key Limitation:** No structural, thermodynamic, or kinetic data. No PDB codes.
