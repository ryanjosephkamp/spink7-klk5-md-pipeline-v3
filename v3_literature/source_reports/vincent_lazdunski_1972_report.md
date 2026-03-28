# Source Report: S-16

## Bibliographic Information
- **ID:** S-16
- **Authors:** Vincent JP, Lazdunski M
- **Title:** Trypsin-Pancreatic Trypsin Inhibitor Association. Dynamics of the Interaction and Role of Disulfide Bridges
- **Journal:** Biochemistry
- **Year:** 1972
- **Volume/Pages:** 11(16):2967–2977
- **DOI:** N/A (ACS, 1972)
- **PDF Filename:** S-16.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category V: BPTI-Trypsin Benchmark

## Summary
Foundational kinetic and thermodynamic characterization of the BPTI-trypsin interaction — the gold-standard protease-inhibitor system. Establishes K_d = 6 × 10⁻¹⁴ M (ΔG ≈ −18 kcal/mol) with complete kon/koff kinetics. Systematically probes the role of each of BPTI's three disulfide bridges by selective reduction (Cys14-Cys38), showing that removal of the reactive-site disulfide increases K_d by ~10⁵-fold. Also characterizes pseudotrypsin and modified trypsin complexes. Provides activation energy for association (10.5 kcal/mol) and pH dependence (k1 maximal pH 8–10, controlled by His46 pK 7.05).

## Extracted Features

### Feature 1: BPTI-Trypsin Dissociation Constant
- **Type:** Thermodynamic
- **Value:** K_d = 6 × 10⁻¹⁴ M (ΔG ≈ −18 kcal/mol at 25°C)
- **Conditions:** 50 mM Tris, 50 mM CaCl₂, 0.1 M NaCl, pH 8.0, 25°C
- **Confidence:** High — universally cited benchmark
- **Benchmarkable:** Yes — gold-standard binding free energy target

### Feature 2: Association/Dissociation Rate Constants
- **Type:** Kinetic
- **Value:** k_on = 1.1 × 10⁶ M⁻¹s⁻¹; k_off = 6.6 × 10⁻⁸ s⁻¹ (complex half-life ~17 weeks)
- **Conditions:** pH 8.0, 25°C
- **Confidence:** High
- **Benchmarkable:** Yes — direct targets for kinetic MD benchmarking

### Feature 3: Disulfide Bridge Variant Kinetics (Table I)
- **Type:** Kinetic/Thermodynamic
- **Value:** R*PTI (reduced Cys14-Cys38): K_d = 1.8 × 10⁻⁹ M, k_on = 3.2 × 10⁵, k_off = 5.7 × 10⁻⁴; RCAM*PTI: K_d = 1.7 × 10⁻¹⁰ M; RAE*PTI: K_d = 9.1 × 10⁻⁹ M; RCOM*PTI: no association
- **Conditions:** pH 8.0, 25°C
- **Confidence:** High
- **Benchmarkable:** Yes — ΔΔG values from disulfide modification directly comparable to FEP

### Feature 4: Activation Energy for Association
- **Type:** Kinetic
- **Value:** E_a = 10.5 kcal/mol (Arrhenius plot, 5–45°C)
- **Conditions:** pH 8.0
- **Confidence:** High
- **Benchmarkable:** Partially — can be compared to PMF barrier heights

### Feature 5: pH Dependence of Association
- **Type:** Kinetic
- **Value:** k_on maximal pH 8–10; controlled by His46 ionization (pK 7.05)
- **Conditions:** Various pH values, 25°C
- **Confidence:** High
- **Benchmarkable:** Partially — constant-pH MD could probe this

### Feature 6: Modified Enzyme Kinetics
- **Type:** Kinetic
- **Value:** Pseudotrypsin-PTI: K_d = 9.0 × 10⁻⁹ M, k_on = 7.0 × 10⁴; RCOM*trypsin-PTI: K_d = 6.0 × 10⁻⁹ M, k_on = 2 × 10⁴
- **Conditions:** pH 8.0, 25°C
- **Confidence:** High
- **Benchmarkable:** Yes — modified enzyme binding provides additional benchmarks

### Feature 7: BPTI-Chymotrypsin Kd
- **Type:** Thermodynamic
- **Value:** K_d = 9 × 10⁻⁹ M
- **Conditions:** pH 8.0, 25°C
- **Confidence:** High
- **Benchmarkable:** Yes — cross-system specificity benchmark

## PDB Codes Referenced
None (1972 pre-dates PDB). References crystallographic data from Huber et al. 1970, 1971.

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** Provides the gold-standard K_d for the BPTI-trypsin system (6 × 10⁻¹⁴ M) with complete kinetic parameters. The disulfide bridge variant data enables systematic structure-activity relationship analysis directly relevant to Kazal-type inhibitors (which also use disulfide-stabilized binding loops). The activation energy and pH dependence provide additional benchmarking targets. This is the primary reference for any computational study of canonical protease inhibitors.
- **Key Limitation:** No structural coordinates (pre-PDB). BPTI is a Kunitz-type inhibitor, not Kazal-type; direct applicability to SPINK7-KLK5 requires analogy arguments.
