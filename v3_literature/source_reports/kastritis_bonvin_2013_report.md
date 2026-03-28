# Source Report: S-23

## Bibliographic Information
- **ID:** S-23
- **Authors:** Kastritis PL, Bonvin AMJJ
- **Title:** On the Binding Affinity of Macromolecular Interactions: Daring to Ask Why Proteins Interact
- **Journal:** Journal of the Royal Society Interface
- **Year:** 2013
- **Volume/Pages:** 10:20120835
- **DOI:** 10.1098/rsif.2012.0835
- **PDF Filename:** S-23.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category VII: PPI Energetics

## Summary
Comprehensive review of macromolecular binding affinity: theory, measurement, and prediction. Establishes empirical relationships between buried surface area (BSA) and binding free energy (ΔG_bond ≈ 0.025 × BSA kcal/mol per Å²). Discusses error sources in Kd measurements (20–50% → 0.1–0.25 kcal/mol ΔG error), temperature effects (2× per 17°C), pH effects (10× over 3 pH units), and the diffusion-limited kon ceiling (~10¹⁰ M⁻¹s⁻¹). Reviews lock-and-key, induced fit, and conformational selection models. 346 references.

## Extracted Features

### Feature 1: BSA-Affinity Relationship
- **Type:** Thermodynamic
- **Value:** ΔG_bond ≈ 0.025 × BSA (kcal/mol per Å²); correlation r = 0.54 for 70 rigid-body complexes
- **Conditions:** Compiled dataset
- **Confidence:** Moderate (r = 0.54 is a weak correlation)
- **Benchmarkable:** Yes — quick sanity check for MD-computed interface area

### Feature 2: Experimental Error Bounds
- **Type:** Methodological
- **Value:** Kd measurement errors: 20–50% → ΔG errors of 0.1–0.25 kcal/mol; temperature variation 18–35°C → 2× Kd change; pH variation 5.5–8.5 → 10× Kd change
- **Conditions:** Literature compilation
- **Confidence:** High
- **Benchmarkable:** Framework — defines expected experimental uncertainty against which to compare MD accuracy

### Feature 3: kon Limits and TransComp Model
- **Type:** Kinetic
- **Value:** Diffusion collision limit ~10¹⁰ M⁻¹s⁻¹; TransComp prediction model r² = 0.72
- **Conditions:** Theory + empirical fitting
- **Confidence:** High
- **Benchmarkable:** Context — upper bound on association rates

### Feature 4: Hot Spot Definition
- **Type:** Thermodynamic
- **Value:** Hot spot: >4 kcal/mol destabilization on Ala mutation; warm spot: 1–2 kcal/mol
- **Conditions:** Standard definition from alanine scanning
- **Confidence:** High (community consensus)
- **Benchmarkable:** Framework — defines thresholds for computational alanine scanning

### Feature 5: van der Waals Energy–koff Correlation
- **Type:** Thermodynamic/Kinetic
- **Value:** vdW energy vs koff: r = −0.45 (N=53, p<0.001); desolvation energy vs koff: r = −0.28 (N=53, p=0.038); rigid binders: vdW vs koff r = −0.56 (N=26, p=0.0024)
- **Conditions:** Compiled dataset, scoring function analysis
- **Confidence:** Moderate
- **Benchmarkable:** Yes — MD-computed vdW energies should correlate with experimental koff

## PDB Codes Referenced
1UBQ, 1UBI, 1F6M (and many others in figures/references)

## Usefulness Assessment
- **Overall Rating:** Medium
- **Rationale:** Essential reference for the theoretical framework of PPI binding energetics. Provides: (1) empirical BSA-ΔG relationship for sanity-checking MD results; (2) experimental error bounds for calibrating expected accuracy; (3) hot/warm spot definitions for computational alanine scanning; (4) scoring function–kinetics correlations. However, this is a review paper with no original data specific to protease-inhibitor systems.
- **Key Limitation:** Review/theory paper — no original experimental data. General PPI framework, not system-specific.
