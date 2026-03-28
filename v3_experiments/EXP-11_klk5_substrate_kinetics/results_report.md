# EXP-11: KLK5 Substrate Kinetics — Results Report

## Abstract

This report documents the outcome of EXP-11, which aimed to compute the catalytic efficiency (kcat/KM) for KLK5 substrate hydrolysis. The experiment is classified as **INCONCLUSIVE** due to two independent blocking factors: (1) no KLK5-substrate peptide complex structure exists in the PDB, and (2) the pipeline lacks the QM/MM capabilities required to model covalent bond breaking during peptide hydrolysis. The AMBER14 force field does not parameterize scissile bond cleavage.

## Introduction / Background

KLK5 is a trypsin-like serine protease that cleaves peptide bonds after arginine and lysine residues. The catalytic efficiency (kcat/KM) is a fundamental kinetic parameter that describes the enzyme's proficiency at converting substrate to product. Debela et al. (2006) reported kcat/KM = 2.5 × 10⁵ M⁻¹ s⁻¹ for KLK5 hydrolysis of synthetic peptide substrates.

Computing kcat from molecular simulations requires modeling the full catalytic cycle, including formation of the Michaelis complex, nucleophilic attack by the catalytic serine, formation of the tetrahedral intermediate, and collapse to the acyl-enzyme intermediate. These steps involve covalent bond making and breaking that cannot be described by classical molecular mechanics force fields.

## Hypothesis

**H-11:** The pipeline can compute kcat/KM for KLK5 substrate hydrolysis within one order of magnitude of the experimental benchmark (2.5 × 10⁵ M⁻¹ s⁻¹).

## Methods

**Planned approach:** Build KLK5-substrate Michaelis complex, perform enhanced sampling simulations along the reaction coordinate to compute activation free energies, and use transition state theory to estimate kcat.

**Pipeline components:** OpenMM 8.1, AMBER14/SB force field, TIP3P water model.

**Status:** Methods could not be executed due to (1) missing KLK5-substrate complex structure and (2) absence of QM/MM methodology in the pipeline.

## Controls

Not applicable — experiment was not executed.

## Results

**Status: INCONCLUSIVE — experiment could not be executed.**

This experiment is blocked by two independent factors:

1. **Missing structural data:** No KLK5-substrate peptide complex structure exists in the PDB. The KLK5 apoenzyme is available (PDB 2PSX), but computing substrate kinetics requires a Michaelis complex with the substrate peptide correctly positioned in the active site. The pipeline's §21 (Pipeline Immutability) prohibits de novo docking of substrates into the active site.

2. **Methodology limitation:** Even if a Michaelis complex structure were available, the pipeline's force field (AMBER14/SB) does not parameterize scissile bond breaking. Computing kcat requires modeling the catalytic mechanism, which involves covalent bond formation and cleavage. This requires QM/MM or reactive force field methods that are outside the current pipeline's scope.

**Feature F-11 (kcat/KM for KLK5 substrate hydrolysis):** Not measured.

**Benchmark comparison:** Not possible.

## Discussion

The inability to compute enzyme kinetics from classical MD simulations is a well-known limitation of force-field-based approaches. Classical force fields treat bonds as harmonic springs with fixed connectivity — they cannot model bond formation or cleavage. Computing kcat for serine protease catalysis requires:

1. **QM/MM methods:** Hybrid quantum mechanics/molecular mechanics approaches that treat the active site quantum mechanically while describing the surrounding protein classically. QM/MM free energy perturbation or umbrella sampling along the reaction coordinate can yield activation barriers from which kcat is estimated via transition state theory.

2. **Reactive force fields:** Empirical reactive force fields (e.g., ReaxFF) that allow bond breaking and formation, though these are less accurate for enzyme catalysis than QM/MM approaches.

3. **Empirical valence bond (EVB):** Warshel's EVB method provides an efficient framework for computing activation free energies in enzymes using classical sampling with empirically calibrated diabatic states.

None of these methods are implemented in the current pipeline, which is designed for equilibrium MD, SMD pulling, and umbrella sampling with classical force fields. This is an inherent scope limitation rather than a pipeline deficiency — the pipeline was designed for binding thermodynamics, not enzyme kinetics.

The KM component could in principle be estimated from binding free energy calculations if a Michaelis complex structure were available, since KM ≈ Kd for simple Michaelis-Menten kinetics. However, kcat estimation fundamentally requires reaction modeling capabilities.

## Conclusions

EXP-11 is INCONCLUSIVE. KLK5 substrate kinetics cannot be computed because (a) no KLK5-substrate co-crystal structure exists in the PDB and the pipeline prohibits de novo substrate docking, and (b) the pipeline lacks QM/MM capabilities required to model covalent bond breaking during peptide hydrolysis. This experiment requires both structural data (a KLK5-substrate Michaelis complex) and methodological extensions (QM/MM or equivalent reactive methods) that are outside the current pipeline scope.

## Figures

No figures generated — experiment not executed due to missing structural data.

## References

1. Debela, M., et al. (2006). Specificity profiling of seven human tissue kallikreins reveals individual subsite preferences. *Journal of Biological Chemistry*, 281(35), 25678–25688.
2. Warshel, A. (1991). *Computer Modeling of Chemical Reactions in Enzymes and Solutions*. Wiley-Interscience.
3. Senn, H. M., & Thiel, W. (2009). QM/MM methods for biomolecular systems. *Angewandte Chemie International Edition*, 48(7), 1198–1229.
4. PDB Entry 2PSX: Crystal structure of KLK5.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp
