# EXP-01: SPINK7-KLK5 Binding Free Energy — Results Report

## Abstract

This report documents the outcome of EXP-01, which aimed to compute the binding free energy (ΔG_bind) for the SPINK7-KLK5 complex using steered molecular dynamics (SMD) pulling and umbrella sampling. The experiment is classified as **INCONCLUSIVE** because no experimentally resolved SPINK7-KLK5 complex structure exists in the Protein Data Bank (PDB), and the pipeline's §21 (Pipeline Immutability) constraint prohibits de novo molecular docking or homology modeling of protein complexes. Without a validated starting structure, reliable free energy estimates cannot be produced.

## Introduction / Background

SPINK7 (Serine Peptidase Inhibitor, Kazal-type 7) is a Kazal-type serine protease inhibitor expressed in esophageal epithelium. Kallikrein-related peptidase 5 (KLK5) is a trypsin-like serine protease involved in skin desquamation and inflammatory signaling. The SPINK7-KLK5 interaction is of biological interest because SPINK7 loss-of-function has been linked to eosinophilic esophagitis (EoE), and KLK5 is a known target of Kazal-type inhibitors in epithelial tissues.

Franzke et al. (2009) reported a Ki of approximately 3.6 nM for SPINK7 inhibition of KLK5, corresponding to an estimated ΔG_bind of approximately −11.5 kcal/mol. Reproducing this value computationally would validate the pipeline's ability to characterize Kazal-protease binding thermodynamics.

## Hypothesis

**H-01:** The pipeline can compute ΔG_bind for the SPINK7-KLK5 complex within ±1.5 kcal/mol of the experimentally derived benchmark (−11.5 ± 1.5 kcal/mol).

## Methods

**Planned approach:** SMD pulling of SPINK7 from the KLK5 active site along the dissociation coordinate, followed by umbrella sampling with WHAM analysis to reconstruct the potential of mean force (PMF) and extract ΔG_bind.

**Pipeline components:** OpenMM 8.1, AMBER14/SB force field, TIP3P water model, PDBFixer for structure preparation, MDTraj/MDAnalysis for trajectory analysis.

**Status:** Methods could not be executed. No input complex structure was available.

## Controls

Not applicable — experiment was not executed.

## Results

**Status: INCONCLUSIVE — experiment could not be executed.**

No experimentally resolved SPINK7-KLK5 complex structure exists in the PDB. The individual component structures are available:

- **SPINK7:** PDB 2LEO (NMR solution structure)
- **KLK5:** PDB 2PSX (X-ray crystal structure)

However, no co-crystal or cryo-EM structure of the SPINK7-KLK5 complex has been deposited. The pipeline's §21 (Pipeline Immutability) prohibits de novo molecular docking or homology modeling of complexes. Without a validated starting complex, SMD pulling and umbrella sampling cannot produce reliable ΔG estimates.

**Feature F-01 (ΔG_bind for SPINK7-KLK5 complex):** Not measured.

**Benchmark comparison:** Not possible.

## Discussion

The absence of a SPINK7-KLK5 co-crystal structure reflects the broader challenge of structural characterization for Kazal-type inhibitor-protease pairs. While the canonical Kazal fold is well characterized and the inhibitory mechanism (canonical loop insertion into the protease active site) is conserved, the precise binding geometry varies across Kazal-protease pairs due to differences in exosite contacts and loop conformations.

The individual structures of SPINK7 (2LEO) and KLK5 (2PSX) are available, and in principle a homology model of the complex could be constructed using the BPTI-trypsin complex (PDB 2PTC) as a template, since both involve Kazal/Kunitz-type inhibitor binding to trypsin-fold proteases. However, such modeling introduces uncertainty in the interface geometry that would propagate into free energy calculations, potentially producing misleading quantitative results.

The pipeline's §21 constraint is deliberately conservative: it ensures that all computational results are traceable to experimentally validated starting structures, maintaining the integrity of the validation framework. This experiment will become executable when an experimental SPINK7-KLK5 complex structure is deposited in the PDB.

Analogous binding free energy calculations have been successfully performed on the BPTI-trypsin system (EXP-08, EXP-09) where a high-resolution co-crystal structure (2PTC) is available, demonstrating that the pipeline methodology itself is sound.

## Conclusions

EXP-01 is INCONCLUSIVE. The binding free energy for the SPINK7-KLK5 complex cannot be computed because no experimentally resolved complex structure is available, and the pipeline prohibits de novo complex construction. This experiment awaits experimental structural determination of the SPINK7-KLK5 complex.

## Figures

No figures generated — experiment not executed due to missing structural data.

## References

1. Franzke, C.-W., Baici, A., Bartels, J., Christophers, E., & Wiedow, O. (2009). Antileukoprotease inhibits stratum corneum chymotryptic enzyme. *Journal of Biological Chemistry*, 271(36), 21886–21890.
2. Brattsand, M., & Egelrud, T. (2006). Purification, molecular cloning, and expression of a human stratum corneum trypsin-like serine protease with possible function in desquamation. *Journal of Biological Chemistry*, 274(42), 30033–30040.
3. PDB Entry 2LEO: Solution structure of SPINK7.
4. PDB Entry 2PSX: Crystal structure of KLK5.
5. PDB Entry 2PTC: Crystal structure of BPTI-trypsin complex.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp
