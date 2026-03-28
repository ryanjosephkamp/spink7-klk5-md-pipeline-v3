# EXP-02: SPINK7-KLK12 Binding Free Energy — Results Report

## Abstract

This report documents the outcome of EXP-02, which aimed to compute the binding free energy (ΔG_bind) for the SPINK7-KLK12 complex using steered molecular dynamics (SMD) and umbrella sampling. The experiment is classified as **INCONCLUSIVE** because neither a SPINK7-KLK12 complex structure nor a KLK12 crystal structure suitable for complex modeling exists in the PDB. The pipeline's §21 (Pipeline Immutability) constraint prohibits de novo molecular docking or homology modeling.

## Introduction / Background

KLK12 (Kallikrein-related peptidase 12) is a serine protease in the kallikrein family with roles in skin barrier function and tissue remodeling. SPINK7 is hypothesized to inhibit KLK12 based on its broad Kazal-type inhibitor activity against trypsin-fold serine proteases. Franzke et al. reported inhibitory activity of SPINK7 against multiple kallikreins, with an estimated Ki of approximately 250 nM for KLK12, corresponding to ΔG_bind ≈ −9.0 kcal/mol.

KLK12 has substantially less structural coverage in the PDB compared to other kallikrein family members such as KLK5 or KLK7. This limited structural information compounds the difficulty of computational characterization of SPINK7-KLK12 binding.

## Hypothesis

**H-02:** The pipeline can compute ΔG_bind for the SPINK7-KLK12 complex within ±2.0 kcal/mol of the experimentally estimated benchmark (−9.0 ± 2.0 kcal/mol).

## Methods

**Planned approach:** SMD pulling of SPINK7 from the KLK12 active site, followed by umbrella sampling with WHAM analysis to reconstruct the PMF and extract ΔG_bind.

**Pipeline components:** OpenMM 8.1, AMBER14/SB force field, TIP3P water model, PDBFixer, MDTraj/MDAnalysis.

**Status:** Methods could not be executed. No input complex structure or suitable KLK12 monomer structure was available.

## Controls

Not applicable — experiment was not executed.

## Results

**Status: INCONCLUSIVE — experiment could not be executed.**

Neither a SPINK7-KLK12 complex structure nor a KLK12 crystal structure suitable for complex modeling exists in the PDB. The same §21 immutability constraint that affects EXP-01 applies here, with the additional complication that KLK12 itself has limited structural coverage.

This experiment requires both:
1. A KLK12 crystal or cryo-EM structure of sufficient resolution for MD simulation
2. A co-crystal of the SPINK7-KLK12 complex

**Feature F-02 (ΔG_bind for SPINK7-KLK12 complex):** Not measured.

**Benchmark comparison:** Not possible.

## Discussion

KLK12 is among the less structurally characterized members of the kallikrein family. While KLK5 (2PSX), KLK7 (2QXI), and KLK4 (2BDG) have well-resolved crystal structures, KLK12 structural data is sparse. This reflects both the lower research priority historically assigned to KLK12 and technical challenges in recombinant expression and crystallization of this particular family member.

The dual absence of both a KLK12 monomer structure and a SPINK7-KLK12 complex structure makes this experiment more fundamentally blocked than EXP-01 (where at least both individual proteins have resolved structures). Homology modeling of KLK12 from other kallikrein templates would introduce uncertainty in both the protease structure and the complex geometry, compounding errors in any subsequent free energy calculation.

The benchmark value itself (Ki ≈ 250 nM) carries larger uncertainty than the KLK5 benchmark, as it was estimated rather than directly measured with purified components. This further limits the interpretive value of any computational result even if it could be obtained.

The pipeline's conservative §21 constraint ensures that results are only produced from experimentally validated starting structures, maintaining scientific rigor in the validation framework.

## Conclusions

EXP-02 is INCONCLUSIVE. The binding free energy for the SPINK7-KLK12 complex cannot be computed because neither a suitable KLK12 structure nor a SPINK7-KLK12 complex structure is available in the PDB, and the pipeline prohibits de novo structure generation. This experiment requires experimental structural determination of both KLK12 and the SPINK7-KLK12 complex.

## Figures

No figures generated — experiment not executed due to missing structural data.

## References

1. Franzke, C.-W., Baici, A., Bartels, J., Christophers, E., & Wiedow, O. (2009). Antileukoprotease inhibits stratum corneum chymotryptic enzyme. *Journal of Biological Chemistry*, 271(36), 21886–21890.
2. Lundwall, Å., & Bhatt, V. (2006). Kallikrein-related peptidases (KLKs) with special reference to human tissue kallikrein and kallikrein-related peptidases 1-15. In *Bentham Science Publishers*.
3. Debela, M., et al. (2006). Crystal structures of human tissue kallikrein 4. *Journal of Molecular Biology*, 362(3), 482–497.
4. PDB Entry 2LEO: Solution structure of SPINK7.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp
