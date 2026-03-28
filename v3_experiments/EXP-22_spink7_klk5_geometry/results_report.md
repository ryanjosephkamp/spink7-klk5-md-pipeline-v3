# EXP-22: SPINK7-KLK5 Complex Geometry — Results Report

## Abstract

This report documents the outcome of EXP-22, which aimed to characterize the interface geometry of the SPINK7-KLK5 complex, including buried surface area (BSA), hydrogen bond count, and shape complementarity. The experiment is classified as **INCONCLUSIVE** because no SPINK7-KLK5 complex structure exists in the PDB, and the pipeline's §21 (Pipeline Immutability) constraint prohibits de novo molecular docking. This experiment depends directly on EXP-01, which is also INCONCLUSIVE for the same structural reason.

## Introduction / Background

Protein-protein interface geometry is a fundamental descriptor of complex stability and specificity. For Kazal-type inhibitor-protease complexes, key geometric features include:

- **Buried Surface Area (BSA):** The total solvent-accessible surface area occluded upon complex formation, typically 1200–1800 Å² for Kazal-protease pairs.
- **Hydrogen bonds:** Inter-molecular hydrogen bonds at the interface, typically ≥8 for stable Kazal-protease complexes.
- **Shape complementarity (Sc):** A measure of geometric fit between the interacting surfaces, typically 0.65–0.80 for protease-inhibitor complexes.

These geometric benchmarks are estimated from characterized Kazal-protease complexes such as BPTI-trypsin (PDB 2PTC) and other Kunitz/Kazal-protease co-crystals. The SPINK7-KLK5 complex is expected to exhibit similar interface geometry given the conserved canonical binding loop mechanism.

## Hypothesis

**H-22:** The SPINK7-KLK5 complex exhibits interface geometry consistent with canonical Kazal-protease interactions: BSA = 1200–1800 Å², H-bonds ≥ 8, and Sc = 0.65–0.80.

## Methods

**Planned approach:** Extract interface geometry metrics from equilibrated MD trajectories of the SPINK7-KLK5 complex using MDTraj and FreeSASA. Compute BSA as the difference in SASA between the complex and isolated components. Count inter-molecular hydrogen bonds using standard geometric criteria (donor-acceptor distance ≤ 3.5 Å, angle ≥ 120°). Evaluate shape complementarity using the Lawrence & Colman algorithm.

**Pipeline components:** OpenMM 8.1, AMBER14/SB force field, TIP3P water model, MDTraj, FreeSASA.

**Status:** Methods could not be executed. No SPINK7-KLK5 complex structure was available as input.

## Controls

Not applicable — experiment was not executed.

## Results

**Status: INCONCLUSIVE — experiment could not be executed.**

No SPINK7-KLK5 complex structure exists in the PDB. This experiment depends directly on EXP-01 (SPINK7-KLK5 binding free energy), which is also INCONCLUSIVE for the identical structural reason. Interface geometry analysis (BSA, H-bond count, shape complementarity) requires a validated complex structure as input — either an experimental co-crystal or a complex equilibrated from an experimental starting structure.

The pipeline cannot create de novo docked complexes per §21 (Pipeline Immutability).

**Feature F-22 (Interface geometry of SPINK7-KLK5 complex):** Not measured.

**Benchmark comparison:** Not possible.

## Discussion

Interface geometry analysis is among the most straightforward computational analyses once a complex structure is available — it requires only geometric calculations on static or equilibrated structures, without the sampling challenges inherent in free energy calculations. The pipeline has successfully performed analogous analyses on available protease-inhibitor complexes:

- **EXP-16, EXP-17, EXP-19:** Interface geometry analyses on the BPTI-trypsin complex (PDB 2PTC) have been completed, validating the pipeline's geometric analysis tools against well-characterized structural benchmarks.

The blocking factor for EXP-22 is entirely the absence of a starting complex structure, not any methodological limitation. Once a SPINK7-KLK5 co-crystal or cryo-EM structure is deposited in the PDB, this experiment can be executed using the same validated tools and protocols applied to BPTI-trypsin.

The estimated benchmarks (BSA = 1200–1800 Å², H-bonds ≥ 8) are derived from the broader family of Kazal/Kunitz-type inhibitor-protease complexes. The actual SPINK7-KLK5 interface may deviate from these estimates if SPINK7 utilizes non-canonical exosite contacts or if the KLK5 surface topology differs significantly from trypsin at the inhibitor binding interface. Such deviations would themselves be scientifically informative, underscoring the value of this experiment when structural data becomes available.

## Conclusions

EXP-22 is INCONCLUSIVE. The interface geometry of the SPINK7-KLK5 complex cannot be characterized because no experimentally resolved complex structure exists in the PDB, and the pipeline prohibits de novo complex construction per §21. This experiment is directly dependent on EXP-01 and will become executable when an experimental SPINK7-KLK5 complex structure is available. The pipeline's geometric analysis tools have been validated on the BPTI-trypsin system (2PTC) and are ready for application to SPINK7-KLK5 when the structural prerequisite is met.

## Figures

No figures generated — experiment not executed due to missing structural data.

## References

1. Lawrence, M. C., & Colman, P. M. (1993). Shape complementarity at protein/protein interfaces. *Journal of Molecular Biology*, 234(4), 946–950.
2. Krissinel, E., & Henrick, K. (2007). Inference of macromolecular assemblies from crystalline state. *Journal of Molecular Biology*, 372(3), 774–797.
3. Franzke, C.-W., et al. (2009). Antileukoprotease inhibits stratum corneum chymotryptic enzyme. *Journal of Biological Chemistry*, 271(36), 21886–21890.
4. PDB Entry 2PTC: Crystal structure of BPTI-trypsin complex.
5. PDB Entry 2LEO: Solution structure of SPINK7.
6. PDB Entry 2PSX: Crystal structure of KLK5.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp
