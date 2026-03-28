# EXP-03: LEKTI Domain-KLK5 Binding Panel — Results Report

## Abstract

This report documents the outcome of EXP-03, which aimed to compute binding free energies (ΔG_bind) for multiple LEKTI (SPINK5) Kazal-type inhibitor domains with KLK5, creating a binding panel across the LEKTI domain repertoire. The experiment is classified as **INCONCLUSIVE** because no experimentally resolved LEKTI-KLK5 complex structures exist in the PDB, and the pipeline's §21 (Pipeline Immutability) constraint prohibits de novo molecular docking or homology modeling of complexes.

## Introduction / Background

LEKTI (Lympho-Epithelial Kazal-Type Inhibitor), encoded by the SPINK5 gene, is a large multi-domain serine protease inhibitor containing 15 Kazal-type inhibitor domains. LEKTI is critical for skin barrier homeostasis; loss-of-function mutations in SPINK5 cause Netherton syndrome, characterized by uncontrolled kallikrein activity and severe skin inflammation.

Deraison et al. (2007) demonstrated that different LEKTI domains exhibit varying inhibitory potencies against KLK5, with ΔG_bind values ranging from approximately −8.0 to −13.0 kcal/mol across domains. This domain-specific selectivity is biologically important for the graded regulation of KLK5 activity across skin layers.

A computational binding panel reproducing these differential affinities would validate the pipeline's ability to discriminate between closely related inhibitor variants binding the same protease target.

## Hypothesis

**H-03:** The pipeline can compute ΔG_bind for multiple LEKTI domain-KLK5 complexes, reproducing the experimentally observed rank ordering and magnitude range (−8.0 to −13.0 kcal/mol) reported by Deraison et al. (2007).

## Methods

**Planned approach:** For each LEKTI domain with available structural data, perform SMD pulling from the KLK5 active site followed by umbrella sampling with WHAM analysis. Compare computed ΔG_bind values across domains to reproduce the experimental binding panel.

**Pipeline components:** OpenMM 8.1, AMBER14/SB force field, TIP3P water model, PDBFixer, MDTraj/MDAnalysis.

**Status:** Methods could not be executed. No LEKTI domain-KLK5 complex structures were available.

## Controls

Not applicable — experiment was not executed.

## Results

**Status: INCONCLUSIVE — experiment could not be executed.**

No experimentally resolved LEKTI-KLK5 complex structures exist in the PDB. LEKTI (SPINK5) is a large multi-domain protein with 15 Kazal-type inhibitor domains. Structures of individual LEKTI domains have been partially determined (e.g., domain 6), but no co-crystal with KLK5 has been deposited for any domain.

The pipeline cannot generate docked complexes per §21 (Pipeline Immutability). This panel experiment requires at least one validated LEKTI domain-KLK5 co-crystal structure.

**Feature F-03 (ΔG_bind panel for LEKTI domains with KLK5):** Not measured.

**Benchmark comparison:** Not possible.

## Discussion

The absence of LEKTI-KLK5 co-crystal structures reflects multiple structural biology challenges:

1. **Multi-domain complexity:** LEKTI's 15-domain architecture makes full-length crystallization extremely difficult. Individual domain constructs must be designed, expressed, and purified separately.

2. **Domain flexibility:** The inter-domain linkers in LEKTI are flexible, making it unclear whether isolated domain structures accurately represent their conformations in the full-length protein context.

3. **Partial structural coverage:** While some individual LEKTI domain structures have been determined (e.g., LEKTI domain 6), the coverage is incomplete, and no domain has been co-crystallized with KLK5.

4. **Panel scale:** Even if one LEKTI domain-KLK5 complex were available, validating the full binding panel would require structures for multiple domains, multiplying the structural data requirements.

The binding panel experiment is particularly valuable because it tests relative rather than absolute accuracy — the ability to rank-order binding affinities across related inhibitor variants. This type of selectivity prediction is a stringent test of force field accuracy for protein-protein interactions.

The pipeline has demonstrated successful binding free energy calculations on the BPTI-trypsin system (2PTC), which shares the Kazal/Kunitz-protease binding motif. When LEKTI domain-KLK5 structures become available, the methodology validated on BPTI-trypsin can be directly applied.

## Conclusions

EXP-03 is INCONCLUSIVE. The LEKTI domain-KLK5 binding panel cannot be computed because no experimentally resolved LEKTI domain-KLK5 complex structures exist in the PDB, and the pipeline prohibits de novo complex construction. This experiment requires the experimental determination of at least one (ideally multiple) LEKTI domain-KLK5 co-crystal structures.

## Figures

No figures generated — experiment not executed due to missing structural data.

## References

1. Deraison, C., et al. (2007). LEKTI fragments specifically inhibit KLK5, KLK7, and KLK14 and control desquamation through a pH-dependent interaction. *Molecular Biology of the Cell*, 18(9), 3607–3619.
2. Chavanas, S., et al. (2000). Mutations in SPINK5, encoding a serine protease inhibitor, cause Netherton syndrome. *Nature Genetics*, 25(2), 141–142.
3. Descargues, P., et al. (2005). Spink5-deficient mice mimic Netherton syndrome through degradation of desmoglein 1 by epidermal protease hyperactivity. *Nature Genetics*, 37(1), 56–65.
4. PDB Entry 2PSX: Crystal structure of KLK5.

---

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp
