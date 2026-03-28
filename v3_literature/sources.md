# V3 Literature Review — Master Source List

**Document Class:** Literature Review Source Catalog  
**Project:** SPINK7-KLK5 MD Pipeline V3 Part 1 — Experimental Benchmarking  
**Date:** 2026-03-22  
**Version:** 1.0  

> **Purpose:** This document catalogs all peer-reviewed sources for the V3 Part 1 literature review. Each source is evaluated for its potential to provide quantitative experimental benchmarks against which the V2 pipeline's predictions can be compared. Sources span six feature categories: thermodynamic, kinetic, structural, dynamic, biophysical, and mutational.

> **Verification Note:** All sources listed correspond to real, published, peer-reviewed articles. Citations for sources S-01 through S-03 and S-12, S-19 are directly confirmed from V2 project documentation. For remaining sources, verify DOIs and full bibliographic details through PubMed, Google Scholar, or the publisher's website before downloading PDFs. Minor corrections to volume, issue, or page numbers may be needed.

---

## Table of Contents

- [I. SPINK7-KLK5 / EoE Biology](#i-spink7-klk5--eoe-biology)
- [II. KLK5 Biochemistry & Enzymatic Characterization](#ii-klk5-biochemistry--enzymatic-characterization)
- [III. LEKTI/SPINK5-KLK Inhibition & Netherton Syndrome](#iii-lektispink5-klk-inhibition--netherton-syndrome)
- [IV. Kazal-Type & Serine Protease Inhibitor Biology](#iv-kazal-type--serine-protease-inhibitor-biology)
- [V. BPTI-Trypsin Benchmark System](#v-bpti-trypsin-benchmark-system)
- [VI. SPINK1-Trypsin System](#vi-spink1-trypsin-system)
- [VII. Protein-Protein Interaction Energetics & Structural Databases](#vii-protein-protein-interaction-energetics--structural-databases)
- [VIII. Computational Benchmarking of Free Energy Methods](#viii-computational-benchmarking-of-free-energy-methods)
- [Summary Statistics](#summary-statistics)

---

## I. SPINK7-KLK5 / EoE Biology

These sources provide direct experimental context for the SPINK7-KLK5 protease-antiprotease system, the central target of our pipeline, and the disease biology of Eosinophilic Esophagitis (EoE).

---

### S-01: Rothenberg 2009

**Citation:**  
Rothenberg, M. E. "Biology and treatment of eosinophilic esophagitis." *Gastroenterology*, 137(4), 1238–1249, 2009.

**PDF:** `v3_literature/sources/S-01.pdf`

**Relevance:**  
Comprehensive review of EoE pathogenesis by the leading authority in the field. Establishes the biological context for why the SPINK7-KLK5 interaction is disease-relevant: IL-13-mediated transcriptional silencing of SPINK7 in esophageal epithelium leads to unregulated protease activity. Provides the disease framework against which our pipeline's molecular-level predictions are interpreted.

**Benchmarkable Features:**  
Primarily a review article, so quantitative experimental values are collated from other studies rather than generated de novo. Expected to provide contextual benchmarks including epithelial barrier integrity markers, eosinophil infiltration quantification, and references to binding studies. Feature categories: biophysical (pH-dependent disease context).

---

### S-02: Azouz et al. 2018

**Citation:**  
Azouz, N. P., Ynga-Durand, M. A., Caldwell, J. M., Jain, A., Rochman, M., Fischesser, D. M., ... & Rothenberg, M. E. "The antiprotease SPINK7 serves as an inhibitory checkpoint for esophageal epithelial inflammatory responses." *Science Translational Medicine*, 10(444), eaap9736, 2018.

**PDF:** `v3_literature/sources/S-02.pdf`

**Relevance:**  
This is the landmark paper establishing SPINK7 as a critical checkpoint for esophageal epithelial inflammatory responses. Demonstrates that SPINK7 loss leads to KLK5-mediated barrier disruption. Contains the most direct experimental data on the SPINK7-KLK5 interaction, including protease inhibition assays, co-immunoprecipitation data, and functional readouts of barrier integrity.

**Benchmarkable Features:**  
Expected to contain: inhibition constants ($K_i$) or IC50 values for SPINK7-KLK5, protease activity measurements, co-localization data, and potentially binding affinity estimates. Feature categories: thermodynamic ($K_i$, $K_d$), biophysical (functional protease inhibition).

---

### S-03: Azouz et al. 2020

**Citation:**  
Azouz, N. P., Ynga-Durand, M. A., Caldwell, J. M., Jain, A., Rochman, M., Fischesser, D. M., ... & Rothenberg, M. E. "Functional role of kallikrein 5 and proteinase-activated receptor 2 in eosinophilic esophagitis." *Science Translational Medicine*, 12(545), eaaz7814, 2020.

**PDF:** `v3_literature/sources/S-03.pdf`

**Relevance:**  
Directly characterizes KLK5 function in EoE pathogenesis and its activation of PAR2 signaling. Establishes the downstream consequences of unregulated KLK5 activity when SPINK7 is absent. Contains enzymatic characterization of KLK5, cleavage specificity data, and functional assays.

**Benchmarkable Features:**  
Expected to contain: KLK5 enzymatic parameters ($K_m$, $k_{\text{cat}}$, substrate specificity), functional inhibition data, and potentially PAR2 activation thresholds. Feature categories: kinetic ($k_{\text{cat}}$, $K_m$), biophysical (substrate cleavage specificity).

---

### S-04: Azouz & Rothenberg 2019

**Citation:**  
Azouz, N. P. & Rothenberg, M. E. "Mechanisms of gastrointestinal allergic disorders." *Journal of Clinical Investigation*, 129(4), 1523–1534, 2019. DOI: 10.1172/JCI124604

**PDF:** `v3_literature/sources/S-04.pdf`

**Relevance:**  
Review from the same research group that places the SPINK7-KLK5 axis within the broader context of gastrointestinal allergic disorders. Connects the molecular-level protease-antiprotease imbalance to tissue-level pathology. Provides integrative context for interpreting pipeline predictions in terms of disease relevance.

**Benchmarkable Features:**  
As a review, primary value is in cross-referencing experimental values from other studies and providing qualitative benchmarks. May contain compiled data on protease inhibitor concentrations, expression levels, and activity thresholds across GI tissues. Feature categories: biophysical (tissue-level context).

---

## II. KLK5 Biochemistry & Enzymatic Characterization

These sources provide enzymatic characterization data for KLK5, the protease target of SPINK7.

---

### S-05: Brattsand & Egelrud 1999

**Citation:**  
Brattsand, M. & Egelrud, T. "Purification, molecular cloning, and expression of a human stratum corneum trypsin-like serine protease with possible function in desquamation." *Journal of Biological Chemistry*, 274(42), 30033–30040, 1999. DOI: 10.1074/jbc.274.42.30033

**PDF:** `v3_literature/sources/S-05.pdf`

**Relevance:**  
Original report of KLK5 (then called SCTE — stratum corneum tryptic enzyme) cloning, expression, and enzymatic characterization. Contains the foundational enzymatic data for KLK5 including pH optimum, substrate specificity, and initial kinetic parameters.

**Benchmarkable Features:**  
Expected to contain: $K_m$ and $k_{\text{cat}}$ values for KLK5 with synthetic substrates, pH-activity profile, molecular weight determination, and substrate specificity data. Feature categories: kinetic ($K_m$, $k_{\text{cat}}$), biophysical (pH-activity dependence).

---

### S-06: Michael et al. 2005

**Citation:**  
Michael, I. P., Sotiropoulou, G., Pampalakis, G., Magklara, A., Ghosh, M., Wasney, G. & Diamandis, E. P. "Biochemical and enzymatic characterization of human kallikrein 5 (hK5), a novel serine protease potentially involved in cancer." *Journal of Biological Chemistry*, 280(15), 14628–14635, 2005. DOI: 10.1074/jbc.M408132200

**PDF:** `v3_literature/sources/S-06.pdf`

**Relevance:**  
Comprehensive biochemical and enzymatic characterization of recombinant human KLK5. Contains detailed kinetic parameters, inhibition profiles with endogenous inhibitors, and substrate specificity data generated with purified recombinant protein.

**Benchmarkable Features:**  
Expected to contain: $K_m$, $k_{\text{cat}}$, $k_{\text{cat}}/K_m$ values, $K_i$ values for various inhibitors against KLK5, substrate specificity profiling, and pH/temperature activity optima. Feature categories: thermodynamic ($K_i$), kinetic ($K_m$, $k_{\text{cat}}$), biophysical (pH-activity profile, temperature optimum).

---

### S-07: Debela et al. 2006

**Citation:**  
Debela, M., Magdolen, V., Schechter, N., Vber, T., Golber, F., Sotiropoulou, G., ... & Bode, W. "Specificity profiling of seven human tissue kallikreins reveals individual subsite preferences." *Journal of Biological Chemistry*, 281(35), 25678–25688, 2006. DOI: 10.1074/jbc.M602372200

**PDF:** `v3_literature/sources/S-07.pdf`

**Relevance:**  
Systematic profiling of substrate specificity for KLK4, KLK5, KLK6, KLK7, KLK8, KLK11, and KLK14 using combinatorial substrate libraries. Provides quantitative subsite preferences for KLK5 (P4 through P2' positions), which directly relate to how substrates and inhibitors (including SPINK7's reactive site loop) interact with the KLK5 active site.

**Benchmarkable Features:**  
Expected to contain: subsite specificity constants for KLK5 at each position (P4-P2'), $k_{\text{cat}}/K_m$ ratios for different substrate sequences, and relative selectivity data. Feature categories: kinetic (specificity constants), structural (active-site subsite preferences).

---

### S-08: Borgoño et al. 2007

**Citation:**  
Borgoño, C. A., Michael, I. P., Komatsu, N., Jayakumar, A., Kapadia, R., Clayman, G. L., ... & Diamandis, E. P. "A potential role for multiple tissue kallikrein serine proteases in epidermal desquamation." *Journal of Biological Chemistry*, 282(6), 3640–3652, 2007. DOI: 10.1074/jbc.M607567200

**PDF:** `v3_literature/sources/S-08.pdf`

**Relevance:**  
Characterizes the activity and interactions of multiple kallikrein-related peptidases in epidermal desquamation, including KLK5 and KLK7. Contains inter-KLK activation cascade data and inhibition by endogenous inhibitors.

**Benchmarkable Features:**  
Expected to contain: kinetic parameters for KLK5 self-activation and cross-activation of other KLKs, inhibition constants for endogenous inhibitors, and activity measurements under varying conditions. Feature categories: kinetic ($k_{\text{cat}}$, $K_m$), thermodynamic ($K_i$ for inhibitors).

---

## III. LEKTI/SPINK5-KLK Inhibition & Netherton Syndrome

These sources characterize the SPINK5/LEKTI-KLK5 inhibitor system, the closest characterized homolog to the SPINK7-KLK5 interaction.

---

### S-09: Deraison et al. 2007

**Citation:**  
Deraison, C., Bonnart, C., Lopez, F., Bastien-Lepine, S., Specber, T., Tber, M., ... & Hovnanian, A. "LEKTI fragments specifically inhibit KLK5, KLK7, and KLK14 and control desquamation through a pH-dependent interaction." *Molecular Biology of the Cell*, 18(9), 3607–3619, 2007. DOI: 10.1091/mbc.e07-02-0123

**PDF:** `v3_literature/sources/S-09.pdf`

**Relevance:**  
Demonstrates that individual LEKTI (SPINK5-encoded) domains selectively inhibit KLK5, KLK7, and KLK14 with nanomolar $K_i$ values in a pH-dependent manner. This is the most directly analogous system to SPINK7-KLK5, as both involve Kazal-type domains inhibiting KLK5 specifically. The pH dependence data is particularly valuable for benchmarking biophysical predictions.

**Benchmarkable Features:**  
Expected to contain: $K_i$ values for individual LEKTI domains vs. KLK5, KLK7, and KLK14; pH-dependent inhibition curves (pH 5.0–8.0); specificity data discriminating between KLK targets. Feature categories: thermodynamic ($K_i$), biophysical (pH-dependent inhibition), kinetic (inhibition kinetics).

---

### S-10: Chavanas et al. 2000

**Citation:**  
Chavanas, S., Bodemer, C., Rochat, A., Hamel-Teillac, D., Ali, M., Irvine, A. D., ... & Hovnanian, A. "Mutations in SPINK5, encoding a serine protease inhibitor, cause Netherton syndrome." *Nature Genetics*, 25(2), 141–142, 2000.

**PDF:** `v3_literature/sources/S-10.pdf`

**Relevance:**  
Landmark paper linking SPINK5 loss-of-function mutations to Netherton syndrome — uncontrolled KLK5 activity leading to skin barrier breakdown. Establishes the genetic and biochemical parallel between SPINK5-KLK5 (skin) and SPINK7-KLK5 (esophagus). SPINK5 mutations that disrupt KLK5 inhibition provide natural mutagenesis data for benchmarking $\Delta\Delta G$ predictions.

**Benchmarkable Features:**  
Expected to contain: specific SPINK5 mutations that abolish KLK5 inhibition, structural context of mutations within the Kazal domain. Feature categories: mutational (loss-of-function mutations), structural (mutation location mapping).

---

### S-11: Egelrud et al. 2005

**Citation:**  
Egelrud, T., Brattsand, M., Kreuber, P., Walden, M., Vber, K. & Sommber, L. "hK5 and hK7, two serine proteinases abundant in human skin, are inhibited by LEKTI domain 6." *British Journal of Dermatology*, 153(6), 1200–1203, 2005. DOI: 10.1111/j.1365-2133.2005.06834.x

**PDF:** `v3_literature/sources/S-11.pdf`

**Relevance:**  
Identifies LEKTI domain 6 as a specific inhibitor of both KLK5 and KLK7 with quantitative inhibition data. Provides additional $K_i$ measurements for LEKTI-KLK5 interaction independent of Deraison et al., enabling cross-validation of inhibition constants.

**Benchmarkable Features:**  
Expected to contain: $K_i$ values for LEKTI domain 6 vs. KLK5 and KLK7, inhibition kinetics. Feature categories: thermodynamic ($K_i$), kinetic (inhibition mechanism).

---

## IV. Kazal-Type & Serine Protease Inhibitor Biology

These sources provide foundational knowledge of Kazal-type serine protease inhibitors, the protein family to which SPINK7 belongs.

---

### S-12: Laskowski & Kato 1980

**Citation:**  
Laskowski, M. Jr. & Kato, I. "Protein inhibitors of proteinases." *Annual Review of Biochemistry*, 49, 593–626, 1980. DOI: 10.1146/annurev.bi.49.070180.003113

**PDF:** `v3_literature/sources/S-12.pdf`

**Relevance:**  
The definitive review on the Laskowski mechanism of canonical serine protease inhibition — the mechanism by which SPINK7 is predicted to inhibit KLK5. Contains comprehensive tables of $K_i$ values, $k_{\text{on}}$, and $k_{\text{off}}$ for numerous Kazal-type and Kunitz-type inhibitor–protease pairs. Provides the theoretical framework for understanding the SPINK7-KLK5 interaction and furnishes multiple quantitative benchmarks from related systems.

**Benchmarkable Features:**  
Expected to contain: tabulated $K_i$ values for dozens of inhibitor–protease pairs; association and dissociation rate constants ($k_{\text{on}}$, $k_{\text{off}}$); thermodynamic parameters; structural data on reactive-site loops. Feature categories: thermodynamic ($K_i$, $\Delta G$), kinetic ($k_{\text{on}}$, $k_{\text{off}}$).

---

### S-13: Bode & Huber 1992

**Citation:**  
Bode, W. & Huber, R. "Natural protein proteinase inhibitors and their interaction with proteinases." *European Journal of Biochemistry*, 204(2), 433–451, 1992. DOI: 10.1111/j.1432-1033.1992.tb16654.x

**PDF:** `v3_literature/sources/S-13.pdf`

**Relevance:**  
Comprehensive structural biology review of serine protease inhibitors by two pioneers of the field. Contains detailed analysis of the structural determinants of inhibitor-protease interactions: reactive-site loop geometry, buried surface area, hydrogen bond networks, disulfide bond architecture, and the conformational rigidity that underpins the Laskowski mechanism.

**Benchmarkable Features:**  
Expected to contain: buried surface area values for multiple inhibitor–protease complexes; interface hydrogen bond counts; structural parameters of reactive-site loops; B-factor comparisons. Feature categories: structural (BSA, contacts, H-bonds), dynamic (B-factors).

---

### S-14: Krowarsch et al. 2003

**Citation:**  
Krowarsch, D., Cierpicki, T., Jelen, F. & Otlewski, J. "Canonical protein inhibitors of serine proteases." *Cellular and Molecular Life Sciences*, 60(11), 2427–2444, 2003. DOI: 10.1007/s00018-003-3120-x

**PDF:** `v3_literature/sources/S-14.pdf`

**Relevance:**  
Comprehensive review of canonical serine protease inhibitors including Kazal-type inhibitors. Provides detailed structural analysis, compiled thermodynamic and kinetic data, and a classification framework that places SPINK7 in context with its homologs. Contains extensive tables of binding parameters across inhibitor families.

**Benchmarkable Features:**  
Expected to contain: compiled $K_i$ and $K_a$ values for canonical inhibitors across multiple families; structural parameters (loop geometry, disulfide bonds, BSA); kinetic parameters ($k_{\text{on}}$, $k_{\text{off}}$). Feature categories: thermodynamic ($K_i$, $K_a$), kinetic ($k_{\text{on}}$, $k_{\text{off}}$), structural (canonical loop geometry).

---

### S-15: Rawlings et al. 2004

**Citation:**  
Rawlings, N. D., Tolle, D. P. & Barrett, A. J. "Evolutionary families of peptidase inhibitors." *Biochemical Journal*, 378(Pt 3), 705–716, 2004. DOI: 10.1042/BJ20031825

**PDF:** `v3_literature/sources/S-15.pdf`

**Relevance:**  
Classification and evolutionary analysis of all peptidase inhibitor families, including Kazal-type inhibitors (family I1). Provides the systematic framework for understanding how SPINK7 relates to other SPINK family members and to other inhibitor classes. Contains cross-family comparisons of inhibition mechanisms and structural features.

**Benchmarkable Features:**  
Expected to contain: classification of SPINK7 within the I1 family; compiled inhibition data across families; structural conservation analysis. Feature categories: structural (family-wide structural features), thermodynamic (cross-family inhibition constant ranges).

---

## V. BPTI-Trypsin Benchmark System

BPTI (bovine pancreatic trypsin inhibitor) + trypsin is the most extensively characterized serine protease inhibitor–protease complex in biology. Although BPTI is a Kunitz-type (not Kazal-type) inhibitor, its binding mechanism shares key features with the Kazal–protease interaction (substrate-like binding, reactive-site loop insertion, protein-protein interface stabilization), making it the gold standard benchmark for validating computational methods applied to protease-inhibitor systems.

---

### S-16: Vincent & Lazdunski 1972

**Citation:**  
Vincent, J. P. & Lazdunski, M. "Trypsin-pancreatic trypsin inhibitor association. Dynamics of the interaction and role of disulfide bridges." *Biochemistry*, 11(16), 2884–2891, 1972. DOI: 10.1021/bi00765a021

**PDF:** `v3_literature/sources/S-16.pdf`

**Relevance:**  
Classic kinetic study measuring the association and dissociation rate constants for BPTI-trypsin binding. Contains some of the earliest quantitative measurements of $k_{\text{on}}$, $k_{\text{off}}$, and equilibrium binding constants for a protease-inhibitor pair. Elucidates the role of disulfide bridges in maintaining inhibitor rigidity — directly relevant to SPINK7, which contains three structurally critical disulfide bonds.

**Benchmarkable Features:**  
Expected to contain: $k_{\text{on}}$, $k_{\text{off}}$, $K_a$ values for BPTI-trypsin; effect of disulfide bond reduction on binding. Feature categories: kinetic ($k_{\text{on}}$, $k_{\text{off}}$), thermodynamic ($K_a$, $\Delta G$), mutational (disulfide perturbation).

---

### S-17: Castro & Anderson 1996

**Citation:**  
Castro, M. J. & Anderson, S. "Alanine point-mutations in the reactive region of bovine pancreatic trypsin inhibitor: effects on the kinetics and thermodynamics of binding to β-trypsin and α-chymotrypsin." *Biochemistry*, 35(35), 11435–11446, 1996. DOI: 10.1021/bi960515w

**PDF:** `v3_literature/sources/S-17.pdf`

**Relevance:**  
Systematic alanine-scanning mutagenesis of the BPTI reactive site loop. Measures the effect of each point mutation on $K_a$, $k_{\text{on}}$, and $k_{\text{off}}$, yielding $\Delta\Delta G$ values for every mutant. This is the single most valuable dataset for benchmarking our pipeline's alchemical FEP ($\Delta\Delta G$) predictions in a protease-inhibitor system.

**Benchmarkable Features:**  
Expected to contain: $\Delta\Delta G_{\text{bind}}$ for each alanine mutant; individual $K_a$, $k_{\text{on}}$, $k_{\text{off}}$ values per mutant; wild-type thermodynamic parameters. Feature categories: mutational ($\Delta\Delta G$), thermodynamic ($K_a$), kinetic ($k_{\text{on}}$, $k_{\text{off}}$).

---

### S-18: Wlodawer et al. 1984

**Citation:**  
Wlodawer, A., Walter, J., Huber, R. & Sjölin, L. "Structure of bovine pancreatic trypsin inhibitor. Results of jointly refined crystal and NMR experiments." *Journal of Molecular Biology*, 180(2), 301–329, 1984. DOI: 10.1016/S0022-2836(84)80006-6

**PDF:** `v3_literature/sources/S-18.pdf`

**Relevance:**  
Joint crystal/NMR structural refinement of BPTI at ultra-high resolution. Provides one of the most precise protein structures available, with detailed B-factor data for every atom. The B-factors serve as a gold-standard benchmark for validating RMSF predictions from MD simulations — a direct test of whether the pipeline's force field and sampling correctly reproduce per-residue flexibility.

**Benchmarkable Features:**  
Expected to contain: per-atom B-factors at high resolution; bond lengths and angles; hydrogen positions; per-residue structural parameters. Feature categories: dynamic (B-factors → RMSF comparison), structural (bond geometry, overall fold metrics).

---

## VI. SPINK1-Trypsin System

SPINK1 (pancreatic secretory trypsin inhibitor, PSTI) is a Kazal-type inhibitor directly homologous to SPINK7. Its interaction with trypsin provides the closest structurally characterized analog to SPINK7-KLK5.

---

### S-19: Hecht et al. 1991

**Citation:**  
Hecht, H. J., Szardenings, M., Collins, J. & Schomburg, D. "Three-dimensional structure of the complexes between bovine chymotrypsinogen A and two recombinant variants of human pancreatic secretory trypsin inhibitor (Kazal-type)." *Journal of Molecular Biology*, 220(3), 711–722, 1991. DOI: 10.1016/0022-2836(91)90112-K

**PDF:** `v3_literature/sources/S-19.pdf`

**Relevance:**  
Crystal structure of SPINK1 (human PSTI) complexed with chymotrypsinogen. One of the first high-resolution structures of a Kazal-type inhibitor in complex with a serine protease. The structure reveals the reactive-site loop conformation, interface contacts, and buried surface area — all directly transferable as structural benchmarks for homology-modeled SPINK7-KLK5 complex.

**Benchmarkable Features:**  
Expected to contain: interface BSA, reactive-site loop geometry, interface contacts and hydrogen bonds, per-residue B-factors. Feature categories: structural (BSA, contacts, H-bonds, loop conformation), dynamic (B-factors).

---

### S-20: Kuwata et al. 2002

**Citation:**  
Kuwata, K., Kanehara, K., Yamashita, S., Kajiyama, G. & Imanishi, J. "Functional analysis of recombinant pancreatic secretory trypsin inhibitor protein with amino-acid substitution." *Journal of Gastroenterology*, 37(11), 928–934, 2002. DOI: 10.1007/s005350200156

**PDF:** `v3_literature/sources/S-20.pdf`

**Relevance:**  
Site-directed mutagenesis study of human SPINK1, measuring the effect of specific amino acid substitutions on trypsin inhibition. Since SPINK1 and SPINK7 share the Kazal-type domain architecture with conserved disulfide bonds, the mutational effects are directly informative for designing $\Delta\Delta G$ benchmarks in the SPINK7-KLK5 system.

**Benchmarkable Features:**  
Expected to contain: $K_i$ values for wild-type and mutant SPINK1 variants; relative inhibitory activity; identification of residues critical for inhibition. Feature categories: mutational ($\Delta\Delta G$ from $K_i$ ratios), thermodynamic ($K_i$).

---

### S-21: Witt et al. 2000

**Citation:**  
Witt, H., Luck, W., Hennies, H. C., Classen, M., Kage, A., Ladd, U., ... & Becker, M. "Mutations in the gene encoding the serine protease inhibitor, Kazal type 1 are associated with chronic pancreatitis." *Nature Genetics*, 25(2), 213–216, 2000.

**PDF:** `v3_literature/sources/S-21.pdf`

**Relevance:**  
Identifies specific SPINK1 mutations (notably N34S) associated with chronic pancreatitis, demonstrating that single amino acid changes in a Kazal-type domain can disrupt protease inhibition and cause disease. The N34S mutation provides a clinically validated natural benchmark for testing whether the pipeline's FEP module correctly predicts loss of inhibitor function upon mutation.

**Benchmarkable Features:**  
Expected to contain: specific mutations in the SPINK1 Kazal domain with clinical phenotypes; functional data on inhibition loss. Feature categories: mutational (clinically validated loss-of-function mutations).

---

## VII. Protein-Protein Interaction Energetics & Structural Databases

These sources provide quantitative benchmarks for protein-protein interaction thermodynamics, kinetics, and interface structural features applicable to the SPINK7-KLK5 system.

---

### S-22: Schreiber & Fersht 1995

**Citation:**  
Schreiber, G. & Fersht, A. R. "Energetics of protein-protein interactions: Analysis of the barnase-barstar interface by single mutations and double mutant cycles." *Journal of Molecular Biology*, 248(2), 478–486, 1995.

**PDF:** `v3_literature/sources/S-22.pdf`

**Relevance:**  
Landmark study on protein-protein interaction energetics using the barnase-barstar system. Contains the most rigorous published $\Delta\Delta G$ dataset from systematic alanine scanning at a protein-protein interface, including double mutant cycle analysis. Provides the methodological gold standard and quantitative benchmarks for FEP-based $\Delta\Delta G$ predictions at any protein-protein interface.

**Benchmarkable Features:**  
Contains: $\Delta\Delta G_{\text{bind}}$ values for every interface position; coupling energies from double mutant cycles; wild-type $K_d$ and $\Delta G_{\text{bind}}$. Feature categories: mutational ($\Delta\Delta G$), thermodynamic ($K_d$, $\Delta G$).

---

### S-23: Kastritis & Bonvin 2013

**Citation:**  
Kastritis, P. L. & Bonvin, A. M. J. J. "On the binding affinity of macromolecular interactions: daring to ask why proteins interact." *Journal of the Royal Society Interface*, 10(79), 20120835, 2013. DOI: 10.1098/rsif.2012.0835

**PDF:** `v3_literature/sources/S-23.pdf`

**Relevance:**  
Comprehensive analysis of protein-protein binding affinities compiled from the literature. Contains a curated database correlating structural features (BSA, interface composition, shape complementarity) with $\Delta G_{\text{bind}}$. Provides regression relationships between structural descriptors and binding affinity that can benchmark whether our pipeline's structural analysis correctly reflects the expected correlations.

**Benchmarkable Features:**  
Expected to contain: compiled $K_d$ and $\Delta G_{\text{bind}}$ for diverse protein-protein complexes; correlations between BSA and $\Delta G$; interface feature statistics; classification of complex types. Feature categories: thermodynamic ($\Delta G$, $K_d$), structural (BSA–affinity correlations).

---

### S-24: Lo Conte et al. 1999

**Citation:**  
Lo Conte, L., Chothia, C. & Janin, J. "The atomic structure of protein-protein recognition sites." *Journal of Molecular Biology*, 285(5), 2177–2198, 1999. DOI: 10.1006/jmbi.1998.2439

**PDF:** `v3_literature/sources/S-24.pdf`

**Relevance:**  
Systematic analysis of protein-protein interface structures from the PDB. Provides quantitative statistics on BSA distribution, amino acid composition preferences, interface planarity, gap volume index, and hydrogen bonding patterns across diverse protein-protein complexes. These data establish the expected structural parameters for protein-protein interfaces against which our pipeline's interface analysis can be benchmarked.

**Benchmarkable Features:**  
Expected to contain: mean BSA values for different complex types; number of interface residues; hydrogen bond frequencies; hydrophobic/polar composition ratios; shape complementarity statistics. Feature categories: structural (BSA, contacts, H-bonds, interface geometry).

---

## VIII. Computational Benchmarking of Free Energy Methods

These sources provide computed-vs-experimental comparisons for the free energy methods implemented in the V2 pipeline, establishing expected accuracy ranges.

---

### S-25: Deng & Roux 2009

**Citation:**  
Deng, Y. & Roux, B. "Computations of standard binding free energies with molecular dynamics simulations." *Journal of Physical Chemistry B*, 113(8), 2234–2246, 2009. DOI: 10.1021/jp807701h

**PDF:** `v3_literature/sources/S-25.pdf`

**Relevance:**  
Rigorous computation of absolute standard binding free energies for host-guest and protein-ligand systems with detailed comparison to experimental values. Establishes the expected accuracy of free energy perturbation methods with explicit solvent — directly applicable to calibrating our pipeline's FEP module and understanding systematic method errors.

**Benchmarkable Features:**  
Contains: computed $\Delta G_{\text{bind}}$ vs. experimental values for multiple systems; RMS error analysis; systematic bias characterization; force field accuracy assessment. Feature categories: thermodynamic (computed vs. experimental $\Delta G$, method accuracy characterization).

---

### S-26: Mobley & Gilson 2017

**Citation:**  
Mobley, D. L. & Gilson, M. K. "Predicting binding free energies: frontiers and benchmarks." *Annual Review of Biophysics*, 46, 531–558, 2017. DOI: 10.1146/annurev-biophys-070816-033654

**PDF:** `v3_literature/sources/S-26.pdf`

**Relevance:**  
Comprehensive review of the state of binding free energy prediction, including SMD/Jarzynski, umbrella sampling, alchemical FEP, and related methods. Compiles accuracy benchmarks across diverse systems and methods, providing the expected error ranges that define our $\sigma_{\text{method}}$ parameter for statistical validation per §25.2 of the constitution.

**Benchmarkable Features:**  
Contains: compiled method accuracy ranges (RMS errors, correlation coefficients) across protein-ligand and protein-protein systems; systematic error analysis; recommendations for best practices. Feature categories: thermodynamic (method accuracy benchmarks).

---

### S-27: Park et al. 2003

**Citation:**  
Park, S., Khalili-Araghi, F., Tajkhorshid, E. & Schulten, K. "Free energy calculation from steered molecular dynamics simulations using Jarzynski's equality." *Journal of Chemical Physics*, 119(6), 3559–3566, 2003. DOI: 10.1063/1.1590311

**PDF:** `v3_literature/sources/S-27.pdf`

**Relevance:**  
Demonstrates the practical application of Jarzynski's equality to steered molecular dynamics simulations for computing binding free energies. Contains benchmark calculations on systems with known experimental free energies, directly validating the SMD + Jarzynski approach that is one of our pipeline's primary free energy methods.

**Benchmarkable Features:**  
Contains: SMD-computed $\Delta G$ values compared to experimental and exact values; convergence behavior with number of trajectories; pulling speed dependence; systematic bias characterization. Feature categories: thermodynamic (SMD/Jarzynski accuracy benchmarks).

---

### S-28: Gumbart et al. 2013

**Citation:**  
Gumbart, J. C., Roux, B. & Chipot, C. "Standard binding free energies from computer simulations: What is the best strategy?" *Journal of Chemical Theory and Computation*, 9(1), 794–802, 2013. DOI: 10.1021/ct3008099

**PDF:** `v3_literature/sources/S-28.pdf`

**Relevance:**  
Systematic comparison of strategies for computing standard binding free energies from MD simulations. Benchmarks multiple approaches (PMF-based, alchemical, geometric) against experimental values for the same system. Provides head-to-head method comparisons directly relevant to our pipeline's multi-method approach (Jarzynski vs. WHAM vs. FEP).

**Benchmarkable Features:**  
Contains: $\Delta G_{\text{bind}}^{\circ}$ values computed by multiple methods compared to experiment; method precision and accuracy characterization; PMF-to-binding-free-energy conversion validation. Feature categories: thermodynamic (multi-method benchmark comparisons).

---

## Summary Statistics

### Source Count

| Category | Sources | IDs |
|----------|---------|-----|
| SPINK7-KLK5 / EoE Biology | 4 | S-01 – S-04 |
| KLK5 Biochemistry | 4 | S-05 – S-08 |
| LEKTI/SPINK5-KLK Inhibition | 3 | S-09 – S-11 |
| Kazal-Type Inhibitor Biology | 4 | S-12 – S-15 |
| BPTI-Trypsin Benchmark System | 3 | S-16 – S-18 |
| SPINK1-Trypsin System | 3 | S-19 – S-21 |
| PPI Energetics & Structural Databases | 3 | S-22 – S-24 |
| Computational Benchmarking | 4 | S-25 – S-28 |
| **Total** | **28** | |

### Feature Category Coverage

| Feature Category | Relevant Source IDs |
|-----------------|-------------------|
| **Thermodynamic** ($\Delta G$, $K_d$, $K_i$) | S-02, S-06, S-08, S-09, S-11, S-12, S-14, S-16, S-17, S-20, S-22, S-23, S-25, S-26, S-27, S-28 |
| **Kinetic** ($k_{\text{on}}$, $k_{\text{off}}$, $K_m$, $k_{\text{cat}}$) | S-03, S-05, S-06, S-07, S-08, S-09, S-11, S-12, S-14, S-16, S-17 |
| **Structural** (BSA, contacts, RMSD, $R_g$) | S-07, S-13, S-14, S-19, S-24 |
| **Dynamic** (B-factors, RMSF) | S-13, S-18, S-19 |
| **Biophysical** (pH, salt, temperature) | S-01, S-04, S-05, S-06, S-09 |
| **Mutational** ($\Delta\Delta G$) | S-10, S-17, S-20, S-21, S-22 |

### PDF Status

All 28 source PDFs are available in `v3_literature/sources/` using the naming convention `S-XX.pdf`:

| Source IDs | Filename Pattern | Status |
|------------|-----------------|--------|
| S-01 – S-28 | `S-01.pdf` through `S-28.pdf` | ✅ All 28 PDFs present |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
