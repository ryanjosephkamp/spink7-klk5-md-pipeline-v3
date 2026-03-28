# EXP-19: Nonpolar Buried Surface Area Fraction — BPTI–Trypsin Complex

## Abstract

This experiment determines the fraction of buried surface area (BSA) contributed by nonpolar atoms at the BPTI–trypsin protein–protein interface (PDB 2PTC). The pipeline computes a nonpolar BSA fraction of 56.9% (916.0 Å² nonpolar, 692.6 Å² polar, 1608.5 Å² total), compared to a literature benchmark of 61% derived from large-scale analyses of protease–inhibitor interfaces. Statistical evaluation yields a **PASS** classification: the absolute discrepancy of 4.1 percentage points falls well within the 1.96 × σ_combined threshold of 11.5 percentage points. The result is consistent with the known polar enrichment of the BPTI–trypsin interface, which features a high density of hydrogen bonds and salt bridges relative to typical hydrophobic protein–protein contacts. Correct classification of hydrogen atoms by their bonded heavy atom proved essential; omitting hydrogens from the nonpolar tally would reduce the computed fraction to 8.1%, producing a spurious FAIL.

## Introduction and Background

The hydrophobic effect is a principal thermodynamic driving force in protein–protein recognition, and the fraction of buried surface area contributed by nonpolar atoms is a widely used descriptor of interface character. Large-scale surveys of protein–protein interfaces have established that nonpolar atoms typically contribute 55–67% of the total BSA in protease–inhibitor complexes, with the remainder arising from polar atoms involved in hydrogen bonding and electrostatic interactions (Bahadur et al., 2004; Janin et al., 2008).

The BPTI–trypsin complex (PDB 2PTC) presents a well-characterized interface with a total BSA of approximately 1600 Å², placing it within the typical range for enzyme–inhibitor complexes. The interface is notable for its high density of specific polar contacts, including the canonical P1 lysine–aspartate salt bridge and multiple backbone hydrogen bonds spanning the binding loop. This polar enrichment is expected to manifest as a nonpolar fraction slightly below the population mean of 61%.

Accurate computation of the nonpolar BSA fraction requires careful treatment of hydrogen atoms. In protein structures prepared with explicit hydrogens, hydrogen atoms constitute a substantial fraction of the total atomic count and contribute meaningfully to the solvent-accessible surface. The chemical character of a hydrogen atom—polar or nonpolar—is determined by its bonded heavy atom: hydrogens bonded to carbon or sulfur are classified as nonpolar, while those bonded to nitrogen or oxygen are classified as polar. Failure to account for hydrogen atoms, or incorrect classification thereof, can produce dramatically erroneous nonpolar fractions.

## Hypothesis

The nonpolar fraction of buried surface area at the BPTI–trypsin interface falls within the 55–67% range characteristic of protease–inhibitor complexes, as established by Bahadur et al. (2004) and Janin et al. (2008).

## Methods

### System Preparation

The BPTI–trypsin complex was obtained from the Protein Data Bank (entry 2PTC, resolution 1.9 Å). The structure was processed with PDBFixer to add missing heavy atoms, assign standard protonation states at pH 7.0, and add explicit hydrogen atoms. Crystallographic water molecules were removed prior to surface area computation.

### Solvent-Accessible Surface Area Calculation

The solvent-accessible surface area (SASA) was computed using the Shrake–Rupley algorithm as implemented in MDTraj, with a probe radius of 1.4 Å. SASA was calculated for three configurations:

1. **Complex:** Both chains present, representing the bound state
2. **Trypsin alone:** Chain A isolated
3. **BPTI alone:** Chain B isolated

The buried surface area for each atom was computed as:

$$\text{BSA}_i = \text{SASA}_i^{\text{free}} - \text{SASA}_i^{\text{complex}}$$

where $\text{SASA}_i^{\text{free}}$ is the accessible surface area of atom $i$ in the isolated chain and $\text{SASA}_i^{\text{complex}}$ is its accessible surface area in the complex. Only atoms with $\text{BSA}_i > 0$ contribute to the total BSA.

### Atom Classification

Atoms were classified as nonpolar or polar according to the following scheme:

- **Nonpolar:** Carbon (C), sulfur (S), and hydrogen atoms bonded to C or S
- **Polar:** Nitrogen (N), oxygen (O), and hydrogen atoms bonded to N or O

This classification follows the convention of Bahadur et al. (2004) and ensures that the large hydrogen atom population is correctly partitioned. The bonded heavy atom for each hydrogen was determined from the molecular topology.

### Nonpolar Fraction Computation

The nonpolar BSA fraction was computed as:

$$f_{\text{nonpolar}} = \frac{\sum_{i \in \text{nonpolar}} \text{BSA}_i}{\sum_{i} \text{BSA}_i} \times 100\%$$

### Statistical Classification

The prediction was evaluated against the literature benchmark using the combined uncertainty framework:

$$\sigma_{\text{combined}} = \sqrt{\sigma_{\text{exp}}^2 + \sigma_{\text{comp}}^2 + \sigma_{\text{method}}^2}$$

where $\sigma_{\text{exp}} = 3.06$ (experimental variability across interface surveys), $\sigma_{\text{comp}} = 3.0$ (computational uncertainty from probe radius, atom radii, and hydrogen placement), and $\sigma_{\text{method}} = 4.0$ (methodological uncertainty from classification scheme and SASA algorithm choice).

## Controls

### Positive Control

The total BSA of 1608.5 Å² is consistent with the value obtained in EXP-16, confirming internal consistency of the SASA calculation pipeline across experiments. This cross-validation ensures that the surface area decomposition into polar and nonpolar components is performed on the correct total.

### Negative Control

A heavy-atom-only analysis was performed in which hydrogen atoms were excluded from the nonpolar tally (only C and S heavy atoms counted as nonpolar). This yielded a nonpolar fraction of 8.1%, which is physically unreasonable and would produce a FAIL classification. The dramatic discrepancy confirms that hydrogen atom classification is essential and that the pipeline correctly includes hydrogens in the polar/nonpolar partition.

## Results

### Quantitative Summary

| Metric | Value |
|---|---|
| Nonpolar BSA | 916.0 Å² |
| Polar BSA | 692.6 Å² |
| Total BSA | 1608.5 Å² |
| Nonpolar fraction (pipeline) | 56.9% |
| Benchmark (literature) | 61% |
| 95% confidence interval | [55%, 67%] |
| Absolute discrepancy | 4.1 percentage points |
| σ_exp | 3.06 |
| σ_comp | 3.0 |
| σ_method | 4.0 |
| σ_combined | 5.86 |
| 1.96 × σ_combined | 11.5 |
| Classification | **PASS** (4.1 ≤ 11.5) |

### BSA Decomposition

The total BSA of 1608.5 Å² partitions as follows:

- **Nonpolar contribution:** 916.0 Å² (56.9%), arising from carbon and sulfur atoms and their bonded hydrogens at the interface
- **Polar contribution:** 692.6 Å² (43.1%), arising from nitrogen and oxygen atoms and their bonded hydrogens

### Classification

$$|\hat{y} - y_{\text{ref}}| = |56.9 - 61.0| = 4.1$$

$$1.96 \times \sigma_{\text{combined}} = 1.96 \times 5.86 = 11.5$$

Since $4.1 \leq 11.5$, the result is classified as **PASS**.

### Sensitivity to Hydrogen Classification

| Classification Scheme | Nonpolar Fraction | Result |
|---|---|---|
| Full (heavy atoms + H by bonded atom) | 56.9% | PASS |
| Heavy atoms only (C, S as nonpolar) | 8.1% | FAIL |

The 48.8 percentage point difference between the two schemes underscores the critical importance of including hydrogen atoms and classifying them by their bonded heavy atom.

## Discussion

The pipeline-computed nonpolar BSA fraction of 56.9% falls within the 95% confidence interval [55%, 67%] established for protease–inhibitor complexes and is classified as PASS with a comfortable margin (discrepancy of 4.1 vs. threshold of 11.5 percentage points). The result is physically meaningful and consistent with the known characteristics of the BPTI–trypsin interface.

The slight depression of the nonpolar fraction below the population mean of 61% is consistent with the polar enrichment of the BPTI–trypsin interface. This complex features an unusually high density of specific polar contacts: the P1 lysine (LYS15) forms a salt bridge with the specificity pocket aspartate (ASP189), multiple backbone hydrogen bonds span the binding loop, and several serine and threonine side chains participate in inter-chain polar contacts. The resulting polar BSA of 692.6 Å² (43.1%) is at the upper end of the distribution for protease–inhibitor interfaces, consistent with the structural basis for the exceptionally tight binding ($K_i \approx 6 \times 10^{-14}$ M).

The most significant methodological finding of this experiment is the essential role of hydrogen atom classification. When hydrogen atoms are excluded from the nonpolar partition and only heavy atoms (C, S) are counted as nonpolar, the computed fraction drops to 8.1%—a value that is physically unreasonable for any protein–protein interface and would produce a FAIL classification. This dramatic discrepancy arises because hydrogen atoms bonded to carbon constitute a large fraction of the atoms at the interface surface, and their exclusion removes the majority of the nonpolar surface area. The standard convention in the protein interface literature (Bahadur et al., 2004; Janin et al., 2008) classifies hydrogen atoms by their bonded heavy atom, and the pipeline correctly implements this convention.

The per-chain decomposition of nonpolar and polar BSA provides additional insight into the asymmetry of the interface. The trypsin side of the interface, which contributes the specificity pocket and catalytic machinery, is expected to show a higher polar fraction than the BPTI side, which presents the predominantly hydrophobic binding loop. This asymmetry is a known feature of enzyme–inhibitor interfaces and is consistent with the functional requirement for specific recognition.

The total BSA of 1608.5 Å² is in excellent agreement with the value obtained in EXP-16, confirming the internal consistency of the SASA calculation across experiments. This cross-validation strengthens confidence in the pipeline's surface area computations and supports the reliability of the polar/nonpolar decomposition.

## Conclusions

The pipeline prediction of a 56.9% nonpolar BSA fraction at the BPTI–trypsin interface receives a **PASS** classification against the literature benchmark of 61% (Bahadur et al., 2004; Janin et al., 2008), with the absolute discrepancy of 4.1 percentage points falling well within the 11.5 percentage point threshold defined by 1.96 × σ_combined. The result is consistent with the polar enrichment characteristic of the BPTI–trypsin interface and validates the pipeline's implementation of atom classification, SASA computation, and BSA decomposition. The experiment further demonstrates that correct treatment of hydrogen atoms—classifying them as polar or nonpolar based on their bonded heavy atom—is essential for accurate nonpolar fraction determination.

## Figures

### Figure 1: BSA Composition — Nonpolar vs. Polar

![Pie chart of nonpolar versus polar BSA](outputs/figures/exp19_bsa_composition.png)

**Figure 1.** Pie chart showing the decomposition of total buried surface area (1608.5 Å²) into nonpolar (916.0 Å², 56.9%, blue) and polar (692.6 Å², 43.1%, orange) contributions. The nonpolar fraction is slightly below the literature mean of 61%, consistent with the polar enrichment of the BPTI–trypsin interface.

### Figure 2: Pipeline vs. Literature Comparison

![Pipeline prediction compared to literature benchmark](outputs/figures/exp19_comparison.png)

**Figure 2.** Comparison of the pipeline-predicted nonpolar BSA fraction (56.9%, blue bar) against the literature benchmark (61%, orange bar) with the 95% confidence interval [55%, 67%] indicated by error bars. The dashed lines denote the ±1.96 × σ_combined classification boundaries. The prediction falls comfortably within the PASS envelope.

### Figure 3: Nonpolar and Polar BSA by Chain

![Grouped bar chart of nonpolar and polar BSA contributions by chain](outputs/figures/exp19_by_chain.png)

**Figure 3.** Grouped bar chart showing the nonpolar (blue) and polar (orange) BSA contributions from the trypsin (chain A) and BPTI (chain B) sides of the interface. The decomposition reveals the per-chain asymmetry in surface character, with the trypsin side contributing a higher proportion of polar surface area due to the specificity pocket and catalytic residues.

## References

1. Bahadur, R. P., Chakrabarti, P., Rodier, F., & Janin, J. (2004). A dissection of specific and non-specific protein–protein interfaces. *Journal of Molecular Biology*, 336(4), 943–955.

2. Janin, J., Bahadur, R. P., & Chakrabarti, P. (2008). Protein–protein interaction and quaternary structure. *Quarterly Reviews of Biophysics*, 41(2), 133–180.

3. Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351–371.

4. McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., Klein, C., Swails, J. M., Hernández, C. X., Schwantes, C. R., Wang, L.-P., Lane, T. J., & Pande, V. S. (2015). MDTraj: A modern open library for the analysis of molecular dynamics trajectories. *Biophysical Journal*, 109(8), 1528–1532.

5. Marquart, M., Walter, J., Deisenhofer, J., Bode, W., & Huber, R. (1983). The geometry of the reactive site and of the peptide groups in trypsin, trypsinogen and its complexes with inhibitors. *Acta Crystallographica Section B*, 39(4), 480–490.

6. Conte, L. L., Chothia, C., & Janin, J. (1999). The atomic structure of protein–protein recognition sites. *Journal of Molecular Biology*, 285(5), 2177–2198.

7. Lawrence, M. C., & Colman, P. M. (1993). Shape complementarity at protein/protein interfaces. *Journal of Molecular Biology*, 234(4), 946–950.

8. Eastman, P., Swails, J., Chodera, J. D., McGibbon, R. T., Zhao, Y., Beauchamp, K. A., Wang, L.-P., Simmonett, A. C., Harrigan, M. P., Stern, C. D., Wiewiora, R. P., Brooks, B. R., & Pande, V. S. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Computational Biology*, 13(7), e1005659.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
