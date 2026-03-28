# EXP-15: Catalytic Contact Distance — Serine Protease Active Site Geometry

## Abstract

This experiment evaluates whether the molecular dynamics pipeline's structural representation of the BPTI–trypsin complex (PDB 2PTC) preserves the catalytic contact distance between the serine protease catalytic residue Ser195 OG and the P1 substrate carbonyl carbon. The pipeline yielded a catalytic contact distance of 2.675 Å, compared to the established crystallographic benchmark of 2.7 Å (Laskowski & Kato, 1980). The discrepancy of 0.025 Å falls well within the 95% confidence interval of [2.4, 3.0] Å, and the combined uncertainty σ_combined = 0.189 Å confirms that the deviation is not statistically significant. The experiment is classified as **PASS**, indicating that the pipeline faithfully reproduces the catalytic geometry essential for serine protease function.

## Introduction and Background

Serine proteases constitute one of the largest and most extensively characterized families of proteolytic enzymes. These enzymes employ a conserved catalytic triad—composed of serine, histidine, and aspartate residues—to hydrolyze peptide bonds via a charge-relay mechanism. In trypsin, the catalytic triad consists of Ser195, His57, and Asp102 (chymotrypsinogen numbering), arranged in a precise spatial configuration that enables nucleophilic attack on substrate peptide bonds.

The distance between the catalytic serine hydroxyl oxygen (OG) and the P1 substrate carbonyl carbon is a critical geometric parameter that reflects the enzyme's readiness for nucleophilic attack. Crystallographic studies have established this distance at approximately 2.7 Å in productive enzyme–substrate and enzyme–inhibitor complexes (Laskowski & Kato, 1980; Marquart et al., 1983). Deviations from this optimal distance indicate either a non-productive binding mode or structural artifacts introduced during computational processing.

The bovine pancreatic trypsin inhibitor (BPTI) in complex with trypsin (PDB 2PTC) serves as a canonical model system for protease–inhibitor interactions. The BPTI reactive-site loop inserts the P1 residue (Lys15) into the trypsin active site, positioning the scissile bond carbonyl carbon in close proximity to Ser195 OG. This complex has been resolved at 1.9 Å resolution and is among the most thoroughly studied protein–protein complexes in structural biology (Marquart et al., 1983).

Accurate preservation of catalytic geometry during system preparation is essential for downstream molecular dynamics simulations. If the pipeline introduces distortions to the active-site distances during structure processing (e.g., through PDBFixer protonation, missing-residue reconstruction, or energy minimization), the resulting trajectories may not faithfully represent the enzyme's catalytic competence. This experiment directly tests whether such distortions occur.

## Hypothesis

The pipeline's structural representation preserves the catalytic contact distance (Ser195 OG to P1 carbonyl C) within 2.4–3.0 Å, consistent with the canonical serine protease mechanism. Formally:

- **H₀**: The pipeline-derived catalytic contact distance falls outside the 95% confidence interval [2.4, 3.0] Å, indicating structural distortion.
- **H₁**: The pipeline-derived catalytic contact distance falls within [2.4, 3.0] Å, indicating faithful preservation of active-site geometry.

## Methods

### System Preparation

The BPTI–trypsin complex was obtained from the Protein Data Bank (PDB ID: 2PTC, resolution 1.9 Å). Structure preparation was performed using PDBFixer (Eastman et al., 2017) with the following protocol:

1. Missing heavy atoms were reconstructed using PDBFixer's template-based rebuilding.
2. Missing hydrogen atoms were added at pH 7.0.
3. Non-standard residues were replaced with standard equivalents where applicable.
4. Crystallographic water molecules were retained.

PDBFixer renumbers residues during processing; consequently, the catalytic triad residues Ser195, His57, and Asp102 (original chymotrypsinogen numbering) correspond to Ser177, His40, and Asp84 in the pipeline-processed structure. Residue identity was verified by spatial proximity analysis rather than sequence numbering alone.

### Distance Measurements

Interatomic distances were measured between the following atom pairs:

| Measurement | Atom 1 | Atom 2 | Description |
|---|---|---|---|
| Catalytic contact | Ser195 OG (Ser177) | P1 carbonyl C | Nucleophilic attack distance |
| Triad H-bond 1 | His57 NE2 (His40) | Ser195 OG (Ser177) | His–Ser hydrogen bond |
| Triad H-bond 2 | Asp102 OD (Asp84) | His57 ND1 (His40) | Asp–His hydrogen bond |

Distances were computed using Euclidean distance in Cartesian coordinate space from the pipeline-processed coordinate file.

### Uncertainty Quantification

Combined uncertainty was calculated by quadrature addition of independent error sources:

$$\sigma_{\text{combined}} = \sqrt{\sigma_{\text{exp}}^2 + \sigma_{\text{comp}}^2 + \sigma_{\text{method}}^2}$$

where σ_exp = 0.153 Å (experimental crystallographic uncertainty), σ_comp = 0.05 Å (computational measurement precision), and σ_method = 0.10 Å (methodological uncertainty from structure preparation).

### Classification Criteria

The experiment is classified as PASS if the pipeline measurement falls within the 95% confidence interval [2.4, 3.0] Å derived from the literature benchmark and combined uncertainty. Otherwise, the experiment is classified as FAIL.

## Controls

**Positive control**: The catalytic triad geometry in the unprocessed PDB 2PTC crystal structure, which reports Ser195 OG–substrate C distances consistent with the 2.7 Å benchmark. Agreement between the raw crystal structure and literature values confirms that the reference system is appropriate.

**Negative control**: A hypothetical structure in which the catalytic serine is displaced beyond 4.0 Å from the substrate, representing a non-productive binding mode. Such a distance would fall outside the confidence interval and produce a FAIL classification, confirming the discriminative power of the metric.

## Results

### Primary Measurement

The pipeline-derived catalytic contact distance (Ser195 OG → P1 substrate carbonyl C) was measured at **2.675 ± 0.189 Å**, compared to the literature benchmark of **2.7 ± 0.153 Å**.

| Quantity | Value |
|---|---|
| Pipeline prediction | 2.675 Å |
| Literature benchmark | 2.7 Å |
| Absolute discrepancy | 0.025 Å |
| Relative discrepancy | 0.93% |
| 95% CI | [2.4, 3.0] Å |
| σ_combined | 0.189 Å |
| **Classification** | **PASS** |

The discrepancy of 0.025 Å represents 0.13σ_combined, well below the 1.96σ threshold for statistical significance at the 95% level.

### Catalytic Triad Distances

All three catalytic triad interatomic distances were measured to assess the overall integrity of the active-site geometry:

| Distance | Pipeline (Å) | Expected (Å) | Status |
|---|---|---|---|
| Ser195 OG → P1 C | 2.675 | 2.7 | Within range |
| His57 NE2 → Ser195 OG | 2.616 | 2.6–2.9 | Within range |
| Asp102 OD → His57 ND1 | 2.625 | 2.5–2.8 | Within range |

All three distances fall within the ranges characteristic of catalytically competent serine protease structures, confirming that the charge-relay network is preserved in its entirety.

### Catalytic Triad Identification

The catalytic triad was identified by spatial proximity analysis in the pipeline-processed structure. Due to PDBFixer residue renumbering, the triad residues correspond to:

| Original (chymotrypsinogen) | Pipeline (renumbered) |
|---|---|
| Ser195 | Ser177 |
| His57 | His40 |
| Asp102 | Asp84 |

Identification was confirmed by verifying that these three residues form a spatially contiguous network of hydrogen-bond-competent distances, consistent with the canonical catalytic triad arrangement.

### Figures

![Catalytic Triad Distances](outputs/figures/exp15_catalytic_distances.png)

**Figure 1.** Bar chart of the three catalytic triad interatomic distances measured from the pipeline-processed BPTI–trypsin structure. All distances fall within the ranges expected for a catalytically competent serine protease. Error bars represent σ_combined for each measurement.

![Pipeline vs Literature Comparison](outputs/figures/exp15_comparison.png)

**Figure 2.** Comparison of the pipeline-derived catalytic contact distance (Ser195 OG → P1 C) with the literature benchmark value of 2.7 Å. Error bars represent the 95% confidence interval. The near-complete overlap confirms quantitative agreement between the pipeline output and established crystallographic data.

![Catalytic Triad Schematic](outputs/figures/exp15_triad_schematic.png)

**Figure 3.** Schematic diagram of the catalytic triad geometry in the BPTI–trypsin complex, showing the spatial arrangement of Ser195 (Ser177), His57 (His40), and Asp102 (Asp84). Dashed lines indicate hydrogen bonds with measured distances annotated. The charge-relay mechanism proceeds from Asp → His → Ser, activating the serine hydroxyl for nucleophilic attack on the substrate P1 carbonyl carbon.

## Discussion

The pipeline-derived catalytic contact distance of 2.675 Å demonstrates excellent agreement with the established benchmark of 2.7 Å. The absolute discrepancy of 0.025 Å is substantially smaller than both the experimental uncertainty (σ_exp = 0.153 Å) and the combined uncertainty (σ_combined = 0.189 Å), indicating that the pipeline's structure preparation protocol preserves catalytic geometry to within crystallographic precision.

The preservation of all three catalytic triad distances is particularly noteworthy. The His–Ser hydrogen bond distance (2.616 Å) and the Asp–His hydrogen bond distance (2.625 Å) both fall within the ranges characteristic of enzymatically active conformations. This demonstrates that the PDBFixer preparation protocol—including protonation state assignment, missing-atom reconstruction, and hydrogen addition—does not introduce distortions into the active-site hydrogen-bonding network.

The residue renumbering introduced by PDBFixer (Ser195→Ser177, His57→His40, Asp102→Asp84) does not affect the structural integrity of the coordinate data but necessitates spatial-proximity-based identification of catalytic residues rather than sequence-number-based lookup. This is an important operational consideration for automated analysis pipelines and underscores the value of geometry-based residue identification algorithms.

The PASS classification for this experiment provides confidence that the pipeline can be trusted to produce structurally faithful representations of enzyme active sites for downstream molecular dynamics simulations of catalytic mechanisms.

## Conclusions

The molecular dynamics pipeline preserves the catalytic contact distance of the BPTI–trypsin serine protease complex to within 0.025 Å of the established crystallographic benchmark, corresponding to a relative discrepancy of 0.93%. All three catalytic triad interatomic distances fall within the ranges characteristic of catalytically competent structures. The experiment is classified as PASS, confirming that the pipeline's structure preparation protocol does not introduce artifacts into the active-site geometry that would compromise the fidelity of subsequent molecular dynamics simulations.

## References

1. Laskowski, M., & Kato, I. (1980). Protein inhibitors of proteinases. *Annual Review of Biochemistry*, 49(1), 593–626.

2. Marquart, M., Walter, J., Deisenhofer, J., Bode, W., & Huber, R. (1983). The geometry of the reactive site and of the peptide groups in trypsin, trypsinogen and its complexes with inhibitors. *Acta Crystallographica Section B*, 39(4), 480–490.

3. Eastman, P., Swails, J., Chodera, J. D., McGibbon, R. T., Zhao, Y., Beauchamp, K. A., ... & Pande, V. S. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Computational Biology*, 13(7), e1005659.

4. Hedstrom, L. (2002). Serine protease mechanism and specificity. *Chemical Reviews*, 102(12), 4501–4524.

5. Blow, D. M. (1976). Structure and mechanism of chymotrypsin. *Accounts of Chemical Research*, 9(4), 145–152.

6. Perona, J. J., & Craik, C. S. (1995). Structural basis of substrate specificity in the serine proteases. *Protein Science*, 4(3), 337–360.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
