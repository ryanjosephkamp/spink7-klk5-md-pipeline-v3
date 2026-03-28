# EXP-21 Results Report: Protease Subsite Specificity Analysis (Crystal Structure)

**Experiment ID:** EXP-21
**Feature:** F-21 — Characterizes the protease subsite architecture and specificity-determining contacts
**Classification:** PASS
**Date:** 2026-03-23

---

## Abstract

This experiment characterizes the subsite specificity architecture of a trypsin-fold serine protease using the bovine trypsin–BPTI crystal structure (PDB: 2PTC) as a same-fold reference for KLK5, whose experimental complex structure is unavailable. A distance-based contact mapping approach ($d \leq 4.5$ Å) identifies protease residues lining each inhibitor subsite from P4 through P3'. The S1 pocket emerges as the dominant specificity determinant with 14 contacts and a confirmed Lys15 NZ–Asp171 OD salt bridge at 3.63 Å, consistent with trypsin's canonical preference for basic P1 residues. These results establish the baseline subsite geometry expected for KLK5 and other trypsin-like kallikreins.

---

## Introduction / Background

Serine proteases of the trypsin fold cleave substrates according to subsite–substrate complementarity, where the nomenclature $P_n \cdots P_1 \!\downarrow\! P_1' \cdots P_n'$ (Schechter–Berger) denotes substrate positions relative to the scissile bond, and $S_n$ denotes the corresponding enzyme pockets. The S1 subsite is the primary specificity determinant for trypsin-like proteases: an aspartate residue (Asp189 in chymotrypsinogen numbering) at the base of the S1 pocket forms a salt bridge with basic residues (Arg or Lys) at the P1 position, conferring the hallmark tryptic specificity.

Kallikrein-related peptidase 5 (KLK5) is a trypsin-fold serine protease implicated in epidermal desquamation and Netherton syndrome pathophysiology. Because no experimentally resolved KLK5–inhibitor complex is available, the well-characterized bovine $\beta$-trypsin–BPTI complex (PDB 2PTC, 1.9 Å resolution) serves as a structural surrogate. Both enzymes share the chymotrypsin-like fold, the catalytic triad (Ser–His–Asp), and a conserved S1 pocket containing an aspartate.

**Residue numbering note:** All residue numbers in this report reflect PDBFixer renumbering. The correspondence to canonical chymotrypsinogen numbering is:

| PDBFixer | Chymotrypsinogen |
|----------|-----------------|
| SER177   | Ser195          |
| HIS40    | His57           |
| ASP84    | Asp102          |
| ASP171   | Asp189          |

---

## Hypothesis

The S1 subsite of trypsin (and, by extension, KLK5) is the deepest and most contact-rich pocket, containing an aspartate residue that forms a salt bridge ($d < 4.0$ Å) with a basic P1 residue, thereby governing primary specificity.

---

## Methods

### Structural data

The crystal structure of bovine $\beta$-trypsin in complex with bovine pancreatic trypsin inhibitor (BPTI) was obtained from the Protein Data Bank (PDB ID: 2PTC). The structure was preprocessed with PDBFixer for residue renumbering and protonation-state assignment.

### P1 identification

The P1 residue was identified as the inhibitor residue whose $C_\alpha$ atom lies closest to the catalytic serine O$\gamma$ (SER177 OG). This yielded **Lys15** of BPTI, consistent with the known reactive-site loop assignment.

### Contact mapping

For each inhibitor position $P_n$ ($n \in \{4, 3, 2, 1, 1', 2', 3'\}$), all protease heavy atoms within a distance cutoff of $d_{\text{cut}} = 4.5$ Å from any heavy atom of the inhibitor residue were recorded. Contacts were aggregated per residue, retaining the minimum heavy-atom distance as the representative contact distance:

$$d_{ij} = \min_{a \in R_i,\, b \in R_j} \| \mathbf{r}_a - \mathbf{r}_b \|$$

where $R_i$ and $R_j$ are the heavy-atom sets of protease residue $i$ and inhibitor residue $j$, respectively.

### S1 salt-bridge validation

The salt bridge was assessed by measuring the distance between Lys15 NZ (P1 side-chain nitrogen) and the nearest carboxylate oxygen of Asp171 (OD1/OD2). A functional salt bridge is defined by $d_{\text{salt}} < 4.0$ Å.

### Software

Analysis was performed with MDAnalysis and BioPython. Figures were generated with Matplotlib and Seaborn.

---

## Controls

| Control | Purpose | Outcome |
|---------|---------|---------|
| P1 assignment via catalytic Ser proximity | Ensures correct reactive-site residue identification | Lys15 identified — agrees with literature |
| Distance cutoff $d = 4.5$ Å | Standard non-bonded contact threshold for crystal structures | Reproduces known subsite assignments |
| Salt-bridge threshold $d < 4.0$ Å | Upper bound for functional ion pairing (literature range: 2.5–4.0 Å) | 3.63 Å — within cutoff |
| Reference structure (2PTC) | High-resolution, well-validated trypsin–inhibitor complex | 1.9 Å resolution; $R_{\text{free}}$ consistent with quality |

---

## Results

### Subsite contact map

The distance-based contact mapping yields the following protease–inhibitor contact profile across the Schechter–Berger positions:

| Position | Residue | # Contacts | Protease residues (distance, Å) |
|----------|---------|-----------|--------------------------------|
| P4  | GLY12  | 1  | GLN174 (3.72) |
| P3  | PRO13  | 3  | GLY194 (3.22), TRP193 (3.48), GLN174 (4.21) |
| P2  | CYS14  | 6  | GLN174 (2.89), HIS40 (3.41), SER192 (3.49), LEU81 (3.71), TRP193 (3.97), SER177 (4.37) |
| **P1**  | **LYS15**  | **14** | SER177 (2.68), GLY175 (2.76), SER172 (3.01), ASP176 (3.06), GLN174 (3.56), SER192 (3.56), TRP193 (3.61), ASP171 (3.63), CYS173 (3.78), HIS40 (3.79), GLY204 (3.81), GLY196 (3.87), GLY194 (3.93), VAL191 (4.14) |
| P1' | ALA16  | 6  | SER177 (3.08), GLY175 (3.42), PHE24 (3.50), GLN174 (3.58), CYS25 (3.66), HIS40 (4.02) |
| P2' | ARG17  | 6  | HIS23 (2.81), PHE24 (2.83), TYR22 (3.31), GLY175 (3.54), TYR131 (3.55), GLN174 (4.30) |
| P3' | ILE18  | 3  | TYR22 (3.89), PHE24 (4.17), HIS40 (4.29) |

The contact distribution shows a clear maximum at P1 (14 contacts), with flanking positions P2, P1', and P2' each contributing 6 contacts. Distal positions (P4, P3, P3') contribute 1–3 contacts each.

### S1 pocket composition

The S1 pocket (contacts to P1 = Lys15) comprises 14 protease residues spanning the catalytic machinery (SER177, HIS40), the specificity-determining aspartate (ASP171), and the pocket-lining backbone (GLY175, GLY194, GLY196, GLY204, SER172, CYS173, ASP176, GLN174, SER192, TRP193, VAL191).

### Salt-bridge validation

$$d(\text{Lys15 NZ} \longleftrightarrow \text{Asp171 OD}) = 3.63 \text{ Å}$$

This distance falls within the accepted range for a functional salt bridge ($2.5 \leq d \leq 4.0$ Å), confirming the electrostatic anchor that defines tryptic specificity.

### Benchmark assessment

| Criterion | Expected | Observed | Status |
|-----------|----------|----------|--------|
| Asp in S1 pocket | Yes | ASP171 present | **PASS** |
| Salt bridge $< 4.0$ Å | Yes | 3.63 Å | **PASS** |
| P1 has most contacts | Yes | 14 (max) | **PASS** |

**Overall classification: PASS**

---

## Discussion

### S1 as the primary specificity determinant

The P1 position (Lys15) engages 14 protease residues — more than any other subsite — reflecting the deep, well-formed S1 pocket that is the hallmark of trypsin-like serine proteases. The pocket buries the P1 side chain and positions it for the specificity-critical salt bridge with ASP171 (Asp189 in chymotrypsinogen numbering). The measured distance of 3.63 Å is within the range for functional ion pairing (typically 2.5–4.0 Å), confirming the electrostatic complementarity that underpins trypsin's strong preference for Arg/Lys at P1.

### Extended subsite architecture

The S2 subsite (contacts to P2 = Cys14) contains 6 contacts, notably including the catalytic HIS40 (His57) and SER192, both of which contribute to the oxyanion hole that stabilizes the tetrahedral intermediate during catalysis. The presence of catalytic residues in S2 contacts is consistent with the positioning of the P2 side chain adjacent to the active-site cleft.

The prime-side subsites S1' and S2' each show 6 contacts, with S2' featuring aromatic residues (PHE24, TYR22, TYR131) that may contribute to hydrophobic packing with the leaving-group segment of the substrate. The distal subsites (S4, S3, S3') are shallow, with only 1–3 contacts, indicating limited specificity enforcement at these positions.

### Implications for KLK5

KLK5 shares the trypsin-like fold and is expected to harbor a similar S1 architecture with an aspartate at the pocket base. The subsite geometry established here — deep S1 with salt-bridge anchoring, moderate S2/S1'/S2' engagement, and shallow distal sites — serves as the structural baseline for interpreting KLK5 specificity. Differences in KLK5's extended subsites (S2–S4, S1'–S3') relative to trypsin are expected to modulate secondary specificity and may contribute to KLK5's restricted substrate repertoire compared to trypsin's broad cleavage profile.

### Limitations

- The analysis is performed on a static crystal structure and does not capture dynamical fluctuations in subsite geometry.
- The 4.5 Å heavy-atom cutoff is a standard but arbitrary threshold; hydrogen-bond versus van der Waals contacts are not distinguished.
- BPTI is a canonical inhibitor, not a substrate; the reactive-site loop geometry may differ from a true Michaelis complex.
- PDBFixer renumbering means residue identifiers differ from the canonical chymotrypsinogen scheme (correspondence table provided in Introduction).

---

## Conclusions

1. The S1 subsite of trypsin (2PTC) is the deepest and most contact-rich pocket, engaging 14 protease residues with the P1 residue (Lys15).
2. A salt bridge between Lys15 NZ and Asp171 OD at 3.63 Å confirms the electrostatic basis of tryptic specificity.
3. The contact gradient — $N_{\text{contacts}}$: P1 (14) $\gg$ P2 = P1' = P2' (6) $>$ P3 = P3' (3) $>$ P4 (1) — quantitatively reflects the dominance of S1 in substrate recognition.
4. These results establish the reference subsite geometry for trypsin-fold proteases and provide the structural baseline for KLK5 specificity analysis.

---

## Figures

### Figure 1: Subsite Contact Heatmap

![Protease–inhibitor contact heatmap](outputs/figures/exp21_contact_map.png)

**Figure 1.** Heatmap of protease–inhibitor contact strengths across all Schechter–Berger subsites (P4–P3'). Color intensity encodes the minimum heavy-atom distance ($d_{ij}$, Å) between each protease residue and the corresponding inhibitor position. The S1 column (P1 = Lys15) shows the largest number and strongest contacts.

### Figure 2: Contact Counts per Inhibitor Position

![Contacts per position bar chart](outputs/figures/exp21_contacts_per_position.png)

**Figure 2.** Bar chart of the number of protease residues within 4.5 Å of each inhibitor position. P1 (Lys15) dominates with 14 contacts, consistent with the deep S1 pocket. Flanking positions P2, P1', and P2' each contribute 6 contacts; distal positions contribute 1–3.

### Figure 3: S1 Salt-Bridge Distance

![S1 salt bridge distance](outputs/figures/exp21_salt_bridge.png)

**Figure 3.** Distance between Lys15 NZ and Asp171 OD (3.63 Å) compared to the functional salt-bridge threshold (4.0 Å, dashed line). The measured distance falls well within the accepted range for ion pairing ($2.5 \leq d \leq 4.0$ Å), confirming the electrostatic anchor of tryptic S1 specificity.

---

## References

1. Schechter, I. & Berger, A. On the size of the active site in proteases. *Biochem. Biophys. Res. Commun.* **27**, 157–162 (1967).
2. Marquart, M., Walter, J., Deisenhofer, J., Bode, W. & Huber, R. The geometry of the reactive site and of the peptide groups in trypsin, trypsinogen and its complexes with inhibitors. *Acta Crystallogr. B* **39**, 480–490 (1983).
3. Perona, J. J. & Craik, C. S. Structural basis of substrate specificity in the serine proteases. *Protein Sci.* **4**, 337–360 (1995).
4. Debela, M. et al. Crystal structures of human tissue kallikreins 5 and 7 and their substrate specificity. *J. Mol. Biol.* **378**, 371–384 (2008).
5. Brattsand, M. et al. A proteolytic cascade of kallikreins in the stratum corneum. *J. Invest. Dermatol.* **124**, 198–203 (2005).
6. Kumar, J. K. Molecular links between skin barrier disruption and Netherton syndrome. *Int. J. Dermatol.* **55**, 239–249 (2016).
7. Hedstrom, L. Serine protease mechanism and specificity. *Chem. Rev.* **102**, 4501–4524 (2002).

---

## Author Block

**Author:** Ryan Kamp
**Affiliation:** Dept. of Computer Science, University of Cincinnati
**Email:** kamprj@mail.uc.edu
**GitHub:** ryanjosephkamp
