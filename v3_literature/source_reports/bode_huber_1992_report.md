# Source Report: S-13

## Bibliographic Information
- **ID:** S-13
- **Authors:** Bode W, Huber R
- **Title:** Natural Protein Proteinase Inhibitors and Their Interaction with Proteinases
- **Journal:** European Journal of Biochemistry (FEBS)
- **Year:** 1992
- **Volume/Pages:** 204:433–451
- **DOI:** N/A
- **PDF Filename:** S-13.pdf

## Author Block
- **Report Author:** Ryan Kamp
- **Affiliation:** Dept. of Computer Science, University of Cincinnati
- **Contact Email:** kamprj@mail.uc.edu
- **GitHub:** ryanjosephkamp

## Source Category
Category IV: Kazal-Type Inhibitor Biology

## Summary
Comprehensive structural review of protein proteinase inhibitors and their interactions, based on available X-ray crystal structures. Defines the canonical binding loop conformation with specific backbone dihedral angles (φ/ψ) for the P3–P1 residues, the sub-van der Waals catalytic contact (P1 carbonyl carbon to Ser195 Oγ ≈ 2.7 Å), and interface geometry (contact surface 600–900 Å²). Covers Kazal, Kunitz (BPTI), SSI, potato inhibitor, squash seed, serpin, and hirudin-thrombin families. Quantifies scaffolding contribution to binding energy (−33 kJ/mol) and thermal stability (BPTI Tm = 95°C, OMTKY3 Tm = 85°C).

## Extracted Features

### Feature 1: Canonical Binding Loop Backbone Angles
- **Type:** Structural
- **Value:** P3: −140° < φ < −120°, 140° < ψ < 170°; P2: −100° < φ < −60°, 139° < ψ < 180°; P1: −120° < φ < −95°, 9° < ψ < 50°
- **Conditions:** Compiled from crystal structures of multiple inhibitor families
- **Confidence:** High — derived from multiple independent structures
- **Benchmarkable:** Yes — directly measurable from MD trajectories of SPINK7 binding loop

### Feature 2: Catalytic Contact Distance
- **Type:** Structural
- **Value:** P1 carbonyl carbon to Ser195 Oγ ≈ 2.7 Å (sub-van der Waals contact)
- **Conditions:** Crystal structure measurements across families
- **Confidence:** High
- **Benchmarkable:** Yes — key geometric parameter for MD validation

### Feature 3: Interface Contact Surface Area
- **Type:** Structural
- **Value:** 600–900 Å² total contact surface; intermolecular area confined to small surface strip
- **Conditions:** Canonical inhibitor-protease complexes
- **Confidence:** High
- **Benchmarkable:** Yes — measurable from MD trajectories via SASA calculations

### Feature 4: Scaffolding Contribution to Binding
- **Type:** Thermodynamic
- **Value:** −33 kJ/mol (−7.9 kcal/mol) from scaffolding (non-contact residues maintaining loop conformation)
- **Conditions:** Estimated from comparison of small peptide vs. full protein inhibitors
- **Confidence:** Moderate — estimate based on indirect comparison
- **Benchmarkable:** Partially — can be probed by alchemical FEP of scaffold mutations

### Feature 5: Thermal Stability of Inhibitors
- **Type:** Thermodynamic
- **Value:** BPTI Tm = 95°C; OMTKY3 Tm = 85°C; both stable in 6M guanidinium chloride
- **Conditions:** Denaturation experiments
- **Confidence:** High
- **Benchmarkable:** No — relates to inhibitor folding stability, not complex binding

### Feature 6: Typical Association Rate Constant
- **Type:** Kinetic
- **Value:** k_on ≈ 10⁶ M⁻¹s⁻¹
- **Conditions:** General for canonical inhibitors
- **Confidence:** High — consistent across families
- **Benchmarkable:** Partially — can be compared to MD association simulations if performed

### Feature 7: Kazal Family Crystal Structures (from Table 1)
- **Type:** Structural
- **Value:** OMJPQ3, PSTI:BTgen, OMTKY3:SGPB (0.18 nm resolution), OMTKY3:CHT (0.18 nm), OMTKY3:HLE, OMSVP3 — all solved and listed
- **Conditions:** X-ray crystallography
- **Confidence:** High
- **Benchmarkable:** Yes — Kazal-protease complex structures are directly relevant benchmark systems

### Feature 8: Thrombin-Hirudin Extended Interface
- **Type:** Structural
- **Value:** Contact surface 18 nm² (1800 Å²) — much larger than canonical complexes
- **Conditions:** Non-canonical lock-and-key mechanism
- **Confidence:** High
- **Benchmarkable:** Context only — demonstrates range of PPI interface sizes

## PDB Codes Referenced
No explicit PDB accession codes listed (structures referenced by literature citations). Key structures discussed: BPTI, OMTKY3, OMSVP3, OMJPQ3, eglin c, CMTI-I, EETI-II, SSI, STI, ETI, MPI, hirudin-thrombin, stefin-papain, PCI-carboxypeptidase A

## Usefulness Assessment
- **Overall Rating:** High
- **Rationale:** Provides the quantitative structural parameters (canonical loop φ/ψ angles, catalytic contact distance, interface area) that serve as direct validation targets for MD simulations. The canonical loop geometry is the single most important structural benchmark for any Kazal-type inhibitor simulation. The scaffolding energy contribution helps decompose binding free energies. The Kazal-specific structural listings enable identification of benchmark systems.
- **Key Limitation:** No PDB accession codes (1992 publication convention). Structures must be identified from literature citations and cross-referenced against the PDB.
