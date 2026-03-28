# EXP-21: KLK5 Subsite Specificity (S1-S4)

**Experiment ID:** EXP-21  
**Feature ID:** F-21 (benchmarks.md)  
**Category:** Structural  
**Status:** QUALITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates that the pipeline's KLK5 model correctly defines the subsites (S1-S4 and S1'-S3') that interact with the inhibitor binding loop residues P1-P4 and P1'-P3'. Trypsin-like serine proteases (including KLK5) have a conserved S1 pocket with Asp189 that selects for basic P1 residues (Arg/Lys). The experiment uses structural analysis of the equilibrated complex to verify correct subsite geometry, pocket depth, and electrostatic character. Cross-reference: Debela et al. 2006 (KLK5 specificity), Michael et al. 2005, Goettig et al. 2010.

---

## 2. Hypothesis

**H₁:** The S1 specificity pocket of KLK5 will contain Asp189 (or equivalent) at the base, forming a salt bridge with the P1 Arg/Lys of the inhibitor.

**H₂:** Subsites S2, S3, and S4 will form contacts with P2, P3, and P4 respectively, with contact distances < 4.5 Å.

**H₃:** The subsite architecture will be consistent with trypsin-like specificity (preference for basic P1).

---

## 3. Background and Rationale

KLK5 is a trypsin-like kallikrein with strong preference for Arg/Lys at P1 (Km = 0.20 mM for Boc-VPR-AMC). Understanding the subsite specificity is essential for interpreting how SPINK7 and LEKTI achieve selective inhibition. The subsites define which residues of the inhibitor binding loop are in direct contact with the protease, and correct subsite mapping is prerequisite for interpreting mutagenesis results (EXP-30) and energy decomposition (EXP-07).

---

## 4. Experimental Protocol

### 4.1 System

SPINK7-KLK5 complex from EXP-01 production trajectory (or BPTI-trypsin from EXP-04 for reference).

### 4.2 Subsite Mapping Protocol

1. **Define subsites structurally:** For each inhibitor residue P4-P1-P1'-P3', identify all protease residues within 4.5 Å (any heavy atom contact).
2. **Classify contacts:** Backbone-backbone, backbone-sidechain, sidechain-sidechain.
3. **S1 pocket analysis:**
   - Identify Asp189 (or equivalent in KLK5).
   - Measure d(Asp189 Oδ, P1 Arg/Lys Nη/Nζ) — salt bridge distance.
   - Measure S1 pocket volume using cavity detection.
4. **S2-S4 pocket analysis:**
   - List all protease residues lining each subsite.
   - Compute pocket-P residue contact area.

### 4.3 Time-Averaged Contacts

From last 50 ns of production MD:
1. Compute per-frame contact maps (P4-P3' vs all protease residues, 4.5 Å cutoff).
2. Persistent contacts: present in >50% of frames.
3. Construct consensus subsite map.

### 4.4 Simulation Parameters (from config.py)

- Force field: AMBER ff14SB
- Duration: 100 ns
- Temperature: 310 K
- Save interval: 10 ps

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Criterion |
|---------------|-----------|
| **PASS** | S1 Asp189 salt bridge maintained (<4.0 Å, >80% occupancy); all subsites correctly mapped |
| **MARGINAL** | S1 correct, but 1–2 subsites ambiguous | 
| **FAIL** | S1 specificity pocket not formed or P1 not in S1 |

---

## 5. Control Conditions

### 5.1 Positive Control: Trypsin Subsite Map

BPTI-trypsin (PDB 2PTC) has extensively characterized subsites. The pipeline's map from EXP-04 MD should match known trypsin subsite assignments from Schechter & Berger nomenclature.

### 5.2 Cross-Validation

Compare KLK5 subsite residue identities with published KLK5 specificity studies (Debela et al. 2006). At minimum, Asp189 (S1), and key S2/S4 residues should match.

---

## 6. Expected Outcomes

| Subsite | Protease Residue(s) | Inhibitor Residue | Key Interaction |
|---------|---------------------|-------------------|-----------------|
| S1 | Asp189 (+ S1 walls) | P1 (Arg/Lys) | Salt bridge |
| S2 | Variable (Leu/His/Gly) | P2 | Hydrophobic/polar |
| S3 | Exposed, less defined | P3 | Backbone H-bonds |
| S4 | Variable | P4 | Peripheral |
| S1' | Adjacent to catalytic triad | P1' | Backbone contacts |
| S2'-S3' | Extended surface | P2'-P3' | Variable |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **KLK5 model quality** | Subsite residues mispositioned | Homology model error | High |
| **Dynamic subsites** | Contact residues change across frames | Flexible protease loops | Medium |
| **Numbering convention** | Residue numbering mismatch with literature | Chymotrypsinogen numbering required | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Asp189 present | KLK5 model contains D189 (or equivalent) | Confirmed in sequence |
| S1 salt bridge | d(D189 Oδ, P1 N) < 4.0 Å in starting structure | Correct docking |
| Contact map | Heatmap shows clear P↔S correlation | Diagonal pattern |
| Trypsin comparison | KLK5 subsite map consistent with trypsin | Same S1, similar S2 |
| Literature match | Key subsite residues match Debela 2006 | Published agreement |

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
