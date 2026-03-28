# EXP-31: Disulfide Bond Ablation (ΔΔG ~ +7 kcal/mol)

**Experiment ID:** EXP-31  
**Feature ID:** F-31 (benchmarks.md)  
**Category:** Mutational  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the contribution of the primary scaffold disulfide bond (C14-C38 in BPTI) to binding free energy by computationally ablating it (C14S/C38S double mutation) and measuring ΔΔG. Experimental data shows that removing this disulfide destabilizes binding by ΔΔG ≈ +7 kcal/mol (Krowarsch et al. 2003, Goldenberg 1988), comparable to the total scaffold energy (EXP-28). This validates the pipeline's FEP methodology for structurally significant mutations and confirms the critical role of disulfide bonds in pre-organizing the binding loop.

---

## 2. Hypothesis

**H₁:** FEP calculation of C14S/C38S double mutation in BPTI-trypsin will yield ΔΔG = +7 ± 3 kcal/mol (destabilization).

**H₂:** The structural effect of disulfide removal will be observable as increased binding loop flexibility (RMSF increase >50%) and partial loss of canonical loop geometry.

---

## 3. Background and Rationale

The C14-C38 disulfide bond in BPTI (and its equivalents in SPINK1, SPINK7, and other Kazal inhibitors) directly constrains the binding loop. Removing it eliminates the covalent crosslink that pre-organizes the P1 residue in the correct conformation for protease recognition. The ΔΔG of ~7 kcal/mol is remarkably close to the scaffold energy (EXP-28: −7.9 kcal/mol), suggesting that the disulfide accounts for the majority of the scaffold's thermodynamic contribution.

---

## 4. Experimental Protocol

### 4.1 Mutation: C14S / C38S

Replace both cysteines with serines (conservative: preserves size and polarity, removes only the S-S bond). This is a double mutation requiring a coupled FEP calculation.

### 4.2 FEP Thermodynamic Cycle

$$\Delta\Delta G = \Delta G_{complex}^{WT \to C14S/C38S} - \Delta G_{free}^{WT \to C14S/C38S}$$

### 4.3 FEP Protocol

**Leg 1: Complex (bound state)**
1. Start from equilibrated BPTI-trypsin (EXP-04).
2. Hybrid topology: C14→S14 + C38→S38 (simultaneous transformation).
3. Lambda windows: n_lambda = 20 (0.0 to 1.0).
4. Per-window: 2.0 ns.
5. Temperature: 310 K.
6. Soft-core: alpha = 0.5, power = 1.
7. Save interval: 1.0 ps.
8. Annihilate electrostatics: True.

**Leg 2: Free inhibitor (unbound state)**
1. Start from equilibrated free BPTI.
2. Same FEP protocol — transform C14/C38 to S14/S38.

### 4.4 Analysis

- MBAR: solver "robust", tolerance 10⁻⁷.
- Bootstrap: n_bootstrap = 200.
- Overlap matrix: verify λ-window overlap.

### 4.5 Structural Analysis Post-Ablation

After FEP at λ=1.0 (fully mutated):
1. Run 50 ns production MD of mutant complex (C14S/C38S-trypsin).
2. Monitor:
   - Binding loop RMSD vs WT (expected increase).
   - P1 K15 position in S1 pocket (expected partial displacement).
   - Interface H-bond count (expected decrease).
   - BSA (expected decrease).

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | ΔΔG (kcal/mol) |
|---------------|----------------|
| **PASS** | [+4, +10] |
| **MARGINAL** | [+2, +13] |
| **FAIL** | < +2 or > +13 or wrong sign |

---

## 5. Control Conditions

### 5.1 C30S/C51S Control

BPTI has three disulfide bonds. C30-C51 is the least important for binding loop geometry. Ablating C30-C51 should yield a much smaller ΔΔG (+1–2 kcal/mol), confirming that C14-C38 is specifically important.

### 5.2 EXP-28 Cross-Reference

ΔΔG(C14S/C38S) should be comparable to the scaffold energy from EXP-28 (−7.9 kcal/mol). If they agree within 3 kcal/mol, the disulfide-scaffold energy connection is validated.

### 5.3 Forward/Reverse FEP

Run reverse FEP (S14/S38 → C14/C38 with disulfide formation). Should yield ΔG = −ΔG_forward within uncertainty.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| ΔΔG (C14S/C38S) | +7 kcal/mol | Goldenberg 1988 |
| 95% CI | [+4, +10] | benchmarks.md |
| Loop RMSD increase (mutant) | >50% vs WT | Structural destabilization |
| H-bond count decrease | −2 to −4 bonds | Interface weakening |
| BSA decrease | −100 to −300 Å² | Partial interface loss |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Double mutation FEP complexity** | Poor overlap at intermediate λ | Two simultaneous perturbations | High |
| **Protein unfolding at λ=1** | BPTI partially unfolds without disulfide | Lost structural restraint | High |
| **Disulfide topology handling** | FEP cannot break covalent bonds | Need special treatment of S-S bond | Critical |
| **Endpoint singularity** | Free energy diverges at λ=0 or λ=1 | Soft-core parameters needed | Medium |

### Special Note on Disulfide Breaking

Breaking a covalent S-S bond in FEP requires careful handling:
- Option A: Restrain S-S distance during transformation, compute restraint correction analytically.
- Option B: Use decoupling approach — decouple C14-C38 interactions while coupling S14-S38.
- Option C: Separate the perturbation: first S-S bond release (restraint release ΔG), then C→S side chain change.

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Hybrid topology | C14S/C38S topology builds cleanly | No atom count errors |
| S-S bond handling | Disulfide breaking method validated | ΔG_restraint computed |
| Lambda overlap | All adjacent windows >10% overlap | MBAR reliable |
| Forward/reverse | ΔG_forward ≈ −ΔG_reverse | Within 2 kcal/mol |
| Mutant MD | Binding loop RMSD increases | Structural effect observed |
| EXP-28 agreement | ΔΔG ≈ scaffold energy | Within 3 kcal/mol |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (BPTI-trypsin complex: ~35,000 atoms; FEP with disulfide bond topology changes) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation (C14S/C38S double mutation) | 0.2 | 0.2 | S-S bond removal + Cys→Ser mutation |
| Equilibration (WT + C14S/C38S mutant, bound + free) | 1 | 0.5 | 4 legs: WT-bound, WT-free, mut-bound, mut-free |
| FEP complex leg (20λ × 2 ns) | 8 | 4 | 40 ns; requires soft-core potentials for S-S |
| FEP solvent leg (20λ × 2 ns) | 8 | 4 | 40 ns |
| BAR/MBAR analysis + hysteresis check | 0.5 | 0.5 | CPU-bound |
| **Total** | **~20–24** | **~10–12** | §10.17: Tier 5 ≈ 20–24 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each λ window; after equilibration of each leg |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-31/` |
| Restart procedure | Mount Drive → identify last completed λ window per leg → verify energy drift < 0.1% → resume |
| Estimated sessions needed | 1–2 (A100) or 1 (H100) |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| FEP window state | OpenMM binary `.chk` | `checkpoint_<leg>_lambda<idx>.chk` |
| Full simulation state | OpenMM XML | `state_<leg>_lambda<idx>.xml` |
| Per-window ΔU samples | NumPy `.npy` | `du_<leg>_lambda<idx>.npy` |
| Disulfide topology (before/after) | OpenMM XML | `system_wt.xml`, `system_c14s_c38s.xml` |
| Combined FEP results | NumPy `.npz` | `fep_disulfide_results.npz` |

**Resume verification protocol:**
1. Reload checkpoint: `simulation.loadCheckpoint('checkpoint_<leg>_lambda<idx>.chk')`.
2. Run 1000 steps; compute energy.
3. Compare to energy at checkpoint save: drift must be < 0.1%.
4. If drift exceeds threshold, discard window and re-equilibrate at that λ for 200 ps.

**Special note:** Disulfide bond removal requires topology modification (remove HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce between Cβ14-Cβ38). Use soft-core Lennard-Jones potentials for alchemical intermediates to avoid singularities.

### 9.5 Platform Configuration

```python
from openmm import Platform
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
# Soft-core potentials recommended for disulfide FEP:
# Use CustomNonbondedForce with alpha=0.5, power=1 for LJ soft-core
```

**Environment verification:**
```bash
!nvidia-smi
python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"
```

### 9.6 Dependency Notes

**Upstream dependency on EXP-04:**

| Required EXP-04 Output | Drive Path | Description |
|-----------------------|------------|-------------|
| Equilibrated BPTI-trypsin complex (2PTC) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/equilibrated_complex.pdb` | Starting structure for C14S/C38S mutation |
| System XML (force field definitions) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/system/system.xml` | Topology template for disulfide modification |

**Pre-execution check:** EXP-04 must be completed. The C14S/C38S mutations are applied to the EXP-04 equilibrated structure.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
