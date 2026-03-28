# EXP-30: BPTI Reactive Region Alanine Scan ΔΔG Panel

**Experiment ID:** EXP-30  
**Feature ID:** F-30 (benchmarks.md)  
**Category:** Thermodynamic / Alchemical FEP  
**Status:** QUANTITATIVE  
**Date:** 2026-07-17  

---

## 1. Abstract

This experiment validates the pipeline's alchemical free energy perturbation (FEP) methodology by computing ΔΔG values for 15 alanine mutations spanning the BPTI reactive region (T11A through G36A) against trypsin. The Castro & Anderson (1996) dataset provides full kinetic (kon, koff, Ki) and thermodynamic (ΔG, ΔH, TΔS) characterization for each mutant. The dominant P1 residue K15A (ΔΔG ≈ +10 kcal/mol) and the null-effect G36A (ΔΔG ≈ 0) bracket the dynamic range. The pipeline must reproduce the rank correlation (Spearman ρ > 0.7) and per-mutant ΔΔG within the 95% CI of ±2.3 kcal/mol.

---

## 2. Hypothesis

**H₁:** FEP-computed ΔΔG values for the 15-mutation alanine scan panel will achieve a Spearman rank correlation ρ > 0.7 against experimental values from Castro & Anderson (1996).

**H₂:** Per-mutant |ΔΔG_comp − ΔΔG_exp| < 2.3 kcal/mol for ≥ 12 of 15 mutations (80%).

**H₃:** The K15A (P1) mutation will be correctly identified as the most destabilizing mutation, with computed ΔΔG > +7 kcal/mol.

**H₄:** The G36A mutation will be correctly identified as neutral, with computed |ΔΔG| < 1.5 kcal/mol.

---

## 3. Background and Rationale

Castro & Anderson (1996) performed the most comprehensive alanine scanning study of any protease-inhibitor system. Key findings:

1. **P1 dominance:** K15A increases Ki by ~10⁷-fold (5 × 10⁻¹⁴ → 1.4 × 10⁻⁶ M), contributing ΔΔG ≈ +10 kcal/mol.
2. **Association vs. dissociation:** Most mutations affect koff more than kon, indicating that binding affinity is determined primarily by complex stability rather than association kinetics.
3. **Enthalpy-entropy compensation:** ΔH and TΔS changes across mutants are anti-correlated.
4. **Structural context:** The reactive region (residues 11–38) encompasses the binding loop, including the P1 residue (K15) and the stabilizing disulfide Cys14-Cys38.

This experiment is the gold standard for validating alchemical FEP on a protein-protein interface, and the 15-point scan provides statistical power for correlation analysis.

---

## 4. Experimental Protocol

### 4.1 System Preparation

| Parameter | Value |
|-----------|-------|
| Complex structure | PDB 2PTC (BPTI-trypsin, 1.9 Å); use EXP-04 equilibrated structure |
| Force field | AMBER ff14SB (`amber14-all.xml`) |
| Water model | TIP3P (`amber14/tip3p.xml`) |
| pH | 7.4 |
| Box padding | 1.2 nm |
| Ionic strength | 0.15 M NaCl |
| Box shape | Cubic |

### 4.2 Mutation Panel

| # | Mutant | Experimental ΔΔG (kcal/mol) | Ki (M) | Category |
|---|--------|----------------------------|--------|----------|
| 1 | T11A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 2 | G12A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 3 | P13A | ≈ +1–2 | ~10⁻¹³ | Peripheral |
| 4 | C14A | ≈ +7 | ~10⁻⁹ | Disulfide (Cys14-Cys38) |
| 5 | K15A (P1) | ≈ +10 | 1.4 × 10⁻⁶ | P1 dominant |
| 6 | A16G | (Ala→Gly control) | varies | Framework |
| 7 | R17A (P2') | ≈ +5 | ~10⁻¹⁰ | Interface contact |
| 8 | I18A | ≈ +3–4 | ~10⁻¹¹ | Hydrophobic packing |
| 9 | I19A | ≈ +3–4 | ~10⁻¹¹ | Hydrophobic packing |
| 10 | R20A | ≈ +2–3 | ~10⁻¹² | Salt bridge |
| 11 | Y21A | ≈ +2–3 | ~10⁻¹² | Aromatic contact |
| 12 | F22A | ≈ +2–3 | ~10⁻¹¹ | Hydrophobic packing |
| 13 | Y23A | ≈ +1–2 | ~10⁻¹³ | Peripheral aromatic |
| 14 | C38A | ≈ +7 | ~10⁻⁹ | Disulfide (Cys14-Cys38) |
| 15 | G36A | ≈ 0 | ~10⁻¹⁴ | Null control |

### 4.3 FEP Protocol (per mutation)

| Parameter | Value |
|-----------|-------|
| Number of λ windows | 20 (evenly spaced 0.0–1.0) |
| Per-window equilibration | 200 ps |
| Per-window production | 2 ns |
| Total per mutation | 2 legs × 20 windows × 2 ns = 80 ns |
| Thermodynamic cycle | ΔΔG = ΔG_mut(complex) − ΔG_mut(solvent) |
| Alchemical method | Single-topology, soft-core Lennard-Jones |
| Electrostatic staging | Separate electrostatic and LJ decoupling |
| Free energy estimator | BAR (primary), MBAR (cross-check) |

### 4.4 Staged Execution Strategy

Given the large number of mutations, execution is staged by priority:

| Stage | Mutations | Est. Total FEP time | Rationale |
|-------|-----------|---------------------|-----------|
| **Stage 1** (pilot) | K15A, G36A, R17A | 240 ns | Extreme + null + moderate; validate pipeline |
| **Stage 2** (core) | C14A, C38A, I18A, I19A | 320 ns | Structural hotspots |
| **Stage 3** (full panel) | Remaining 8 mutations | 640 ns | Complete correlation dataset |

### 4.5 Statistical Framework

| Component | Estimate | Justification |
|-----------|----------|---------------|
| σ_exp | 0.3 kcal/mol | Precise Ki measurements with replicates (Castro & Anderson 1996) |
| σ_comp | 0.5 kcal/mol | FEP statistical uncertainty from lambda windows |
| σ_method | 1.0 kcal/mol | Alchemical FEP for point mutations (§25.4: ±0.5–1.5) |
| σ_combined | 1.2 kcal/mol | √(0.3² + 0.5² + 1.0²) |
| 95% CI per mutant | ΔΔG_exp ± 2.3 kcal/mol | 1.96 × 1.2 |

### 4.6 Acceptance Criteria

| Metric | PASS | MARGINAL | FAIL |
|--------|------|----------|------|
| Spearman ρ (ΔΔG panel) | > 0.7 | 0.5–0.7 | < 0.5 |
| RMSE (kcal/mol) | < 2.0 | 2.0–3.0 | > 3.0 |
| K15A ΔΔG | > +7.0 | +4.0 to +7.0 | < +4.0 |
| G36A |ΔΔG| | < 1.5 | 1.5–3.0 | > 3.0 |
| Per-mutant within CI | ≥ 80% | 60–80% | < 60% |

---

## 5. Control Conditions

### 5.1 Positive Control (Internal)

**K15A (P1 → Ala):** The most drastic mutation, ΔΔG ≈ +10 kcal/mol. If the pipeline cannot detect this massive effect, it fails fundamentally.

### 5.2 Negative Control (Internal)

**G36A:** No experimental effect on binding. The pipeline should predict |ΔΔG| < 1.5 kcal/mol.

### 5.3 Cross-Experiment Control

**EXP-31 (disulfide ablation):** C14S/C38S double mutation via topology modification should be consistent with the C14A + C38A single mutation ΔΔG values from this experiment.

### 5.4 Convergence Control

Forward vs. reverse FEP ΔG should agree within 1.0 kcal/mol per mutation (hysteresis check).

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| ΔΔG (K15A) | ≈ +10 kcal/mol | Castro & Anderson 1996, Table 2 |
| ΔΔG (G36A) | ≈ 0 kcal/mol | Castro & Anderson 1996 |
| ΔΔG (C14A) | ≈ +7 kcal/mol | Disulfide disruption |
| ΔΔG (R17A) | ≈ +5 kcal/mol | Interface contact loss |
| Spearman ρ (15-mutant panel) | > 0.7 | Pipeline validation target |
| RMSE | < 2.0 kcal/mol | Typical for alchemical FEP |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Large K15 side chain** | Alchemical K→A transformation involves large charge + volume change | Soft-core LJ may not fully converge | High |
| **Disulfide mutations (C14A, C38A)** | Removing Cys from S-S pair requires topology change in addition to FEP | Hybrid topology complexity | High |
| **Sampling insufficiency** | ΔΔG estimates oscillate with more λ windows | 20 windows insufficient for large perturbations | Medium |
| **Enthalpy-entropy compensation missed** | ΔΔG correct but ΔΔH and TΔΔS individually wrong | Single-temperature FEP only gives ΔΔG | Low |
| **Force field bias** | Systematic over/under-estimation of alanine mutation effects | AMBER ff14SB parameterization | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| Structure quality | 2PTC from EXP-04 loads correctly | All atoms present, disulfides intact |
| Pilot K15A | ΔΔG > +5 kcal/mol | Correct direction and magnitude |
| Pilot G36A | |ΔΔG| < 2.0 kcal/mol | Near-null effect detected |
| Hysteresis (per mutation) | Forward-reverse ΔG agreement < 1.0 kcal/mol | Converged sampling |
| Histogram overlap | Adjacent λ windows overlap > 10% | Adequate phase space sampling |
| Stage 1 correlation (3 mutations) | Rank order: K15A > R17A > G36A | Correct qualitative ordering |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (BPTI-trypsin complex: ~35,000 atoms; FEP lambda windows with alchemical transformations) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation (15 mutant topologies) | 0.5 | 0.5 | CPU-bound; build alchemical systems |
| Equilibration (WT + 15 mutants, bound + free) | 5 | 2.5 | ~32 systems × 1.5 ns each |
| Stage 1: K15A, G36A, R17A (3 × 80 ns FEP) | 15 | 7.5 | 240 ns total (complex + solvent) |
| Stage 2: C14A, C38A, I18A, I19A (4 × 80 ns FEP) | 20 | 10 | 320 ns total |
| Stage 3: Remaining 8 mutations (8 × 80 ns FEP) | 45 | 22 | 640 ns total |
| BAR/MBAR analysis + correlation | 1 | 1 | CPU-bound |
| **Total** | **~90–100** | **~45–50** | §10.17: Tier 5 ≈ 90–100 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each λ window per mutation per leg |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-30/` |
| Restart procedure | Mount Drive → read `mutation_progress.json` → identify next incomplete mutation/λ → verify energy drift < 0.1% → resume |
| Estimated sessions needed | 5–6 (A100) or 3–4 (H100) |

**Session planning (staged approach):**

| Session | Content | Est. hours |
|---------|---------|------------|
| 1 | Equilibration of all systems + Stage 1 pilot (K15A, G36A, R17A) | 20 |
| 2 | Stage 2 (C14A, C38A, I18A, I19A) | 20 |
| 3–4 | Stage 3 part 1 (4 of 8 remaining) | 20–24 |
| 5–6 | Stage 3 part 2 (4 of 8 remaining) + analysis | 20–24 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| FEP window state | OpenMM binary `.chk` | `checkpoint_<mutation>_<leg>_lambda<idx>.chk` |
| Full simulation state | OpenMM XML | `state_<mutation>_<leg>_lambda<idx>.xml` |
| Per-window ΔU samples | NumPy `.npy` | `du_<mutation>_<leg>_lambda<idx>.npy` |
| Per-mutation ΔΔG results | JSON | `ddg_<mutation>.json` |
| Mutation completion log | JSON | `mutation_progress.json` |
| Correlation results | JSON | `alanine_scan_correlation.json` |

**Resume verification protocol:**
1. Read `mutation_progress.json` to determine next incomplete mutation.
2. Reload checkpoint: `simulation.loadCheckpoint('checkpoint_<mutation>_<leg>_lambda<idx>.chk')`.
3. Run 1000 steps; compute energy.
4. Compare to energy at checkpoint save: drift must be < 0.1%.
5. If drift exceeds threshold, discard window and re-equilibrate at that λ for 200 ps.

### 9.5 Platform Configuration

```python
from openmm import Platform
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
# For alanine scan FEP, automate across the mutation panel:
# for mutation in mutation_panel:
#     for leg in ['complex', 'solvent']:
#         for lam in np.linspace(0, 1, 20):
#             run_fep_window(mutation, leg, lam)
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
| Equilibrated BPTI-trypsin complex (2PTC) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/equilibrated_complex.pdb` | Starting structure for all 15 mutations |
| System XML | `/content/drive/MyDrive/v3_gpu_results/EXP-04/system/system.xml` | Force field parameters for topology modification |

**Pre-execution check:** EXP-04 must be completed. All 15 mutations are applied to the EXP-04 equilibrated BPTI structure.

**Cross-validation with EXP-31:** C14A and C38A ΔΔG values should be consistent with the C14S/C38S disulfide ablation result from EXP-31.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.0 — Created with GPU execution requirements for Step 5A Colab execution (2026-07-17).
