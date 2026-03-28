# EXP-33: Barnase–Barstar Double-Mutant Cycle Analysis (~45 Pairs, FEP)

**Experiment ID:** EXP-33  
**Feature ID:** F-33 (benchmarks.md)  
**Category:** Mutational  
**Status:** QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates the pipeline's FEP methodology on the barnase-barstar double-mutant cycle dataset, the largest and most rigorous experimental dataset for protein-protein binding energetics. Schreiber & Fersht (1995) measured ΔΔG for ~45 double-mutant pairs, enabling quantification of coupling energies (ΔΔΔ G) between pairs of interface residues. The pipeline should reproduce the single-mutation ΔΔG values with R > 0.7 and correctly identify coupled pairs (|ΔΔΔG| > 1.5 kcal/mol). This is the ultimate test of computational mutagenesis accuracy.

---

## 2. Hypothesis

**H₁:** FEP-calculated single-mutation ΔΔG values for barnase-barstar will correlate with experimental values with R > 0.7.

**H₂:** FEP-calculated coupling energies (ΔΔΔG) will correctly identify the top 5 most coupled residue pairs.

**H₃:** RMSE for single-mutation ΔΔG will be < 2.0 kcal/mol.

---

## 3. Background and Rationale

Double-mutant cycle analysis is the gold standard for measuring pairwise interaction energies between residues at protein-protein interfaces. For a pair of mutations A→A' and B→B':

$$\Delta\Delta\Delta G = \Delta\Delta G_{AB \to A'B'} - \Delta\Delta G_{A \to A'} - \Delta\Delta G_{B \to B'}$$

A non-zero ΔΔΔG indicates energetic coupling between positions A and B (they interact directly or through a shared structural element). Barnase-barstar has ~10 interfacial mutations characterized in all combinations, yielding ~45 double-mutant cycles. This provides a high-dimensional test of FEP accuracy.

---

## 4. Experimental Protocol

### 4.1 System Preparation

- PDB: 1BRS (barnase-barstar complex, 2.0 Å resolution)
- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P (`amber14/tip3p.xml`)
- pH: 7.4
- Box padding: 1.2 nm, ionic strength 0.15 M NaCl
- Equilibration: standard protocol (10,000 step min, 500 ps NVT, 1000 ps NPT at 310 K)

### 4.2 Mutation Panel

Key barnase mutations (Schreiber & Fersht 1995):
- Barnase: K27A, R59A, H102A, R83A, R87A, D54A, E73A, and others
- Barstar: D35A, D39A, Y29A, E76A, T42A, G43A, and others
- Total single mutations: ~10–15
- Double-mutant combinations: ~45

### 4.3 FEP Protocol (per mutation)

For each X→A mutation:

**Leg 1: Complex**
1. Hybrid topology: X → Ala
2. Lambda windows: n_lambda = 20
3. Per-window: 2.0 ns
4. Temperature: 310 K
5. Soft-core: alpha = 0.5, power = 1
6. Save: 1.0 ps
7. n_equil: 5000

**Leg 2: Free protein (barnase or barstar alone)**
1. Same protocol on the free chain

$$\Delta\Delta G = \Delta G_{complex}^{X \to A} - \Delta G_{free}^{X \to A}$$

### 4.4 Double-Mutant Cycle Computation

For each pair (X_i, X_j):
1. Compute ΔΔG(X_i → A) — single mutation i
2. Compute ΔΔG(X_j → A) — single mutation j
3. Compute ΔΔG(X_i,X_j → A,A) — double mutation (both simultaneously)
4. ΔΔΔG = ΔΔG(double) − ΔΔG(single_i) − ΔΔG(single_j)
5. |ΔΔΔG| > 1.5 kcal/mol → coupled pair

### 4.5 Staged Implementation

Due to the large number of FEP calculations (~45 double + ~15 single = ~60 mutations × 2 legs × 20 λ × 2 ns = ~4800 ns), this experiment may be staged:

**Stage 1:** 10 single mutations (10 × 2 × 20 × 2 = 800 ns)
**Stage 2:** 5 double mutations (most coupled pairs from experiment) (5 × 2 × 20 × 2 = 400 ns)
**Stage 3:** Remaining doubles if Stage 1-2 pass (35 × 2 × 20 × 2 = 2800 ns)

### 4.6 Acceptance Criteria (from benchmarks.md)

| Classification | Criterion |
|---------------|-----------|
| **PASS** | R > 0.7, RMSE < 2.0 kcal/mol, top 5 coupled pairs correctly identified |
| **MARGINAL** | R > 0.5, RMSE < 3.0 kcal/mol, top 3 coupled pairs correct |
| **FAIL** | R < 0.5 or RMSE > 3.0 or coupled pairs not identified |

---

## 5. Control Conditions

### 5.1 Self-Perturbation Control

A→A identity transformation on 2 residues: ΔΔG should be 0.0 ± 0.3 kcal/mol.

### 5.2 Additivity Control

For non-interacting pairs (residues on opposite sides of the interface), ΔΔΔG should be ~0 (mutations are independent). This validates the cycle thermodynamics.

### 5.3 Cross-System Validation

Compare barnase-barstar FEP accuracy to BPTI-trypsin FEP accuracy (EXP-30). Both systems should show similar R and RMSE values.

### 5.4 EXP-13 Cross-Reference

The barnase-barstar ΔG from US/WHAM (EXP-13 / structural validation) should be consistent with the FEP-derived stability.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| R (single mutations) | >0.7 | benchmarks.md |
| RMSE (single mutations) | < 2.0 kcal/mol | FEP standard |
| Top coupled pairs | R59-D35, K27-D39, H102-D35 | Schreiber & Fersht 1995 |
| ΔΔΔG for coupled pairs | 1.5–5.0 kcal/mol | Direct interaction |
| ΔΔΔG for non-coupled | < 0.5 kcal/mol | Independent sites |
| Total compute (full) | ~4800 ns | Resource estimate |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Charge-change FEP errors** | K27A, R59A, D35A all involve charge changes | Net charge ΔQ ≠ 0 | Critical |
| **Computational cost** | 4800 ns exceeds available resources | Need staged approach | High |
| **Error accumulation in DMC** | ΔΔΔG computed from 3 independent FEP → error compounds | √3× individual error | High |
| **Barnase flexibility** | Active-site loop more flexible than BPTI | Different dynamics | Medium |
| **Double mutation FEP** | Simultaneous large perturbation → poor overlap | Need more λ windows | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| 1BRS quality | Structure loads, all interface atoms present | No missing residues |
| Equilibration | Complex RMSD < 2.5 Å | Stable |
| Self-FEP | A→A ΔΔG ≈ 0 | |ΔΔG| < 0.3 kcal/mol |
| Stage 1 (singles) | 10 ΔΔG values computed | R > 0.5 preliminary |
| Lambda overlap | All windows show adequate overlap | MBAR converges |
| Top coupling | R59-D35 ΔΔΔG > 1.5 kcal/mol | Key pair identified |
| Non-coupled control | Two peripheral mutations ΔΔΔG ≈ 0 | Independence confirmed |
| Cross-system | Barnase FEP accuracy ≈ BPTI FEP accuracy | Robust methodology |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~4–6 GB (barnase-barstar complex: ~50,000 atoms; FEP with ~45 double-mutant cycles) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation (~10 single mutations) | 0.5 | 0.5 | CPU-bound; build mutant topologies |
| Equilibration (WT + mutants, bound + free) | 5 | 2.5 | ~20 systems × 1.5 ns each |
| FEP pilot mutations (5 singles, 20λ × 2 ns each) | 15 | 7.5 | 200 ns total (complex + solvent) |
| FEP remaining singles (5 more, 20λ × 2 ns each) | 15 | 7.5 | 200 ns total |
| Double-mutant FEP (~15 key pairs, 20λ × 2 ns each) | 60 | 30 | ~600 ns total; staged approach |
| BAR/MBAR analysis + DMC ΔΔG_int calculation | 2 | 2 | CPU-bound |
| Correlation analysis (R, RMSE vs. Schreiber & Fersht) | 0.5 | 0.5 | CPU-bound |
| **Total** | **~100–120** | **~50–60** | §10.17: Tier 5 ≈ 100–120 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each λ window per mutation per leg; after each equilibration |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-33/` |
| Restart procedure | Mount Drive → identify last completed mutation/λ → verify energy drift < 0.1% → resume |
| Estimated sessions needed | 5–6 (A100) or 3–4 (H100) |

**Session planning (staged approach):**
| Session | Content | Est. hours |
|---------|---------|------------|
| 1 | Equilibration of all systems + pilot 2–3 single mutations | 20 |
| 2–3 | Remaining single mutations | 20–30 |
| 4–6 | Key double-mutant pairs + analysis | 40–60 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| FEP window state | OpenMM binary `.chk` | `checkpoint_<mutation>_<leg>_lambda<idx>.chk` |
| Full simulation state | OpenMM XML | `state_<mutation>_<leg>_lambda<idx>.xml` |
| Per-window ΔU samples | NumPy `.npy` | `du_<mutation>_<leg>_lambda<idx>.npy` |
| Per-mutation ΔΔG results | JSON | `ddg_<mutation>.json` |
| DMC interaction energies | CSV | `dmc_interaction_energies.csv` |
| Mutation completion log | JSON | `mutation_progress.json` |

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
# For DMC-FEP, automate lambda scheduling across mutations:
# for mutation in mutation_list:
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

**Upstream dependency on EXP-13:**

| Required EXP-13 Output | Drive Path | Description |
|-----------------------|------------|-------------|
| Equilibrated barnase-barstar complex (1BRS) | `/content/drive/MyDrive/v3_gpu_results/EXP-13/production/` | Starting structure for mutation panel |
| System XML | `/content/drive/MyDrive/v3_gpu_results/EXP-13/system/system.xml` | Force field parameters |

**Pre-execution check:** EXP-13 must be completed (equilibrated 1BRS structure on Drive). The ~45 double-mutant cycles use pairs from Schreiber & Fersht (1995) Table 2.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
