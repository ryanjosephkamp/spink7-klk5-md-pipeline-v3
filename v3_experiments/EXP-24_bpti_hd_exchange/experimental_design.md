# EXP-24: BPTI Hydrogen/Deuterium Exchange – Amide Protection Pattern

**Experiment ID:** EXP-24  
**Feature ID:** F-24 (benchmarks.md)  
**Category:** Dynamic  
**Status:** SEMI-QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment validates the pipeline's ability to reproduce the hydrogen/deuterium (H/D) exchange protection pattern of BPTI. Experimental NMR studies identify 11 backbone amide protons that are highly protected from exchange (slow exchange, protection factor >10⁴), corresponding to residues buried in stable secondary structures or shielded by disulfide bonds (Dempsey 2001, Roder & Wüthrich 1986). The pipeline should identify these protected amides from MD-derived metrics: low RMSF, high H-bond persistence, and low solvent accessibility.

---

## 2. Hypothesis

**H₁:** MD analysis will correctly identify the 11 highly protected amide protons (protection factor >10⁴) as those with: (a) backbone NH solvent accessibility < 10%, AND (b) persistent intramolecular H-bond occupancy > 80%.

**H₂:** The correlation between MD-derived protection proxies and experimental protection factors will be r > 0.5.

---

## 3. Background and Rationale

H/D exchange rates report on local structural stability and solvent accessibility. Amides buried in stable secondary structures exchange slowly (days to years), while exposed loop amides exchange within minutes. The 11 protected amides in BPTI are among the best-characterized H/D exchange data for any protein, making them an ideal benchmark for MD-predicted dynamics. Matching the protection pattern validates the force field's representation of local flexibility, H-bonding, and solvent penetration.

---

## 4. Experimental Protocol

### 4.1 System

Free BPTI (PDB 4PTI) from EXP-23 production trajectory.

### 4.2 Protected Amide Identification from Experiment

The 11 highly protected amides in BPTI (from NMR H/D exchange):
- Core β-sheet and α-helix residues: ~Y21, F22, Y23, N24, A25 (β-sheet); ~L49, F50, A51, C51, D53 (α-helix); and C14 (disulfide-stabilized).
- Exact residue list will be confirmed from Dempsey 2001.

### 4.3 MD-Derived Protection Proxies

For each backbone amide N-H in BPTI (58 residues minus Pro and N-term = ~52 amides):

**4.3.1 Solvent Accessibility:**
1. Per-frame SASA of each backbone NH hydrogen.
2. Time-averaged: low SASA (<10% of max) → protected.

**4.3.2 Intramolecular H-Bond Occupancy:**
1. For each backbone N-H, check if it is H-bonded to any backbone C=O.
2. H-bond criterion: d(N···O) < 3.5 Å, angle > 135°.
3. Occupancy = fraction of frames with H-bond.
4. High occupancy (>80%) → protected.

**4.3.3 RMSF:**
1. Per-residue Cα RMSF from MD.
2. Low RMSF (< 0.5 Å) → rigid, likely protected.

**4.3.4 Combined Protection Score:**
Protection_proxy = f(SASA, H-bond occupancy, RMSF) — a composite metric.
Binary classification: protected if SASA < 10% AND HB occupancy > 80%.

### 4.4 Production MD Parameters (from config.py)

- Force field: AMBER ff14SB
- Duration: 100 ns
- Temperature: 310 K
- Save interval: 10 ps (10,000 frames)

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Correct Identification of 11 Protected Amides |
|---------------|-----------------------------------------------|
| **PASS** | ≥8 of 11 correctly identified; ≤3 false positives |
| **MARGINAL** | 6–7 of 11 correctly identified; ≤5 false positives |
| **FAIL** | <6 of 11 or >5 false positives |

---

## 5. Control Conditions

### 5.1 Positive Control: Crystal B-Factors

Crystallographic B-factors should correlate with MD RMSF (r > 0.5). Low B-factor residues should overlap with protected amides.

### 5.2 Negative Control: Terminal Residues

N-terminal and C-terminal residues should be identified as unprotected (fast exchange). If the proxy incorrectly classifies them as protected, the metric is flawed.

### 5.3 Disulfide-Adjacent Control

Residues adjacent to C5-C55, C14-C38, C30-C51 should show moderate protection due to disulfide stabilization, even if not in regular secondary structure.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Protected amides (experimental) | 11 | Dempsey 2001 |
| True positives (from MD proxy) | ≥8 | benchmarks.md |
| False positives | ≤3 | Specificity |
| Protection-RMSF anticorrelation | r < −0.5 | Expected relationship |
| H-bond occupancy for protected | >80% | By definition |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Proxy too crude** | Poor correlation with experimental PF | Linear proxy insufficient | Medium |
| **TIP3P over-penetration** | Over-estimates solvent access → false unprotected | Water model artifact | Medium |
| **100 ns too short** | Miss rare opening events | H/D exchange occurs on ms-h timescale | Medium |
| **Missing slow motions** | Underestimate flexibility → false protected | Insufficient sampling | Medium |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| RMSF profile | Core rigid, termini flexible | Matches B-factor pattern |
| H-bond persistence | Key β-sheet/α-helix H-bonds > 80% | Correct secondary structure |
| SASA analysis | Buried amides have SASA ≈ 0 | Physical expectation |
| Binary classification | ≥8/11 true positives | Above chance (~3/11) |
| False positive check | Terminal residues unprotected | Correct |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~2–4 GB (free BPTI, PDB 4PTI: ~15,000 atoms solvated in TIP3P with 1.2 nm padding) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Structure preparation | < 0.1 | < 0.1 | CPU-bound |
| Equilibration (NVT + NPT) | 0.3 | 0.2 | 1.5 ns total; small system |
| Production MD (100 ns) | 3 | 1.5 | ~2000 ns/day (A100) for ~15k atoms |
| H/D exchange analysis (protection factors) | 1 | 1 | CPU-bound analysis over trajectory |
| Solvent accessibility + H-bond analysis | 0.5 | 0.5 | Per-frame backbone amide analysis |
| **Total** | **~8** | **~4** | §10.17: Tier 3 ≈ 8 GPU-hrs (A100) |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | Every 10 ns of production MD; after equilibration |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-24/` |
| Restart procedure | Mount Drive → load latest checkpoint → verify energy drift < 0.1% → resume production MD |
| Estimated sessions needed | 1 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Positions, velocities, box vectors, RNG state | OpenMM binary `.chk` | `checkpoint_production_<ns>.chk` |
| Full simulation state (portable) | OpenMM XML serialization | `state_production_<ns>.xml` |
| Trajectory (DCD) | MDTraj `.dcd` | `bpti_production_100ns.dcd` |
| Protection factor results | JSON | `protection_factors.json` |

**Resume verification protocol:**
1. Reload checkpoint: `simulation.loadCheckpoint('checkpoint_production_<ns>.chk')`.
2. Run 1000 steps; compute energy.
3. Compare to energy at checkpoint save: drift must be < 0.1%.
4. If drift exceeds threshold, discard last 100 ps and re-equilibrate for 200 ps.

### 9.5 Platform Configuration

```python
from openmm import Platform
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
```

**Environment verification:**
```bash
!nvidia-smi
python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"
```

### 9.6 Dependency Notes

EXP-24 has no upstream dependencies. EXP-24 outputs are required by:

| Downstream Experiment | Required Output | Drive Path |
|-----------------------|----------------|------------|
| EXP-25 (Conformational variability) | 100 ns production MD trajectory | `/content/drive/MyDrive/v3_gpu_results/EXP-24/production/bpti_production_100ns.dcd` |

The EXP-24 trajectory file must be persisted to Google Drive immediately after production MD completes.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
