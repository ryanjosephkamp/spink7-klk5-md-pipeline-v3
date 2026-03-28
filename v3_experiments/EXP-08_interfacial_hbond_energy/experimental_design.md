# EXP-08: Interfacial Hydrogen Bond Energy

**Experiment ID:** EXP-08  
**Feature ID:** F-08 (benchmarks.md)  
**Category:** Thermodynamic  
**Status:** SEMI-QUANTITATIVE  
**Date:** 2026-03-22  

---

## 1. Abstract

This experiment quantifies the average energetic contribution per interfacial hydrogen bond in protease-inhibitor complexes. Literature consensus places this at ~1.5 kcal/mol per H-bond (Krowarsch et al. 2003), with a 95% pipeline confidence interval of [0.7, 2.3] kcal/mol. The experiment uses thermodynamic integration or energy decomposition on the BPTI-trypsin complex, comparing bound vs. unbound states to extract per-H-bond energetic contributions. This validates the pipeline's electrostatic description of the binding interface.

---

## 2. Hypothesis

**H₁:** The average per-hydrogen-bond energetic contribution to binding free energy in the BPTI-trypsin interface is 1.5 ± 0.8 kcal/mol (95% CI: [0.7, 2.3] kcal/mol).

**H₂:** The total H-bond energetic contribution correlates with the number of persistent interfacial H-bonds (EXP-17), predicting ~15 kcal/mol total from ~10 H-bonds.

---

## 3. Background and Rationale

Hydrogen bonds at the protease-inhibitor interface stabilize the Laskowski loop in the canonical binding orientation. Unlike bulk solvent H-bonds (compensated by water→water H-bonds upon desolvation), interfacial H-bonds between complementary donor-acceptor pairs provide net stabilization. The ~1.5 kcal/mol value reflects the difference between the protein-protein H-bond and the water-mediated alternative, consistent with mutagenesis data across multiple protease families.

---

## 4. Experimental Protocol

### 4.1 System Preparation

Use the equilibrated BPTI-trypsin complex from EXP-04 (PDB 2PTC).

### 4.2 Hydrogen Bond Identification

1. Identify all persistent interfacial H-bonds from the last 50 ns of production MD.
2. H-bond criteria: donor-acceptor distance ≤ 3.5 Å, D-H···A angle ≥ 135°.
3. Persistent = present in >50% of frames.
4. Expected count: ~10 persistent H-bonds (cross-reference EXP-17).

### 4.3 Per-H-Bond Energy Estimation

**Method A: Direct energy decomposition:**
1. For each persistent H-bond pair, compute time-averaged electrostatic + vdW interaction energy.
2. Subtract the equivalent energy of each donor/acceptor with water in the unbound state.
3. Net ΔE_HB = E_protein-protein − E_protein-water.

**Method B: Correlation approach:**
1. From the total ΔG_bind (EXP-04) and total persistent H-bond count (EXP-17), compute average per-H-bond contribution.
2. This is a consistency check, not an independent measurement.

### 4.4 Simulation Parameters

Production MD parameters (from config.py):
- Force field: AMBER ff14SB (`amber14-all.xml`)
- Water model: TIP3P (`amber14/tip3p.xml`)
- Temperature: 310 K (Langevin thermostat, friction 1.0 ps⁻¹)
- Timestep: 2 fs
- Trajectory save interval: 10 ps
- Analysis window: last 50 ns of 100 ns production run

### 4.5 Acceptance Criteria (from benchmarks.md)

| Classification | Per-H-Bond Energy |
|---------------|-------------------|
| **PASS** | [0.7, 2.3] kcal/mol |
| **MARGINAL** | [0.3, 3.5] kcal/mol |
| **FAIL** | Outside marginal range |

---

## 5. Control Conditions

### 5.1 Positive Control

- **BPTI C14-C38 H-bonds:** These backbone H-bonds are well-characterized in neutron diffraction studies and provide a structural anchor.
- **EXP-07 cross-check:** The P1 residue forms 2–3 strong H-bonds; their individual energies should be at the high end (~2 kcal/mol each).

### 5.2 Negative Control

- **Solvent-exposed H-bonds:** Residues not at the interface should show near-zero net H-bond stabilization (compensated by water).
- **H-bond to bulk water:** Protein surface→water H-bonds should average ~0 net contribution.

---

## 6. Expected Outcomes

| Metric | Expected Value | Source |
|--------|---------------|--------|
| Per-H-bond energy | ~1.5 kcal/mol | Krowarsch 2003 |
| 95% CI | [0.7, 2.3] kcal/mol | benchmarks.md |
| Total H-bond contribution | ~15 kcal/mol | ~10 bonds × 1.5 |
| Strongest individual H-bond | P1 backbone N-H···O=C | Structural consensus |

---

## 7. Potential Failure Modes

| Failure Mode | Manifestation | Limitation | Severity |
|-------------|--------------|-----------|----------|
| **Fixed-charge overestimation** | Per-H-bond energy >3 kcal/mol | No polarization in ff14SB | Medium |
| **Desolvation penalty underestimate** | Net energy too favorable | Implicit continuum effects missing | Medium |
| **H-bond geometry fluctuation** | Persistent count unstable | 3.5 Å cutoff sensitivity | Low |

---

## 8. Intermediate Verification Tests

| Step | Verification | Pass Criterion |
|------|-------------|----------------|
| H-bond identification | Persistent count matches EXP-17 | 8–12 persistent H-bonds |
| Energy decomposition | All per-H-bond values negative | Stabilizing interactions |
| Donor/acceptor balance | ~equal D and A contributions | Symmetric pairing |
| Method A vs B consistency | Per-H-bond estimates within 0.5 | Agreement between methods |

---

## 9. GPU Execution Requirements (Step 5A)

> **Added:** v1.1 — GPU experiment execution via Google Colab (Step 5A, Task 5).

### 9.1 GPU Hardware Requirements

| Requirement | Specification |
|-------------|---------------|
| Minimum GPU | NVIDIA A100 40 GB or H100 80 GB |
| CUDA version | ≥ 12.0 |
| OpenMM version | ≥ 8.1, with CUDA platform |
| Estimated VRAM | ~2–4 GB (H-bond energy analysis on pre-computed trajectory; no live simulation) |

### 9.2 Runtime Estimates

| Phase | A100 (hours) | H100 (hours) | Notes |
|-------|-------------|-------------|-------|
| Load EXP-04 trajectory (last 50 ns, 500 frames) | 0.1 | 0.1 | I/O from Drive |
| H-bond identification and energy estimation | 0.8 | 0.5 | Per-frame H-bond detection + energy |
| Bootstrap analysis per H-bond | 0.3 | 0.3 | CPU-bound statistical analysis |
| **Total** | **~1–2** | **~1** | §10.17: Tier 2 ≈ 1–2 GPU-hrs |

### 9.3 Colab Session Management

| Parameter | Value |
|-----------|-------|
| Maximum session duration | 24 hours |
| Checkpoint frequency | After each trajectory batch (every 100 frames) |
| Google Drive mount path | `/content/drive/MyDrive/v3_gpu_results/EXP-08/` |
| Restart procedure | Mount Drive → load partial H-bond results → resume from next unprocessed frame |
| Estimated sessions needed | 1 |

### 9.4 Checkpoint Strategy

| State Component | Format | Naming Convention |
|----------------|--------|-------------------|
| Per-frame H-bond occupancy (partial) | NumPy `.npy` | `hbond_occupancy_frames_<start>_<end>.npy` |
| Per-H-bond energy estimates | NumPy `.npz` | `hbond_energies.npz` |

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

**CRITICAL — Tier 2 upstream dependency on EXP-04:**

| Required EXP-04 Output | Drive Path | Description |
|-----------------------|------------|-------------|
| Production MD trajectory (last 50 ns) | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/trajectory_last50ns.dcd` | 500 frames at 100 ps intervals |
| Equilibrated complex structure | `/content/drive/MyDrive/v3_gpu_results/EXP-04/production/equilibrated_complex.pdb` | Topology file for H-bond geometry |
| System XML | `/content/drive/MyDrive/v3_gpu_results/EXP-04/system/system.xml` | Force field parameters |

**Pre-execution check:** Verify all EXP-04 outputs exist on Drive before starting EXP-08. If missing, EXP-04 must be completed first.

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp  
Revision: v1.1 — Added §9 (GPU Execution Requirements for Step 5A Colab execution).
