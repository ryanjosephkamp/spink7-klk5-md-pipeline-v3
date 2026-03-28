# EXP-08: Interfacial Hydrogen Bond Energy — Implementation Guide

**Experiment ID:** EXP-08  
**Feature ID:** F-08 (benchmarks.md)  
**Category:** Thermodynamic (Semi-Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Quantifies average energetic contribution per interfacial hydrogen bond in protease-inhibitor complexes. Literature: ~1.5 kcal/mol per H-bond (Krowarsch et al. 2003), 95% CI [0.7, 2.3] kcal/mol. Uses energy decomposition on BPTI-trypsin (PDB 2PTC) comparing bound vs. unbound states. Validates pipeline's electrostatic description of binding interface.

### 2. Hypothesis

**H₁:** Average per-H-bond energy = 1.5 ± 0.8 kcal/mol (95% CI [0.7, 2.3]). **H₂:** Total H-bond contribution correlates with persistent count (EXP-17), predicting ~15 kcal/mol from ~10 H-bonds.

### 3. Protocol

Use equilibrated BPTI-trypsin from EXP-04. Identify persistent interfacial H-bonds (>50% occupancy, D-A ≤ 3.5 Å, angle ≥ 135°). Method A: direct per-H-bond energy decomposition. Method B: ΔG/n_hbond correlation. PASS: [0.7, 2.3]; MARGINAL: [0.3, 3.5]; FAIL: outside.

### 4. Controls

Positive: P1 H-bonds (~2 kcal/mol each). Negative: solvent-exposed H-bonds should show ~0 net.

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-08_interfacial_hbond_energy"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import KCAL_TO_KJ
from src.analyze.structural import compute_rmsd
from src.analyze.contacts import compute_hbonds
from src.analyze.trajectory import load_trajectory
from src.physics.units import kj_to_kcal
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load EXP-04 Production Trajectory

```python
EXP04_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind" / "outputs"
solvated_pdb = EXP04_DIR / "solvated_complex.pdb"
prod_traj_path = EXP04_DIR / "production" / "production_trajectory.dcd"

traj = load_trajectory(prod_traj_path, solvated_pdb, stride=10)
n_frames = traj.n_frames
# Last 50 ns
frames_50ns = int(50.0 / 0.1)
traj_analysis = traj[-frames_50ns:]
print(f"Analysis: {traj_analysis.n_frames} frames")
```

### Step 3: Identify Interfacial Hydrogen Bonds

```python
topology = traj_analysis.topology
chain_a_sel = "chainid 0"  # trypsin
chain_b_sel = "chainid 1"  # BPTI

# Use pipeline's compute_hbonds function
hbond_results = compute_hbonds(traj_analysis, chain_a_sel, chain_b_sel)

# Also use mdtraj's baker_hubbard for robust identification
hbonds_bh = md.baker_hubbard(traj_analysis, freq=0.5)
print(f"Total persistent H-bonds (Baker-Hubbard, >50% occupancy): {len(hbonds_bh)}")

# Filter to interfacial only
bpti_atoms = set(topology.select(chain_b_sel).tolist())
trypsin_atoms = set(topology.select(chain_a_sel).tolist())

interfacial_hbonds = []
for hb in hbonds_bh:
    d, h, a = hb
    d_chain = d in bpti_atoms
    a_chain = a in bpti_atoms
    if d_chain != a_chain:  # one in BPTI, one in trypsin
        interfacial_hbonds.append(hb)

print(f"Interfacial H-bonds: {len(interfacial_hbonds)}")

# Compute per-H-bond occupancy
hbond_occupancies = {}
for hb in interfacial_hbonds:
    d, h, a = hb
    d_name = f"{topology.atom(d).residue.name}{topology.atom(d).residue.resSeq}-{topology.atom(d).name}"
    a_name = f"{topology.atom(a).residue.name}{topology.atom(a).residue.resSeq}-{topology.atom(a).name}"
    key = f"{d_name}...{a_name}"
    # Compute occupancy
    distances = md.compute_distances(traj_analysis, [[d, a]])
    occupancy = np.mean(distances < 0.35)
    hbond_occupancies[key] = occupancy
    print(f"  {key}: occupancy = {occupancy:.2f}")
```

### Step 4: Per-H-Bond Energy Estimation

```python
# Method A: Electrostatic interaction energy between D-A pairs
from openmm import XmlSerializer
system_xml_path = EXP04_DIR / "equilibration" / "system.xml"

per_hbond_energies = []
for hb in interfacial_hbonds:
    d, h, a = hb
    # Compute time-averaged D...A distance
    distances = md.compute_distances(traj_analysis, [[d, a]])[0]
    # Approximate Coulombic energy: E ≈ -332 * q_d * q_a / (ε * r) kcal/mol
    # For backbone N-H...O=C: q_donor ≈ +0.27, q_acceptor ≈ -0.57
    # ε_eff ≈ 4 for protein interior
    mean_dist_ang = np.mean(distances) * 10  # nm to Å
    e_hbond_kcal = -332 * 0.27 * 0.57 / (4.0 * mean_dist_ang)

    # Subtract desolvation penalty (~0.5 kcal/mol per D or A)
    net_e = abs(e_hbond_kcal) - 1.0  # rough desolvation correction
    per_hbond_energies.append(max(0.1, net_e))

avg_per_hbond = np.mean(per_hbond_energies)
std_per_hbond = np.std(per_hbond_energies) / np.sqrt(len(per_hbond_energies))
total_hbond_energy = np.sum(per_hbond_energies)

print(f"\nAverage per-H-bond energy: {avg_per_hbond:.2f} ± {std_per_hbond:.2f} kcal/mol")
print(f"Total H-bond contribution: {total_hbond_energy:.1f} kcal/mol ({len(interfacial_hbonds)} bonds)")

# Method B: ΔG from EXP-04 / n_hbonds
exp04_results_path = PROJECT_ROOT / "v3_experiments" / "EXP-04_bpti_trypsin_dg_bind" / "results.json"
if exp04_results_path.exists():
    with open(exp04_results_path) as f:
        exp04 = json.load(f)
    dg_04 = exp04.get("dg_bind_wham_kcal", -18.0)
    method_b = abs(dg_04) / len(interfacial_hbonds) if len(interfacial_hbonds) > 0 else 0
    print(f"Method B (ΔG/n_hbond): {method_b:.2f} kcal/mol")
```

### Step 5: Classification

```python
dg_exp_per_hbond = 1.5
if 0.7 <= avg_per_hbond <= 2.3:
    classification = "PASS"
elif 0.3 <= avg_per_hbond <= 3.5:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-08", "feature_id": "F-08",
    "avg_per_hbond_kcal": float(avg_per_hbond),
    "std_per_hbond_kcal": float(std_per_hbond),
    "n_interfacial_hbonds": len(interfacial_hbonds),
    "total_hbond_contribution_kcal": float(total_hbond_energy),
    "exp_per_hbond_kcal": 1.5,
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-08 Classification: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Per-H-Bond Energy Distribution

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.bar(range(len(per_hbond_energies)), sorted(per_hbond_energies, reverse=True),
       color='steelblue', edgecolor='black')
ax.axhline(y=1.5, color='red', linestyle='--', linewidth=2, label='Literature (1.5 kcal/mol)')
ax.axhspan(0.7, 2.3, alpha=0.1, color='green', label='95% CI [0.7, 2.3]')
ax.set_xlabel('H-bond index (ranked)', fontsize=12)
ax.set_ylabel('Energy (kcal/mol)', fontsize=12)
ax.set_title('EXP-08: Per-H-Bond Energetic Contributions', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-08_per_hbond_energy.png", dpi=300)
plt.close(fig)
```

### Figure 2: H-Bond Occupancy

```python
fig, ax = plt.subplots(1, 1, figsize=(14, 6))
names = list(hbond_occupancies.keys())
occs = list(hbond_occupancies.values())
ax.barh(names, occs, color='coral', edgecolor='black')
ax.axvline(x=0.5, color='gray', linestyle='--', label='50% threshold')
ax.set_xlabel('Occupancy', fontsize=12)
ax.set_title('EXP-08: Interfacial H-Bond Occupancy', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-08_hbond_occupancy.png", dpi=300)
plt.close(fig)
```

### Figure 3: Predicted vs Expected Per-H-Bond Energy

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.bar(['Literature', 'Pipeline (Method A)'], [1.5, avg_per_hbond],
       color=['gold', 'steelblue'], edgecolor='black')
ax.errorbar([1], [avg_per_hbond], yerr=[std_per_hbond*1.96], fmt='none', color='black', capsize=5)
ax.axhspan(0.7, 2.3, alpha=0.1, color='green', label='95% CI')
ax.set_ylabel('Per-H-bond Energy (kcal/mol)', fontsize=12)
ax.set_title(f'EXP-08: H-Bond Energy — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-08_predicted_vs_expected.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-08: Interfacial Hydrogen Bond Energy — Results Report

**Experiment ID:** EXP-08  **Feature ID:** F-08  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## 1. Abstract
## 2. Introduction
## 3. Hypothesis
## 4. Methods
## 5. Controls
## 6. Results
### 6.1 Per-H-bond energy: [value] ± [σ] kcal/mol
### 6.2 N interfacial H-bonds: [value]
### 6.3 Total H-bond contribution: [value] kcal/mol
### 6.4 Method A vs B consistency
## 7. Discussion
## 8. Conclusions
## 9. Figures
## 10. References
## 11. Author Block
---
Author: Ryan Kamp
Affiliation: Dept. of Computer Science, University of Cincinnati
Contact:
Email: kamprj@mail.uc.edu
GitHub: ryanjosephkamp
```

---

## Part 5 — GPU/Colab Execution Procedures

> **Added:** Step 5A GPU documentation update. These sections augment the existing Part 2 steps. EXP-08 loads the EXP-04 production trajectory — GPU accelerates H-bond energy estimation.

### §5.1 Colab Environment Setup

```python
# §5.1 Colab Environment Setup
!nvidia-smi

from google.colab import drive
drive.mount('/content/drive')

!pip install openmm mdtraj parmed matplotlib numpy scipy pymbar openmmtools

import os, sys
from pathlib import Path

EXP_ID = "EXP-08"
DRIVE_BASE = Path(f"/content/drive/MyDrive/v3_gpu_results/{EXP_ID}")
for subdir in ["checkpoints", "outputs", "figures"]:
    (DRIVE_BASE / subdir).mkdir(parents=True, exist_ok=True)

!ln -sf /content/drive/MyDrive/medium_project_2/src /content/src
sys.path.insert(0, "/content")

EXP04_DRIVE = Path("/content/drive/MyDrive/v3_gpu_results/EXP-04")
assert (EXP04_DRIVE / "outputs" / "production").exists(), \
    "EXP-04 production outputs required — run EXP-04 first"

PROJECT_ROOT = Path("/content")
EXP_DIR = DRIVE_BASE
OUTPUT_DIR = DRIVE_BASE / "outputs"
FIGURES_DIR = DRIVE_BASE / "figures"

print(f"Colab environment ready. EXP-04 data: {EXP04_DRIVE}")
```

### §5.2 GPU Platform Selection

```python
# §5.2 GPU Platform Selection
import openmm
from src.simulate.platform import select_platform

platform = select_platform("CUDA")
properties = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
print(f"Platform: {platform.getName()}")

# EXP-08 uses GPU for H-bond energy re-evaluation across trajectory frames.
```

### §5.3 Checkpoint and Resume Integration

```python
# §5.3 Checkpoint and Resume Integration
import json, time
import numpy as np

CHECKPOINT_DIR = DRIVE_BASE / "checkpoints"

def save_analysis_checkpoint(data, phase_name):
    """Save analysis progress to Drive."""
    ckpt = CHECKPOINT_DIR / f"{phase_name}.npz"
    np.savez(ckpt, **{k: v for k, v in data.items() if isinstance(v, np.ndarray)})
    meta = {k: v for k, v in data.items() if not isinstance(v, np.ndarray)}
    meta["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(CHECKPOINT_DIR / f"{phase_name}_meta.json", 'w') as f:
        json.dump(meta, f, indent=2, default=str)
    print(f"  Analysis checkpoint saved: {phase_name}")

# ─── Checkpoint points for EXP-08 ───
# After H-bond identification: save_analysis_checkpoint(hbonds, "hbond_id")
# After energy estimation: save_analysis_checkpoint(energies, "hbond_energy")
```

### §5.4 Progress Monitoring

```python
# §5.4 Progress Monitoring
import time, subprocess

def report_gpu_status():
    result = subprocess.run(
        ['nvidia-smi', '--query-gpu=name,memory.used,memory.total,utilization.gpu',
         '--format=csv,noheader'], capture_output=True, text=True)
    print(f"GPU: {result.stdout.strip()}")

# ─── EXP-08 runtime estimates (A100 40GB) ───
# H-bond identification + energy: ~30–60 minutes
# Total: ~1 hour
```

### §5.5 Results Persistence

```python
# §5.5 Results Persistence
import shutil

def sync_to_drive(local_dir, drive_subdir="outputs"):
    drive_dir = DRIVE_BASE / drive_subdir
    drive_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for f in Path(local_dir).rglob("*"):
        if f.is_file():
            dest = drive_dir / f.relative_to(local_dir)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dest)
            count += 1
    print(f"  Synced {count} files → {drive_dir}")

def verify_drive_sync():
    critical = ["results.json"]
    for f in critical:
        path = DRIVE_BASE / f
        status = "✓" if path.exists() else "✗ MISSING"
        print(f"  {status} {f}")
```

### §5.6 GPU-Optimized Figure Generation

```python
# §5.6 GPU-Optimized Figure Generation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save_figure(fig, name):
    fig_path = DRIVE_BASE / "figures" / f"{name}.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(fig_path), dpi=150, bbox_inches='tight')
    print(f"  Figure saved: {fig_path}")
    plt.close(fig)
```

### §5.7 Error Recovery Procedures

```python
# §5.7 Error Recovery Procedures
# ─── EXP-08–specific recovery ───
# 1. EXP-04 trajectory not found: verify EXP-04 completed on Drive
# 2. H-bond occupancy calculation failure:
#    - Ensure mdtraj.baker_hubbard() receives aligned trajectory
#    - Check donor-acceptor distance/angle thresholds
# 3. Method A vs B inconsistency:
#    - Not a GPU error — document both estimates
```

---

Revision: v1.1 — Added GPU/Colab execution sections (Part 5, §5.1–§5.7) for Step 5A.

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
