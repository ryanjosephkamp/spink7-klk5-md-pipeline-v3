# EXP-21: KLK5 Subsite Specificity — Implementation Guide

**Experiment ID:** EXP-21  
**Feature ID:** F-21 (benchmarks.md)  
**Category:** Structural (Qualitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Validates KLK5 subsite architecture (S1–S4, S1'–S3') by mapping contacts with inhibitor binding loop (P4–P3'). Trypsin-like S1 containing Asp189 for basic P1 selectivity. Uses EXP-01 (SPINK7–KLK5) and EXP-04 (BPTI–trypsin) trajectories for cross-comparison.

### 2. Classification (§25.1)

- PASS: S1 Asp salt bridge >80% occupancy; all subsites correctly mapped
- MARGINAL: S1 correct but 1–2 subsites ambiguous
- FAIL: S1 not formed or P1 not in S1

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-21_klk5_subsite_specificity"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.analyze.trajectory import load_trajectory
from src.analyze.structural import compute_interface_contacts
import matplotlib.pyplot as plt
import mdtraj as md

print("All imports successful.")
```

### Step 2: Load Trajectories

```python
EXP_ROOT = PROJECT_ROOT / "v3_experiments"

SYSTEMS = {
    "SPINK7-KLK5": {
        "exp": "EXP-01_spink7_klk5_binding",
        "protease_chain": 0,
        "inhibitor_chain": 1,
        "p1_resname": "ARG",
    },
    "BPTI-trypsin": {
        "exp": "EXP-04_bpti_trypsin_dg_bind",
        "protease_chain": 0,
        "inhibitor_chain": 1,
        "p1_resname": "LYS",
    },
}

trajectories = {}
for name, info in SYSTEMS.items():
    traj_path = EXP_ROOT / info["exp"] / "outputs/production/production_trajectory.dcd"
    top_path = EXP_ROOT / info["exp"] / "outputs/solvated_complex.pdb"
    if traj_path.exists():
        traj = load_trajectory(str(traj_path), str(top_path), stride=10)
        n_skip = traj.n_frames // 2
        trajectories[name] = traj[n_skip:]
        print(f"{name}: {trajectories[name].n_frames} frames")
```

### Step 3: Subsite Contact Mapping

```python
cutoff = 0.45  # 4.5 Å in nm

subsite_maps = {}

for name, traj in trajectories.items():
    info = SYSTEMS[name]
    topo = traj.topology

    # Find P1 and neighboring inhibitor residues
    inh_residues = [r for r in topo.residues if r.chain.index == info["inhibitor_chain"]]
    p1_candidates = [r for r in inh_residues if r.name == info["p1_resname"]]

    if not p1_candidates:
        print(f"{name}: P1 not found")
        continue

    p1_res = p1_candidates[0]
    p1_idx = p1_res.index

    # Define P4-P3' positions relative to P1
    positions = {"P4": p1_idx-4, "P3": p1_idx-3, "P2": p1_idx-2, "P1": p1_idx,
                 "P1'": p1_idx+1, "P2'": p1_idx+2, "P3'": p1_idx+3}

    protease_residues = [r for r in topo.residues if r.chain.index == info["protease_chain"]]

    contact_map = {}
    for pos_name, res_idx in positions.items():
        if res_idx < 0 or res_idx >= topo.n_residues:
            continue
        inh_r = topo.residue(res_idx)
        inh_atoms = [a.index for a in inh_r.atoms if a.element.symbol != 'H']

        contacts = {}
        for prot_r in protease_residues:
            prot_atoms = [a.index for a in prot_r.atoms if a.element.symbol != 'H']

            # Compute minimum distance across all frames
            occupancy = 0
            for frame_idx in range(0, traj.n_frames, max(1, traj.n_frames // 50)):
                frame = traj[frame_idx]
                pairs = [[ia, pa] for ia in inh_atoms for pa in prot_atoms]
                if pairs:
                    dists = md.compute_distances(frame, pairs)[0]
                    if np.min(dists) < cutoff:
                        occupancy += 1

            n_sampled = len(range(0, traj.n_frames, max(1, traj.n_frames // 50)))
            occ_pct = occupancy / n_sampled * 100
            if occ_pct > 20:
                contacts[f"{prot_r.name}{prot_r.resSeq}"] = occ_pct

        contact_map[pos_name] = contacts

    subsite_maps[name] = contact_map
    print(f"\n{name} subsite contacts:")
    for pos, cts in contact_map.items():
        if cts:
            top_cts = sorted(cts.items(), key=lambda x: -x[1])[:5]
            print(f"  {pos}: {', '.join(f'{r}({o:.0f}%)' for r,o in top_cts)}")
```

### Step 4: S1 Salt Bridge Analysis

```python
s1_results = {}

for name, traj in trajectories.items():
    info = SYSTEMS[name]
    topo = traj.topology

    inh_residues = [r for r in topo.residues if r.chain.index == info["inhibitor_chain"]]
    p1_candidates = [r for r in inh_residues if r.name == info["p1_resname"]]
    if not p1_candidates:
        continue

    p1_res = p1_candidates[0]
    # P1 basic nitrogen (Arg: CZ or NH1/NH2, Lys: NZ)
    if info["p1_resname"] == "ARG":
        p1_n = [a.index for a in p1_res.atoms if a.name in ("NH1", "NH2")]
    else:
        p1_n = [a.index for a in p1_res.atoms if a.name == "NZ"]

    # Find Asp in S1 pocket (closest Asp to P1)
    prot_asp = [r for r in topo.residues
                if r.chain.index == info["protease_chain"] and r.name == "ASP"]
    asp_od = []
    for asp_r in prot_asp:
        for a in asp_r.atoms:
            if a.name in ("OD1", "OD2"):
                asp_od.append((asp_r, a.index))

    if p1_n and asp_od:
        # Find closest Asp across trajectory
        best_asp = None
        best_occ = 0
        for asp_r, od_idx in asp_od:
            occ = 0
            for p1_idx in p1_n:
                dists = md.compute_distances(traj, [[p1_idx, od_idx]])[:, 0] * 10
                occ_val = np.mean(dists < 4.0) * 100
                if occ_val > best_occ:
                    best_occ = occ_val
                    best_asp = f"{asp_r.name}{asp_r.resSeq}"

        s1_results[name] = {"asp": best_asp, "occupancy": best_occ}
        print(f"{name}: S1 salt bridge to {best_asp}, occupancy = {best_occ:.1f}%")
```

### Step 5: Classification

```python
spink7_occ = s1_results.get("SPINK7-KLK5", {}).get("occupancy", 0)
spink7_subsites = subsite_maps.get("SPINK7-KLK5", {})
n_mapped = sum(1 for v in spink7_subsites.values() if v)

if spink7_occ > 80 and n_mapped >= 5:
    classification = "PASS"
elif spink7_occ > 50 and n_mapped >= 3:
    classification = "MARGINAL"
else:
    classification = "FAIL"

results = {
    "experiment_id": "EXP-21", "feature_id": "F-21",
    "s1_salt_bridge": s1_results,
    "subsite_maps": subsite_maps,
    "n_subsites_mapped": n_mapped,
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-21: {classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: Subsite Contact Heatmap

```python
fig, axes = plt.subplots(1, len(subsite_maps), figsize=(8*len(subsite_maps), 10))
if len(subsite_maps) == 1: axes = [axes]

for ax, (name, contact_map) in zip(axes, subsite_maps.items()):
    positions = list(contact_map.keys())
    all_prot_residues = set()
    for cts in contact_map.values():
        all_prot_residues.update(cts.keys())
    prot_res_list = sorted(all_prot_residues)

    matrix = np.zeros((len(positions), len(prot_res_list)))
    for i, pos in enumerate(positions):
        for j, pres in enumerate(prot_res_list):
            matrix[i, j] = contact_map[pos].get(pres, 0)

    im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=100)
    ax.set_yticks(range(len(positions)))
    ax.set_yticklabels(positions, fontsize=10)
    ax.set_xticks(range(len(prot_res_list)))
    ax.set_xticklabels(prot_res_list, fontsize=8, rotation=90)
    ax.set_title(f'{name}: Subsite Contacts', fontsize=12)
    plt.colorbar(im, ax=ax, label='Occupancy (%)')

plt.suptitle(f'EXP-21: Subsite Specificity — {classification}', fontsize=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-21_subsite_heatmap.png", dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Figure 2: S1 Salt Bridge Distance

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
# Plot S1 salt bridge distance time series for primary system
for name, info in s1_results.items():
    ax.text(0.5, 0.5, f'{name}\nS1: {info["asp"]}\nOccupancy: {info["occupancy"]:.1f}%',
            transform=ax.transAxes, fontsize=16, va='center', ha='center',
            bbox=dict(boxstyle='round', facecolor='lightgreen' if info['occupancy'] > 80 else 'lightyellow'))
ax.axis('off')
ax.set_title(f'EXP-21: S1 Salt Bridge Summary — {classification}', fontsize=14)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-21_s1_summary.png", dpi=300)
plt.close(fig)
```

### Figure 3: Subsite Mapping Diagram

```python
fig, ax = plt.subplots(1, 1, figsize=(14, 6))
positions = ["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'"]
subsites = ["S4", "S3", "S2", "S1", "S1'", "S2'", "S3'"]

for i, (pos, sub) in enumerate(zip(positions, subsites)):
    # Inhibitor residue (top)
    ax.add_patch(plt.Rectangle((i*2, 3), 1.5, 0.8, facecolor='lightblue', edgecolor='black'))
    ax.text(i*2 + 0.75, 3.4, pos, ha='center', va='center', fontsize=12, fontweight='bold')
    # Protease subsite (bottom)
    color = 'lightgreen' if sub == 'S1' else 'lightyellow'
    ax.add_patch(plt.Rectangle((i*2, 1), 1.5, 0.8, facecolor=color, edgecolor='black'))
    ax.text(i*2 + 0.75, 1.4, sub, ha='center', va='center', fontsize=12, fontweight='bold')
    # Arrow
    ax.annotate('', xy=(i*2+0.75, 1.8), xytext=(i*2+0.75, 3.0),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))

ax.set_xlim(-0.5, 14.5)
ax.set_ylim(0, 5)
ax.axis('off')
ax.set_title(f'EXP-21: Subsite–Substrate Mapping — {classification}', fontsize=15)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-21_subsite_diagram.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-21: KLK5 Subsite Specificity — Results Report

**Experiment ID:** EXP-21  **Feature ID:** F-21  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## S1 Salt Bridge
| System | Asp Residue | Occupancy (%) |
|--------|-------------|---------------|
| SPINK7-KLK5 | [res] | [val] |
| BPTI-trypsin | [res] | [val] |

## Subsite Map (SPINK7-KLK5)
| Position | Top Protease Contacts | Occupancy |
|----------|----------------------|-----------|
| P1 | [residues] | [val]% |
| ... | ... | ... |

## Figures
1. Subsite contact heatmap
2. S1 summary
3. Subsite mapping diagram

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
