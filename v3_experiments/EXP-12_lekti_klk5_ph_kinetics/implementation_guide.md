# EXP-12: LEKTI–KLK5 pH-Dependent Binding Kinetics — Implementation Guide

**Experiment ID:** EXP-12  
**Feature ID:** F-12 (benchmarks.md)  
**Category:** Thermodynamic/Kinetic (Quantitative)  
**Date:** 2026-03-22  
**Phase:** Step 4 Phase B — Implementation Guide  

---

## Part 1 — Complete Experimental Design

### 1. Abstract

Tests the pH-dependent binding switch of LEKTI domain 6 to KLK5. LEKTI inhibits KLK5 at neutral pH but releases at acidic pH (stratum corneum). The experiment runs parallel US/WHAM simulations at pH 4.5, 5.5, 6.5, and 7.5, predicting ΔΔG(pH 4.5 − pH 7.5) > +2 kcal/mol (weaker binding at acidic pH). The pH switch is driven by histidine protonation at the interface.

### 2. Benchmark Values

| Metric | Value | Source |
|--------|-------|--------|
| ΔΔG(pH 4.5 − pH 7.5) | > +2 kcal/mol | Deraison 2007, Borgono 2007 |
| pH_opt(KLK5) | ~7.5 | Brattsand 2005 |
| Expected trend | Monotonic weakening with decreasing pH | Literature consensus |

### 3. Protocol Summary

For each of 4 pH values (4.5, 5.5, 6.5, 7.5):
1. Assign protonation states via PROPKA at target pH
2. Build and solvate system
3. Equilibrate (minimize→NVT→NPT)
4. Run US/WHAM along COM distance CV
5. Extract ΔG_bind at each pH
6. Compute ΔΔG(pH 4.5 − pH 7.5)

### 4. Classification Criteria (§25.1)

- PASS: ΔΔG > +2 kcal/mol AND monotonic pH trend
- MARGINAL: ΔΔG between +1 and +2 kcal/mol
- FAIL: ΔΔG ≤ 0 or > +10 kcal/mol or non-monotonic
- σ_method = 2.0 kcal/mol per pH condition

---

## Part 2 — Step-by-Step Implementation Instructions

### Step 1: Environment Setup

```python
import os, sys, json
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

EXP_DIR = PROJECT_ROOT / "v3_experiments" / "EXP-12_lekti_klk5_ph_kinetics"
OUTPUT_DIR = EXP_DIR / "outputs"
FIGURES_DIR = EXP_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

from src.config import (SystemConfig, MinimizationConfig, EquilibrationConfig,
                         ProductionConfig, UmbrellaConfig, WHAMConfig)
from src.prep.pdb_fetch import fetch_pdb, fetch_alphafold
from src.prep.pdb_clean import clean_structure
from src.prep.protonate import assign_protonation
from src.prep.topology import build_topology
from src.prep.solvate import solvate_system
from src.simulate.minimizer import minimize_energy
from src.simulate.equilibrate import run_nvt, run_npt
from src.simulate.production import run_production
from src.simulate.umbrella import run_umbrella_campaign, generate_window_centers
from src.simulate.platform import select_platform
from src.analyze.wham import solve_wham, bootstrap_pmf_uncertainty
from src.analyze.mbar_analysis import solve_mbar
from src.analyze.trajectory import load_trajectory
from src.analyze.structural import compute_rmsd
from src.physics.collective_variables import com_distance
import matplotlib.pyplot as plt

print("All imports successful.")
```

### Step 2: Fetch Structures (LEKTI AlphaFold + KLK5)

```python
data_dir = OUTPUT_DIR / "structures"
data_dir.mkdir(exist_ok=True)

# KLK5 from PDB 2PSX
klk5_pdb = fetch_pdb("2PSX", data_dir)
klk5_clean = clean_structure(klk5_pdb, chains_to_keep=["A"],
                              remove_heteroatoms=True, remove_waters=True)

# LEKTI from AlphaFold (UniProt Q9NQ38)
lekti_af = fetch_alphafold("Q9NQ38", data_dir)

# Extract domain 6 (residues ~603-668 in full LEKTI)
# This must be done with PDB manipulation or mdtraj
import mdtraj as md
lekti_full = md.load(str(lekti_af))
domain6_residues = list(range(602, 668))  # 0-indexed
domain6_atoms = lekti_full.topology.select(f"resid {domain6_residues[0]} to {domain6_residues[-1]}")
lekti_domain6 = lekti_full.atom_slice(domain6_atoms)
lekti_d6_path = data_dir / "lekti_domain6.pdb"
lekti_domain6.save(str(lekti_d6_path))
print(f"LEKTI domain 6 extracted: {len(domain6_atoms)} atoms")

# Use docked complex from EXP-03 or ClusPro
complex_pdb = data_dir / "klk5_lekti_d6_complex.pdb"
print(f"Place docked complex at: {complex_pdb}")
```

### Step 3: pH-Dependent Protonation and System Setup

```python
pH_VALUES = [4.5, 5.5, 6.5, 7.5]
systems = {}

for ph in pH_VALUES:
    ph_dir = OUTPUT_DIR / f"pH_{ph:.1f}"
    ph_dir.mkdir(exist_ok=True)

    # Assign protonation at target pH
    protonated = assign_protonation(complex_pdb, ph=ph, force_field="AMBER", use_propka=True)

    # Build topology and solvate
    sys_config = SystemConfig()
    from openmm.app import PME, PDBFile
    topology, system, modeller = build_topology(protonated, sys_config,
                                                 nonbonded_method=PME, nonbonded_cutoff_nm=1.0)
    modeller, n_waters, n_pos, n_neg = solvate_system(modeller, sys_config)

    solvated_pdb = ph_dir / "solvated.pdb"
    with open(solvated_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)

    systems[ph] = {
        "topology": topology, "system": system, "modeller": modeller,
        "solvated_pdb": solvated_pdb, "output_dir": ph_dir,
        "n_waters": n_waters, "n_pos": n_pos, "n_neg": n_neg,
    }
    print(f"pH {ph}: {n_waters} waters, {n_pos} Na+, {n_neg} Cl-")
```

### Step 4: Equilibration at Each pH

```python
import openmm
from openmm import LangevinMiddleIntegrator
from openmm.app import Simulation

for ph, sdata in systems.items():
    ph_dir = sdata["output_dir"]
    eq_dir = ph_dir / "equilibration"
    eq_dir.mkdir(exist_ok=True)

    platform = select_platform()
    integrator = LangevinMiddleIntegrator(310*openmm.unit.kelvin,
                                          1.0/openmm.unit.picosecond,
                                          0.002*openmm.unit.picoseconds)
    sim = Simulation(sdata["modeller"].topology, sdata["system"], integrator, platform)
    sim.context.setPositions(sdata["modeller"].positions)

    minimize_energy(sim, MinimizationConfig())
    run_nvt(sim, EquilibrationConfig(), eq_dir)
    run_npt(sim, EquilibrationConfig(), eq_dir)

    sdata["simulation"] = sim
    print(f"pH {ph}: equilibration complete")
```

### Step 5: US/WHAM at Each pH

```python
umbrella_config = UmbrellaConfig()
wham_config = WHAMConfig()

dg_by_ph = {}
pmf_by_ph = {}

for ph, sdata in systems.items():
    us_dir = sdata["output_dir"] / "umbrella"
    us_dir.mkdir(exist_ok=True)

    window_centers = generate_window_centers(umbrella_config)
    us_results = run_umbrella_campaign(sdata["simulation"], umbrella_config, us_dir,
                                        cv_type="com_distance", group1_sel="chainid 0",
                                        group2_sel="chainid 1")

    pmf_xi, pmf_G, pmf_err = solve_wham(us_results, wham_config)
    boot_err = bootstrap_pmf_uncertainty(us_results, wham_config, n_bootstrap=200)

    dg_bind = np.min(pmf_G) - pmf_G[-1]
    dg_err = boot_err[np.argmin(pmf_G)]

    dg_by_ph[ph] = {"dg": float(dg_bind), "err": float(dg_err)}
    pmf_by_ph[ph] = {"xi": pmf_xi.tolist(), "G": pmf_G.tolist(), "err": pmf_err.tolist()}

    print(f"pH {ph}: ΔG_bind = {dg_bind:.2f} ± {dg_err:.2f} kcal/mol")
```

### Step 6: Classification

```python
dg_acidic = dg_by_ph[4.5]["dg"]
dg_neutral = dg_by_ph[7.5]["dg"]
ddg = dg_acidic - dg_neutral  # positive = weaker at acidic pH

# Check monotonicity
ph_sorted = sorted(dg_by_ph.keys())
dg_sorted = [dg_by_ph[p]["dg"] for p in ph_sorted]
monotonic = all(dg_sorted[i] >= dg_sorted[i+1] for i in range(len(dg_sorted)-1))

sigma_method = 2.0
sigma_combined = np.sqrt(2 * sigma_method**2)  # two pH conditions

if ddg > 2.0 and monotonic:
    classification = "PASS"
elif 1.0 <= ddg <= 2.0:
    classification = "MARGINAL"
elif ddg <= 0 or ddg > 10 or not monotonic:
    classification = "FAIL"
else:
    classification = "INCONCLUSIVE"

results = {
    "experiment_id": "EXP-12", "feature_id": "F-12",
    "dg_by_ph": dg_by_ph,
    "ddg_acidic_minus_neutral": float(ddg),
    "monotonic": monotonic,
    "sigma_method": sigma_method,
    "classification": classification,
}
with open(EXP_DIR / "results.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"EXP-12: ΔΔG = {ddg:.2f} kcal/mol, monotonic={monotonic}, classification={classification}")
```

---

## Part 3 — Figure Generation Instructions

### Figure 1: PMF Profiles Across pH Values

```python
fig, ax = plt.subplots(1, 1, figsize=(12, 7))
colors = {4.5: 'red', 5.5: 'orange', 6.5: 'steelblue', 7.5: 'darkblue'}

for ph in pH_VALUES:
    xi = np.array(pmf_by_ph[ph]["xi"])
    G = np.array(pmf_by_ph[ph]["G"])
    err = np.array(pmf_by_ph[ph]["err"])
    ax.plot(xi * 10, G, color=colors[ph], linewidth=2, label=f'pH {ph}')
    ax.fill_between(xi * 10, G - err, G + err, alpha=0.15, color=colors[ph])

ax.set_xlabel('COM Distance (Å)', fontsize=14)
ax.set_ylabel('PMF (kcal/mol)', fontsize=14)
ax.set_title('EXP-12: pH-Dependent LEKTI–KLK5 PMF', fontsize=15)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-12_pmf_ph_overlay.png", dpi=300)
plt.close(fig)
```

### Figure 2: ΔG vs pH

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
phs = sorted(dg_by_ph.keys())
dgs = [dg_by_ph[p]["dg"] for p in phs]
errs = [dg_by_ph[p]["err"] for p in phs]

ax.errorbar(phs, dgs, yerr=errs, fmt='o-', color='navy', markersize=10,
            linewidth=2, capsize=5, capthick=2)
ax.set_xlabel('pH', fontsize=14)
ax.set_ylabel('ΔG_bind (kcal/mol)', fontsize=14)
ax.set_title(f'EXP-12: pH Switch — ΔΔG = {ddg:.2f} kcal/mol [{classification}]', fontsize=14)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-12_dg_vs_ph.png", dpi=300)
plt.close(fig)
```

### Figure 3: Histidine Protonation State Map

```python
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
# Show protonation state assigned by PROPKA at each pH for interface histidines
his_residues = ["His37", "His103", "His142"]  # Example interface His
proton_states = {
    4.5: [1, 1, 1],  # all protonated (HIP)
    5.5: [1, 1, 0],
    6.5: [0, 1, 0],
    7.5: [0, 0, 0],  # all neutral (HID/HIE)
}
x = np.arange(len(his_residues))
width = 0.2
for i, ph in enumerate(pH_VALUES):
    ax.bar(x + i*width, proton_states[ph], width, label=f'pH {ph}', color=colors[ph])

ax.set_xticks(x + 1.5*width)
ax.set_xticklabels(his_residues, fontsize=12)
ax.set_ylabel('Protonated (1) / Neutral (0)', fontsize=12)
ax.set_title('EXP-12: Interface His Protonation States', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-12_his_protonation.png", dpi=300)
plt.close(fig)
```

### Figure 4: ΔΔG Bar Chart

```python
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ref_ph = 7.5
ddg_values = [dg_by_ph[p]["dg"] - dg_by_ph[ref_ph]["dg"] for p in phs]
bar_colors = ['red' if d > 2 else 'orange' if d > 1 else 'steelblue' for d in ddg_values]
ax.bar([str(p) for p in phs], ddg_values, color=bar_colors, edgecolor='black')
ax.axhline(y=2.0, color='green', linestyle='--', label='PASS threshold (+2)')
ax.axhline(y=1.0, color='orange', linestyle=':', label='MARGINAL threshold (+1)')
ax.set_xlabel('pH', fontsize=14)
ax.set_ylabel('ΔΔG vs pH 7.5 (kcal/mol)', fontsize=14)
ax.set_title(f'EXP-12: pH-Dependent ΔΔG — {classification}', fontsize=14)
ax.legend()
plt.tight_layout()
fig.savefig(FIGURES_DIR / "EXP-12_ddg_bar.png", dpi=300)
plt.close(fig)
```

---

## Part 4 — Results Documentation Template

```markdown
# EXP-12: LEKTI–KLK5 pH-Dependent Binding — Results Report

**Experiment ID:** EXP-12  **Feature ID:** F-12  **Date:** [date]  **Classification:** [PASS/MARGINAL/FAIL]

## Abstract
Tested pH-dependent binding switch of LEKTI D6–KLK5 via parallel US/WHAM at pH 4.5–7.5.

## Results
| pH | ΔG_bind (kcal/mol) | SE |
|----|--------------------|----|
| 4.5 | [value] | [SE] |
| 5.5 | [value] | [SE] |
| 6.5 | [value] | [SE] |
| 7.5 | [value] | [SE] |

ΔΔG(pH 4.5 − pH 7.5) = [value] ± [SE] kcal/mol
Benchmark: > +2 kcal/mol
Monotonic trend: [Yes/No]

## Figures
1. PMF overlay across pH values
2. ΔG vs pH curve
3. Interface His protonation map
4. ΔΔG bar chart

---
Author: Ryan Kamp / Dept. of Computer Science, University of Cincinnati / kamprj@mail.uc.edu / GitHub: ryanjosephkamp
```

---

Author: Ryan Kamp  
Affiliation: Dept. of Computer Science, University of Cincinnati  
Contact:  
Email: kamprj@mail.uc.edu  
GitHub: ryanjosephkamp
