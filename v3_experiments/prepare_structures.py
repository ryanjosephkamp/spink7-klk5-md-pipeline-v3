#!/usr/bin/env python3
"""Prepare 1BRS (barnase-barstar) using PDBFixer for pipeline validation.

PDBFixer resolves alternate conformations, missing atoms, and missing residues
that cause topology building failures with raw PDB files.
This is a V3 analysis script (§21.3 permitted action), not a pipeline modification.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path("/Users/noir/visual_studio/Visual_Studio__UC_Spring_26/"
                    "CS_RES_SELF_STUDY/medium_projects/medium_project_2")
sys.path.insert(0, str(PROJECT_ROOT))

from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Fix 1BRS
raw_pdb = PROJECT_ROOT / "data" / "pdb" / "raw" / "1BRS.pdb"
output_dir = PROJECT_ROOT / "v3_experiments" / "structures"
output_dir.mkdir(parents=True, exist_ok=True)

print("Fixing 1BRS with PDBFixer...")
fixer = PDBFixer(str(raw_pdb))

# Keep chains A (barnase) and D (barstar) only
chains_to_remove = []
for chain in fixer.topology.chains():
    if chain.id not in ("A", "D"):
        chains_to_remove.append(chain.id)
fixer.removeChains(chainIds=chains_to_remove)

# Fix standard operations
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.4)

output_path = output_dir / "1BRS_AD_fixed.pdb"
with open(output_path, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

# Count atoms and chains
n_atoms = fixer.topology.getNumAtoms()
n_residues = sum(1 for _ in fixer.topology.residues())
n_chains = sum(1 for _ in fixer.topology.chains())
print(f"Fixed 1BRS: {n_atoms} atoms, {n_residues} residues, {n_chains} chains")
print(f"Saved: {output_path}")

# Also fix 2PSX (KLK5) for experiments that need it
raw_2psx = PROJECT_ROOT / "data" / "pdb" / "raw" / "2PSX.pdb"
if raw_2psx.exists():
    print("\nFixing 2PSX (KLK5)...")
    fixer2 = PDBFixer(str(raw_2psx))
    chains_to_remove = [c.id for c in fixer2.topology.chains() if c.id != "A"]
    fixer2.removeChains(chainIds=chains_to_remove)
    fixer2.findMissingResidues()
    fixer2.findNonstandardResidues()
    fixer2.replaceNonstandardResidues()
    fixer2.removeHeterogens(keepWater=False)
    fixer2.findMissingAtoms()
    fixer2.addMissingAtoms()
    fixer2.addMissingHydrogens(pH=7.4)
    
    output_2psx = output_dir / "2PSX_A_fixed.pdb"
    with open(output_2psx, "w") as f:
        PDBFile.writeFile(fixer2.topology, fixer2.positions, f)
    
    n2 = fixer2.topology.getNumAtoms()
    nr2 = sum(1 for _ in fixer2.topology.residues())
    print(f"Fixed 2PSX: {n2} atoms, {nr2} residues")
    print(f"Saved: {output_2psx}")

print("\nStructure preparation complete.")
