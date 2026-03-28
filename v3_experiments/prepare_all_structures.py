#!/usr/bin/env python3
"""Prepare all needed PDB structures with PDBFixer for V3 experiments."""
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from pathlib import Path

STRUCT_DIR = Path(__file__).parent / "structures"

configs = [
    ("2PTC", ["E", "I"], "2PTC_EI_fixed.pdb"),   # E=trypsin, I=BPTI
    ("4PTI", None, "4PTI_fixed.pdb"),             # all chains
    ("1TGS", ["Z", "I"], "1TGS_ZI_fixed.pdb"),   # Z=chymotrypsinogen, I=PSTI
    ("4EIK", ["A", "B"], "4EIK_fixed.pdb"),       # A=SH3, B=p41
]

for pdb_id, chains, out_name in configs:
    pdb_path = STRUCT_DIR / f"{pdb_id}.pdb"
    if not pdb_path.exists():
        print(f"SKIP {pdb_id}: {pdb_path} not found")
        continue

    fixer = PDBFixer(str(pdb_path))

    if chains:
        chain_ids = [c.id for c in fixer.topology.chains()]
        to_remove = [i for i, cid in enumerate(chain_ids) if cid not in chains]
        fixer.removeChains(to_remove)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    out_path = STRUCT_DIR / out_name
    with open(out_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    n_atoms = fixer.topology.getNumAtoms()
    n_res = sum(1 for _ in fixer.topology.residues())
    n_chains = sum(1 for _ in fixer.topology.chains())
    print(f"{pdb_id} -> {out_name}: {n_atoms} atoms, {n_res} residues, {n_chains} chains")
