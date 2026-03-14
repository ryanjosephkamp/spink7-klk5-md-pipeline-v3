"""Tests for OpenMM topology construction utilities."""

from __future__ import annotations

from pathlib import Path

import pytest
from openmm import unit

from src import PhysicalValidityError
from src.config import SystemConfig
from src.prep.topology import build_topology


def _write_protonated_peptide_pdb(pdb_path: Path) -> Path:
    """Write a compact protonated peptide PDB for AMBER topology mapping."""

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """ATOM      1  N   ALA A   1      -0.648   1.210   0.000  1.00 20.00           N
ATOM      2  H   ALA A   1      -1.300   1.500   0.730  1.00 20.00           H
ATOM      3  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      4  C   ALA A   1       1.526   0.000   0.000  1.00 20.00           C
ATOM      5  O   ALA A   1       2.046  -1.143   0.000  1.00 20.00           O
ATOM      6  CB  ALA A   1      -0.507  -0.774  -1.206  1.00 20.00           C
ATOM      7  N   GLY A   2       2.260   1.105   0.000  1.00 20.00           N
ATOM      8  H   GLY A   2       1.836   2.019   0.000  1.00 20.00           H
ATOM      9  CA  GLY A   2       3.714   1.084   0.000  1.00 20.00           C
ATOM     10  C   GLY A   2       4.234  -0.337   0.000  1.00 20.00           C
ATOM     11  O   GLY A   2       5.430  -0.596   0.000  1.00 20.00           O
ATOM     12  OXT GLY A   2       3.451  -1.314   0.000  1.00 20.00           O
TER
END
""",
        encoding="utf-8",
    )
    return pdb_path


def test_build_topology_maps_all_atoms_without_parameter_exceptions(tmp_path: Path) -> None:
    """Chunk 7 gate: AMBER mapping should produce a fully parameterized OpenMM system."""

    protonated_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide_protonated.pdb"
    )

    topology, system, modeller = build_topology(protonated_path, SystemConfig())

    hydrogen_count = sum(
        1 for atom in topology.atoms() if atom.element is not None and atom.element.symbol == "H"
    )
    masses = [system.getParticleMass(index).value_in_unit(unit.dalton) for index in range(system.getNumParticles())]

    assert topology.getNumAtoms() == modeller.topology.getNumAtoms()
    assert system.getNumParticles() == topology.getNumAtoms()
    assert hydrogen_count > 0
    assert all(mass > 0.0 for mass in masses)


def test_build_topology_rejects_zero_atom_input(tmp_path: Path) -> None:
    """Topology generation should fail fast for physically invalid empty structures."""

    empty_pdb = tmp_path / "data" / "pdb" / "prepared" / "empty_protonated.pdb"
    empty_pdb.parent.mkdir(parents=True, exist_ok=True)
    empty_pdb.write_text("END\n", encoding="utf-8")

    with pytest.raises(PhysicalValidityError, match="valid non-empty PDB|zero atoms"):
        build_topology(empty_pdb, SystemConfig())


def test_build_topology_rejects_non_pdb_input(tmp_path: Path) -> None:
    """Public API should reject unsupported structure-file extensions."""

    invalid_path = tmp_path / "data" / "pdb" / "prepared" / "peptide_protonated.cif"
    invalid_path.parent.mkdir(parents=True, exist_ok=True)
    invalid_path.write_text("data_test\n", encoding="utf-8")

    with pytest.raises(ValueError, match="pdb_path must have a .pdb extension"):
        build_topology(invalid_path, SystemConfig())