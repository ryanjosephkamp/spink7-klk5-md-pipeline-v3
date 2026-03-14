"""Tests for PDB cleaning utilities."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.prep.pdb_clean import clean_structure


def _write_test_pdb(pdb_path: Path) -> Path:
    """Write a compact PDB fixture containing protein, ligand, and water records."""

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """ATOM      1  N   ALA A   1      11.104  13.207   9.331  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.000  13.000   9.000  1.00 20.00           C
HETATM    3  C1  LIG A 101      14.000  12.000   8.000  1.00 20.00           C
HETATM    4  O   HOH A 201      15.000  12.500   8.500  1.00 20.00           O
TER
ATOM      5  N   GLY B   1      21.104  23.207  19.331  1.00 20.00           N
ATOM      6  CA  GLY B   1      22.000  23.000  19.000  1.00 20.00           C
HETATM    7  O   HOH B 202      25.000  22.500  18.500  1.00 20.00           O
TER
ATOM      8  N   SER C   1      31.104  33.207  29.331  1.00 20.00           N
END
""",
        encoding="utf-8",
    )
    return pdb_path


def _atom_and_hetatm_lines(pdb_path: Path) -> list[str]:
    """Collect structural record lines from a PDB file."""

    return [
        line
        for line in pdb_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]


def test_clean_structure_keeps_only_requested_chains_by_default(tmp_path: Path) -> None:
    """Chunk 5 gate: default cleaning should keep only requested protein chains."""

    input_pdb = _write_test_pdb(tmp_path / "data" / "pdb" / "raw" / "complex.pdb")

    cleaned_path = clean_structure(input_pdb, chains_to_keep=["A", "B"])
    cleaned_lines = _atom_and_hetatm_lines(cleaned_path)

    assert cleaned_path == tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb"
    assert all(line[21].strip() in {"A", "B"} for line in cleaned_lines)
    assert all(line.startswith("ATOM") for line in cleaned_lines)
    assert len(cleaned_lines) == 4


def test_clean_structure_keeps_ligands_when_heteroatom_removal_disabled(tmp_path: Path) -> None:
    """Disabling heteroatom removal should preserve non-water HETATM records."""

    input_pdb = _write_test_pdb(tmp_path / "data" / "pdb" / "raw" / "complex.pdb")

    cleaned_path = clean_structure(
        input_pdb,
        chains_to_keep=["A"],
        remove_heteroatoms=False,
        remove_waters=True,
    )
    cleaned_lines = _atom_and_hetatm_lines(cleaned_path)

    assert [line[17:20].strip() for line in cleaned_lines] == ["ALA", "ALA", "LIG"]
    assert all(line[17:20].strip() != "HOH" for line in cleaned_lines)


def test_clean_structure_keeps_waters_when_water_removal_disabled(tmp_path: Path) -> None:
    """Disabling water removal should preserve water records in retained chains."""

    input_pdb = _write_test_pdb(tmp_path / "data" / "pdb" / "raw" / "complex.pdb")

    cleaned_path = clean_structure(
        input_pdb,
        chains_to_keep=["A"],
        remove_heteroatoms=False,
        remove_waters=False,
    )
    cleaned_lines = _atom_and_hetatm_lines(cleaned_path)

    assert [line[17:20].strip() for line in cleaned_lines] == ["ALA", "ALA", "LIG", "HOH"]


def test_clean_structure_rejects_empty_chain_selection(tmp_path: Path) -> None:
    """Public API should reject empty chain selection input."""

    input_pdb = _write_test_pdb(tmp_path / "data" / "pdb" / "raw" / "complex.pdb")

    with pytest.raises(ValueError, match="chains_to_keep must contain at least one non-empty chain identifier"):
        clean_structure(input_pdb, chains_to_keep=[])
