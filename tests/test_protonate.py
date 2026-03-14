"""Tests for protonation-state assignment utilities."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.prep.protonate import assign_protonation


def _write_cleaned_test_pdb(pdb_path: Path) -> Path:
    """Write a compact cleaned PDB containing titratable residues."""

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """ATOM      1  N   HIS A   1      11.104  13.207   9.331  1.00 20.00           N
ATOM      2  CA  HIS A   1      12.000  13.000   9.000  1.00 20.00           C
ATOM      3  CG  HIS A   1      12.600  14.200   8.500  1.00 20.00           C
ATOM      4  ND1 HIS A   1      13.900  14.300   8.300  1.00 20.00           N
ATOM      5  CE1 HIS A   1      14.100  15.400   7.700  1.00 20.00           C
ATOM      6  NE2 HIS A   1      13.000  16.000   7.500  1.00 20.00           N
ATOM      7  CD2 HIS A   1      11.700  15.200   8.000  1.00 20.00           C
TER
ATOM      8  N   ASP A   2      15.104  16.207  10.331  1.00 20.00           N
ATOM      9  CA  ASP A   2      16.000  16.000  10.000  1.00 20.00           C
ATOM     10  CG  ASP A   2      16.900  17.200   9.700  1.00 20.00           C
ATOM     11  OD1 ASP A   2      18.100  17.100   9.400  1.00 20.00           O
ATOM     12  OD2 ASP A   2      16.300  18.300   9.800  1.00 20.00           O
END
""",
        encoding="utf-8",
    )
    return pdb_path


def _atom_lines(pdb_path: Path) -> list[str]:
    """Collect ATOM lines from a PDB file."""

    return [line for line in pdb_path.read_text(encoding="utf-8").splitlines() if line.startswith("ATOM")]


def test_assign_protonation_adds_hydrogens_and_logs_decisions(tmp_path: Path) -> None:
    """Chunk 6 gate: protonation should add hydrogens and log titratable-state decisions."""

    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb")

    protonated_path = assign_protonation(input_pdb, ph=7.4, force_field="AMBER")
    atom_lines = _atom_lines(protonated_path)
    protonation_log = protonated_path.with_name("complex_cleaned_protonation.log").read_text(encoding="utf-8")

    assert protonated_path == tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned_protonated.pdb"
    assert any(line[12:16].strip() == "HE2" and line[17:20].strip() == "HIE" for line in atom_lines)
    assert any(line[17:20].strip() == "ASP" for line in atom_lines)
    assert "HIS A 1 -> HIE" in protonation_log
    assert "ASP A 2 -> ASP" in protonation_log


def test_assign_protonation_protonates_aspartate_below_pka(tmp_path: Path) -> None:
    """Low-pH protonation should switch ASP to ASH and add the sidechain proton."""

    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb")

    protonated_path = assign_protonation(input_pdb, ph=3.0, force_field="AMBER")
    atom_lines = _atom_lines(protonated_path)

    assert any(line[17:20].strip() == "ASH" for line in atom_lines)
    assert any(line[12:16].strip() == "HD2" and line[17:20].strip() == "ASH" for line in atom_lines)


def test_assign_protonation_rejects_unsupported_force_field(tmp_path: Path) -> None:
    """Public API should reject non-AMBER force-field identifiers."""

    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb")

    with pytest.raises(ValueError, match="force_field must be 'AMBER'"):
        assign_protonation(input_pdb, force_field="CHARMM")


def test_assign_protonation_rejects_non_positive_ph(tmp_path: Path) -> None:
    """Public API should reject non-physical pH inputs."""

    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb")

    with pytest.raises(ValueError, match="ph must be positive"):
        assign_protonation(input_pdb, ph=0.0)
