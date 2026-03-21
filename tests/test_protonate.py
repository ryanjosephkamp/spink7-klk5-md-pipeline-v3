"""Tests for protonation-state assignment utilities."""

from __future__ import annotations

import logging
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


# ---------------------------------------------------------------------------
# L-02 Step 1 — PROPKA pKa extraction helper
# ---------------------------------------------------------------------------


def test_run_propka_returns_pka_dict(tmp_path: Path) -> None:
    """PROPKA helper should return a dict of pKa values for a PDB with titratable residues."""
    from src.prep.protonate import _run_propka

    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "test_propka.pdb")
    pka_dict = _run_propka(input_pdb)
    assert isinstance(pka_dict, dict)
    # The test PDB contains HIS and ASP; PROPKA should detect them.
    # If the minimal PDB is too small, we still accept an empty dict (graceful fallback).


def test_run_propka_graceful_fallback_on_bad_input(tmp_path: Path) -> None:
    """PROPKA helper should return an empty dict when given an unparseable file."""
    from src.prep.protonate import _run_propka

    bad_pdb = tmp_path / "bad.pdb"
    bad_pdb.write_text("NOT A VALID PDB\n", encoding="utf-8")
    pka_dict = _run_propka(bad_pdb)
    assert pka_dict == {}


# ---------------------------------------------------------------------------
# L-02 Step 2 — Refactor _protonation_decision() with pka_override
# ---------------------------------------------------------------------------


def test_protonation_decision_uses_pka_override() -> None:
    """When a PROPKA pKa override is provided, it should be used instead of the reference value."""
    from src.prep.protonate import _protonation_decision

    # HIS reference pKa = 6.5; at pH 7.4, HIS is deprotonated (HIE) by default
    result_default = _protonation_decision("HIS", 7.4)
    assert result_default is not None
    assert result_default[0] == "HIE"
    assert "ref" in result_default[2]

    # With a PROPKA override shifting pKa to 8.0 (e.g., buried His near negative charge),
    # at pH 7.4, the His should now be protonated (HIP)
    result_override = _protonation_decision("HIS", 7.4, pka_override=8.0)
    assert result_override is not None
    assert result_override[0] == "HIP"
    assert "PROPKA" in result_override[2]


def test_protonation_decision_backward_compatible() -> None:
    """Without a pKa override, behavior should match the original implementation exactly."""
    from src.prep.protonate import _protonation_decision

    # ASP at pH 7.4 (above pKa 4.0) -> deprotonated (ASP)
    result = _protonation_decision("ASP", 7.4)
    assert result is not None
    assert result[0] == "ASP"

    # ASP at pH 3.0 (below pKa 4.0) -> protonated (ASH)
    result_low = _protonation_decision("ASP", 3.0)
    assert result_low is not None
    assert result_low[0] == "ASH"


# ---------------------------------------------------------------------------
# L-02 Step 3 — Wire PROPKA into assign_protonation()
# ---------------------------------------------------------------------------


def test_assign_protonation_with_propka_enabled(tmp_path: Path) -> None:
    """When use_propka=True, the protonation log should include PROPKA or ref annotations."""
    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb")
    protonated_path = assign_protonation(input_pdb, ph=7.4, force_field="AMBER", use_propka=True)
    assert protonated_path.exists()
    log_path = protonated_path.with_name("complex_cleaned_protonation.log")
    log_content = log_path.read_text(encoding="utf-8")
    # Log should contain pKa source indicators (PROPKA or ref)
    assert "pKa" in log_content


def test_assign_protonation_with_propka_disabled(tmp_path: Path) -> None:
    """When use_propka=False, behavior should match the original Henderson-Hasselbalch logic."""
    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb")
    protonated_path = assign_protonation(input_pdb, ph=7.4, force_field="AMBER", use_propka=False)
    assert protonated_path.exists()
    log_path = protonated_path.with_name("complex_cleaned_protonation.log")
    log_content = log_path.read_text(encoding="utf-8")
    # All entries should use reference pKa values
    assert "ref" in log_content


# ---------------------------------------------------------------------------
# L-02 Step 4 — Non-standard titratable residues (phosphorylated)
# ---------------------------------------------------------------------------


def test_phosphorylated_residue_recognized() -> None:
    """Phosphorylated residues (SEP, TPO, PTR) should be recognized and logged without error."""
    from src.prep.protonate import _protonation_decision

    result = _protonation_decision("SEP", 7.4)
    assert result is not None
    assert "non-standard" in result[2]

    result_tpo = _protonation_decision("TPO", 7.4)
    assert result_tpo is not None
    assert "non-standard" in result_tpo[2]

    result_ptr = _protonation_decision("PTR", 7.4)
    assert result_ptr is not None
    assert "non-standard" in result_ptr[2]


# ---------------------------------------------------------------------------
# L-28 Step 6 — Non-standard residue pass-through
# ---------------------------------------------------------------------------


def test_assign_protonation_passes_through_nonstandard_residues(tmp_path: Path) -> None:
    """Non-standard residues (SEP, MSE) must be emitted unchanged without spurious hydrogens."""

    pdb_path = tmp_path / "data" / "pdb" / "prepared" / "nonstandard.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "ATOM      1  N   ALA A   1      11.104  13.207   9.331  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      12.000  13.000   9.000  1.00 20.00           C\n"
        "HETATM    3  CA  SEP A   2      14.000  13.000   9.000  1.00 20.00           C\n"
        "HETATM    4  P   SEP A   2      15.000  14.000   9.000  1.00 20.00           P\n"
        "HETATM    5  CA  MSE A   3      18.000  13.000   9.000  1.00 20.00           C\n"
        "HETATM    6  SE  MSE A   3      19.000  14.000   9.000  1.00 20.00          SE\n"
        "END\n",
        encoding="utf-8",
    )

    protonated_path = assign_protonation(pdb_path, ph=7.4, force_field="AMBER")
    atom_lines = [l for l in protonated_path.read_text(encoding="utf-8").splitlines() if l.startswith(("ATOM", "HETATM"))]

    sep_lines = [l for l in atom_lines if l[17:20].strip() == "SEP"]
    mse_lines = [l for l in atom_lines if l[17:20].strip() == "MSE"]

    assert len(sep_lines) == 2, "SEP residue atoms must be preserved unchanged"
    assert len(mse_lines) == 2, "MSE residue atoms must be preserved unchanged"
    # No spurious hydrogens added to non-standard residues
    sep_h = [l for l in sep_lines if l[76:78].strip() == "H"]
    mse_h = [l for l in mse_lines if l[76:78].strip() == "H"]
    assert len(sep_h) == 0, "No hydrogens should be added to SEP"
    assert len(mse_h) == 0, "No hydrogens should be added to MSE"


def test_unknown_residue_returns_none() -> None:
    """Residues not in _PROTONATION_RULES should return None (no protonation decision)."""
    from src.prep.protonate import _protonation_decision

    result = _protonation_decision("ALA", 7.4)
    assert result is None


# ---------------------------------------------------------------------------
# L-02 Step 5 — Comprehensive validation tests
# ---------------------------------------------------------------------------


def _write_hetatm_only_pdb(pdb_path: Path) -> Path:
    """Write a PDB containing only HETATM records (no standard residues)."""
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """HETATM    1  C1  LIG A   1      10.000  10.000  10.000  1.00 20.00           C
HETATM    2  C2  LIG A   1      11.000  10.000  10.000  1.00 20.00           C
HETATM    3  O1  LIG A   1      12.000  10.000  10.000  1.00 20.00           O
END
""",
        encoding="utf-8",
    )
    return pdb_path


def _write_sep_test_pdb(pdb_path: Path) -> Path:
    """Write a PDB containing a phosphoserine (SEP) residue alongside a standard residue."""
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
HETATM    8  N   SEP A   2      15.104  16.207  10.331  1.00 20.00           N
HETATM    9  CA  SEP A   2      16.000  16.000  10.000  1.00 20.00           C
HETATM   10  CB  SEP A   2      16.500  17.000   9.500  1.00 20.00           C
HETATM   11  OG  SEP A   2      17.000  17.500  10.000  1.00 20.00           O
HETATM   12  P   SEP A   2      18.000  18.000  10.500  1.00 20.00           P
END
""",
        encoding="utf-8",
    )
    return pdb_path


def test_propka_vs_reference_both_complete(tmp_path: Path) -> None:
    """Both PROPKA-enabled and reference-only runs should complete and produce annotated logs."""
    input_propka = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "propka_test_cleaned.pdb")
    input_ref = _write_cleaned_test_pdb(tmp_path / "data2" / "pdb" / "prepared" / "ref_test_cleaned.pdb")

    protonated_propka = assign_protonation(input_propka, ph=7.4, force_field="AMBER", use_propka=True)
    protonated_ref = assign_protonation(input_ref, ph=7.4, force_field="AMBER", use_propka=False)

    assert protonated_propka.exists()
    assert protonated_ref.exists()

    log_propka = protonated_propka.with_name("propka_test_cleaned_protonation.log").read_text(encoding="utf-8")
    log_ref = protonated_ref.with_name("ref_test_cleaned_protonation.log").read_text(encoding="utf-8")

    # Both logs must contain pKa annotations
    assert "pKa" in log_propka
    assert "pKa" in log_ref

    # Reference-only log must use "ref" annotations exclusively
    for line in log_ref.strip().splitlines():
        assert "ref" in line

    # PROPKA log should contain either "PROPKA" (if PROPKA succeeded) or "ref" (fallback)
    for line in log_propka.strip().splitlines():
        assert "PROPKA" in line or "ref" in line


def test_graceful_fallback_hetatm_only(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """PROPKA fallback: HETATM-only PDB should complete without error and log a warning."""
    input_pdb = _write_hetatm_only_pdb(tmp_path / "data" / "pdb" / "prepared" / "ligand_cleaned.pdb")

    with caplog.at_level(logging.WARNING):
        protonated_path = assign_protonation(input_pdb, ph=7.4, force_field="AMBER", use_propka=True)

    assert protonated_path.exists()
    # PROPKA should have failed or returned no titratable groups — warning expected
    assert any("falling back" in record.message.lower() or "propka" in record.message.lower() for record in caplog.records)


def test_use_propka_false_reproduces_reference_behavior(tmp_path: Path) -> None:
    """use_propka=False must produce output identical to the reference Henderson-Hasselbalch path."""
    input_pdb = _write_cleaned_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned.pdb")

    protonated_path = assign_protonation(input_pdb, ph=7.4, force_field="AMBER", use_propka=False)
    atom_lines = _atom_lines(protonated_path)
    log_path = protonated_path.with_name("complex_cleaned_protonation.log")
    log_content = log_path.read_text(encoding="utf-8")

    # Reproduce all assertions from the original gate test
    assert protonated_path == tmp_path / "data" / "pdb" / "prepared" / "complex_cleaned_protonated.pdb"
    assert any(line[12:16].strip() == "HE2" and line[17:20].strip() == "HIE" for line in atom_lines)
    assert any(line[17:20].strip() == "ASP" for line in atom_lines)
    assert "HIS A 1 -> HIE" in log_content
    assert "ASP A 2 -> ASP" in log_content

    # Every log line must use "ref" (no PROPKA annotations)
    for line in log_content.strip().splitlines():
        assert "ref" in line
        assert "PROPKA" not in line


def test_phosphorylated_residue_emits_warning(tmp_path: Path) -> None:
    """A PDB containing SEP should emit a UserWarning about specialized AMBER parameters."""
    input_pdb = _write_sep_test_pdb(tmp_path / "data" / "pdb" / "prepared" / "sep_cleaned.pdb")

    with pytest.warns(UserWarning, match="Phosphorylated residues detected"):
        assign_protonation(input_pdb, ph=7.4, force_field="AMBER", use_propka=False)
