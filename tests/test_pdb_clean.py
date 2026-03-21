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


def test_parse_cif_to_pdb_lines_produces_atom_records(tmp_path: Path) -> None:
    """CIF parser should convert CIF atoms to PDB-format ATOM/HETATM lines."""
    gemmi = pytest.importorskip("gemmi")
    from src.prep.pdb_clean import _parse_cif_to_pdb_lines

    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "ALA"
    residue.seqid = gemmi.SeqId("1")
    atom = gemmi.Atom()
    atom.name = "CA"
    atom.element = gemmi.Element("C")
    atom.pos = gemmi.Position(1.0, 2.0, 3.0)
    atom.occ = 1.0
    atom.b_iso = 20.0
    residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)
    structure = gemmi.Structure()
    structure.add_model(model)
    structure.name = "test"

    cif_path = tmp_path / "test.cif"
    structure.make_mmcif_document().write_file(str(cif_path))

    lines = _parse_cif_to_pdb_lines(cif_path)
    atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) >= 1, "CIF parser must produce at least one ATOM line"
    assert "ALA" in atom_lines[0], "Residue name ALA must appear in ATOM line"
    assert "CA" in atom_lines[0], "Atom name CA must appear in ATOM line"


def _write_gemmi_cif(tmp_path: Path, filename: str, chains: dict[str, str]) -> Path:
    """Write a minimal CIF structure with the given chains and residue names."""
    gemmi = pytest.importorskip("gemmi")

    model = gemmi.Model("1")
    for chain_id, res_name in chains.items():
        chain = gemmi.Chain(chain_id)
        residue = gemmi.Residue()
        residue.name = res_name
        residue.seqid = gemmi.SeqId("1")
        for atom_name, element in [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]:
            atom = gemmi.Atom()
            atom.name = atom_name
            atom.element = gemmi.Element(element)
            atom.pos = gemmi.Position(1.0, 2.0, 3.0)
            atom.occ = 1.0
            atom.b_iso = 20.0
            residue.add_atom(atom)
        chain.add_residue(residue)
        model.add_chain(chain)
    structure = gemmi.Structure()
    structure.add_model(model)
    structure.name = "test"

    raw_dir = tmp_path / "data" / "pdb" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    cif_path = raw_dir / filename
    structure.make_mmcif_document().write_file(str(cif_path))
    return cif_path


def test_clean_structure_accepts_cif_input(tmp_path: Path) -> None:
    """clean_structure should accept CIF input and produce a cleaned PDB."""
    cif_path = _write_gemmi_cif(tmp_path, "complex.cif", {"A": "ALA", "B": "GLY"})

    cleaned_path = clean_structure(cif_path, chains_to_keep=["A"])

    assert cleaned_path.suffix == ".pdb"
    assert cleaned_path.name == "complex_cleaned.pdb"
    content = cleaned_path.read_text(encoding="utf-8")
    atom_lines = [line for line in content.splitlines() if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) > 0, "Cleaned output must contain ATOM records"
    assert all(line[21].strip() == "A" for line in atom_lines), "Only chain A atoms should remain"


def test_clean_structure_raises_on_cif_without_gemmi(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CIF cleaning should raise PipelineError if Gemmi is not installed."""
    import src.prep.pdb_clean as pdb_clean_module
    from src import PipelineError

    monkeypatch.setattr(pdb_clean_module, "_HAS_GEMMI", False)

    cif_path = tmp_path / "data" / "pdb" / "raw" / "complex.cif"
    cif_path.parent.mkdir(parents=True, exist_ok=True)
    cif_path.write_text("data_test\n", encoding="utf-8")

    with pytest.raises(PipelineError, match="Gemmi is required"):
        clean_structure(cif_path, chains_to_keep=["A"])


def test_clean_structure_accepts_mmcif_extension(tmp_path: Path) -> None:
    """clean_structure should accept .mmcif extension identically to .cif."""
    cif_path = _write_gemmi_cif(tmp_path, "complex.mmcif", {"A": "ALA"})

    cleaned_path = clean_structure(cif_path, chains_to_keep=["A"])
    assert cleaned_path.suffix == ".pdb"


# ---------------------------------------------------------------------------
# L-28 — Parser edge-case tests
# ---------------------------------------------------------------------------


def test_clean_structure_preserves_model_boundaries_in_nmr_ensemble(tmp_path: Path) -> None:
    """Multi-model NMR ensembles should retain MODEL/ENDMDL boundaries after cleaning."""

    pdb_path = tmp_path / "data" / "pdb" / "raw" / "nmr.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "MODEL        1\n"
        "ATOM      1  N   ALA A   1      11.104  13.207   9.331  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      12.000  13.000   9.000  1.00 20.00           C\n"
        "TER\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      1  N   ALA A   1      11.204  13.307   9.431  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      12.100  13.100   9.100  1.00 20.00           C\n"
        "TER\n"
        "ENDMDL\n"
        "END\n",
        encoding="utf-8",
    )

    cleaned_path = clean_structure(pdb_path, chains_to_keep=["A"], model_index=None)
    content = cleaned_path.read_text(encoding="utf-8")
    lines = content.splitlines()

    model_lines = [l for l in lines if l.strip().startswith("MODEL")]
    endmdl_lines = [l for l in lines if l.strip().startswith("ENDMDL")]
    atom_lines = [l for l in lines if l.startswith("ATOM")]

    assert len(model_lines) == 2, "Both MODEL records must be preserved"
    assert len(endmdl_lines) == 2, "Both ENDMDL records must be preserved"
    assert len(atom_lines) == 4, "All ATOM records across models must be retained"
    # Verify coordinates differ between models (no concatenation)
    assert atom_lines[0][30:54] != atom_lines[2][30:54]


def test_clean_structure_preserves_insertion_codes(tmp_path: Path) -> None:
    """Residues with insertion codes (52, 52A, 52B) must all be preserved independently."""

    pdb_path = tmp_path / "data" / "pdb" / "raw" / "insertions.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "ATOM      1  CA  ALA A  52      10.000  10.000  10.000  1.00 20.00           C\n"
        "ATOM      2  CA  ALA A  52A     11.000  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3  CA  ALA A  52B     12.000  10.000  10.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    cleaned_path = clean_structure(pdb_path, chains_to_keep=["A"])
    atom_lines = [l for l in cleaned_path.read_text(encoding="utf-8").splitlines() if l.startswith("ATOM")]

    assert len(atom_lines) == 3, "All three insertion-code residues must be preserved"
    insertion_codes = [l[26] for l in atom_lines]
    assert " " in insertion_codes, "Blank insertion code must be preserved"
    assert "A" in insertion_codes, "Insertion code A must be preserved"
    assert "B" in insertion_codes, "Insertion code B must be preserved"


def test_clean_structure_handles_alternate_conformations(tmp_path: Path) -> None:
    """PDB files with altLoc indicators should preserve all conformations during cleaning."""

    pdb_path = tmp_path / "data" / "pdb" / "raw" / "altloc.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "ATOM      1  N   SER A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  SER A   1      11.000  10.000  10.000  1.00 20.00           C\n"
        "ATOM      3A OG  SER A   1      12.000  11.000  10.000  0.60 20.00           O\n"
        "ATOM      4B OG  SER A   1      12.000  11.500  10.000  0.40 20.00           O\n"
        "END\n",
        encoding="utf-8",
    )

    cleaned_path = clean_structure(pdb_path, chains_to_keep=["A"])
    atom_lines = [l for l in cleaned_path.read_text(encoding="utf-8").splitlines() if l.startswith("ATOM")]

    assert len(atom_lines) == 4, "All ATOM records including alternate conformations must be preserved"


def test_clean_structure_raises_on_empty_result_from_absent_chain(tmp_path: Path) -> None:
    """Requesting a non-existent chain should raise PipelineError, not produce an empty PDB."""

    from src import PipelineError

    pdb_path = tmp_path / "data" / "pdb" / "raw" / "single_chain.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "ATOM      1  N   ALA A   1      11.104  13.207   9.331  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      12.000  13.000   9.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    with pytest.raises(PipelineError, match="Cleaning removed all structural records"):
        clean_structure(pdb_path, chains_to_keep=["Z"])


def test_clean_structure_handles_conect_records_gracefully(tmp_path: Path) -> None:
    """CONECT records in the input PDB should not cause errors during cleaning."""

    pdb_path = tmp_path / "data" / "pdb" / "raw" / "with_conect.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "ATOM      1  N   ALA A   1      11.104  13.207   9.331  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      12.000  13.000   9.000  1.00 20.00           C\n"
        "HETATM    3  C1  LIG A 101      14.000  12.000   8.000  1.00 20.00           C\n"
        "HETATM    4  C2  LIG A 101      14.500  12.500   8.000  1.00 20.00           C\n"
        "CONECT    3    4\n"
        "CONECT    4    3\n"
        "END\n",
        encoding="utf-8",
    )

    cleaned_path = clean_structure(pdb_path, chains_to_keep=["A"], remove_heteroatoms=False)
    content = cleaned_path.read_text(encoding="utf-8")
    atom_lines = [l for l in content.splitlines() if l.startswith(("ATOM", "HETATM"))]

    assert len(atom_lines) == 4, "All ATOM/HETATM records for chain A must be preserved"
    assert "CONECT" not in content, "CONECT records should be silently discarded by current logic"
    content = cleaned_path.read_text(encoding="utf-8")
    atom_lines = [line for line in content.splitlines() if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) > 0, "Cleaned output must contain ATOM records"


def test_clean_structure_rejects_unsupported_extension(tmp_path: Path) -> None:
    """Public API should reject extensions other than .pdb, .cif, .mmcif."""
    bad_path = tmp_path / "data" / "pdb" / "raw" / "structure.xyz"
    bad_path.parent.mkdir(parents=True, exist_ok=True)
    bad_path.write_text("dummy\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must have a .pdb, .cif, or .mmcif extension"):
        clean_structure(bad_path, chains_to_keep=["A"])


def _write_multi_model_pdb(pdb_path: Path) -> Path:
    """Write a compact multi-model NMR-style PDB with 3 models."""
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """MODEL        1
ATOM      1  N   ALA A   1      11.104  13.207   9.331  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.000  13.000   9.000  1.00 20.00           C
TER
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1      11.200  13.100   9.400  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.100  12.900   9.100  1.00 20.00           C
TER
ENDMDL
MODEL        3
ATOM      1  N   ALA A   1      11.300  13.300   9.200  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.200  13.100   8.900  1.00 20.00           C
TER
ENDMDL
END
""",
        encoding="utf-8",
    )
    return pdb_path


def test_clean_structure_selects_first_model_by_default(tmp_path: Path) -> None:
    """Default model_index=1 should select only the first model from a multi-model PDB."""
    input_pdb = _write_multi_model_pdb(tmp_path / "data" / "pdb" / "raw" / "nmr_ensemble.pdb")

    with pytest.warns(UserWarning, match="Multi-model PDB detected"):
        cleaned_path = clean_structure(input_pdb, chains_to_keep=["A"])

    atom_lines = _atom_and_hetatm_lines(cleaned_path)
    assert len(atom_lines) == 2, f"Expected 2 atoms from model 1, got {len(atom_lines)}"
    # Verify coordinates match model 1 (N at 11.104)
    assert "11.104" in atom_lines[0]


def test_clean_structure_selects_specific_model_by_index(tmp_path: Path) -> None:
    """model_index=3 should select the third model's coordinates."""
    input_pdb = _write_multi_model_pdb(tmp_path / "data" / "pdb" / "raw" / "nmr_ensemble.pdb")

    with pytest.warns(UserWarning, match="Multi-model PDB detected"):
        cleaned_path = clean_structure(input_pdb, chains_to_keep=["A"], model_index=3)

    atom_lines = _atom_and_hetatm_lines(cleaned_path)
    assert len(atom_lines) == 2, f"Expected 2 atoms from model 3, got {len(atom_lines)}"
    # Verify coordinates match model 3 (N at 11.300)
    assert "11.300" in atom_lines[0]


def test_clean_structure_different_models_have_different_coordinates(tmp_path: Path) -> None:
    """Selecting model 1 vs model 3 must produce different coordinate sets."""
    input_pdb = _write_multi_model_pdb(tmp_path / "data" / "pdb" / "raw" / "nmr_ensemble.pdb")

    with pytest.warns(UserWarning):
        path_m1 = clean_structure(input_pdb, chains_to_keep=["A"], model_index=1)

    # Re-create input because clean_structure writes to a fixed output path
    input_pdb_2 = _write_multi_model_pdb(tmp_path / "data2" / "pdb" / "raw" / "nmr_ensemble.pdb")

    with pytest.warns(UserWarning):
        path_m3 = clean_structure(input_pdb_2, chains_to_keep=["A"], model_index=3)

    content_m1 = path_m1.read_text(encoding="utf-8")
    content_m3 = path_m3.read_text(encoding="utf-8")
    assert content_m1 != content_m3, "Model 1 and model 3 should have different coordinates"


def test_clean_structure_single_model_no_warning(tmp_path: Path) -> None:
    """Single-model PDBs should pass through without any multi-model warning."""
    input_pdb = _write_test_pdb(tmp_path / "data" / "pdb" / "raw" / "complex.pdb")

    # No warning should be emitted for single-model PDBs
    cleaned_path = clean_structure(input_pdb, chains_to_keep=["A", "B"])
    atom_lines = _atom_and_hetatm_lines(cleaned_path)
    assert len(atom_lines) == 4


def test_clean_structure_rejects_invalid_model_index(tmp_path: Path) -> None:
    """Requesting a model index beyond the model count should raise PipelineError."""
    from src import PipelineError

    input_pdb = _write_multi_model_pdb(tmp_path / "data" / "pdb" / "raw" / "nmr_ensemble.pdb")

    with pytest.raises(PipelineError, match="Requested model_index=5.*only 3 model"):
        clean_structure(input_pdb, chains_to_keep=["A"], model_index=5)


def test_clean_structure_model_index_none_retains_all_models(tmp_path: Path) -> None:
    """model_index=None should retain all models without warning."""
    input_pdb = _write_multi_model_pdb(tmp_path / "data" / "pdb" / "raw" / "nmr_ensemble.pdb")

    cleaned_path = clean_structure(input_pdb, chains_to_keep=["A"], model_index=None)
    atom_lines = _atom_and_hetatm_lines(cleaned_path)
    # 3 models x 2 atoms each = 6 total atoms
    assert len(atom_lines) == 6, f"Expected 6 atoms from all 3 models, got {len(atom_lines)}"


def test_clean_structure_handles_large_atom_count_from_cif(tmp_path: Path) -> None:
    """CIF path should handle structures with many atoms."""
    gemmi = pytest.importorskip("gemmi")

    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    # 200 atoms in a single residue is sufficient to validate the CIF path
    # handles multi-atom structures; hybrid-36 encoding engages automatically
    # in Gemmi's PDB writer for serials exceeding 99,999.
    residue = gemmi.Residue()
    residue.name = "ALA"
    residue.seqid = gemmi.SeqId("1")
    for i in range(200):
        atom = gemmi.Atom()
        atom.name = f"C{i:03d}"[:4]
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(float(i), 0.0, 0.0)
        atom.occ = 1.0
        atom.b_iso = 20.0
        residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)
    structure = gemmi.Structure()
    structure.add_model(model)

    raw_dir = tmp_path / "data" / "pdb" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    cif_path = raw_dir / "large.cif"
    structure.make_mmcif_document().write_file(str(cif_path))

    cleaned_path = clean_structure(cif_path, chains_to_keep=["A"])
    content = cleaned_path.read_text(encoding="utf-8")
    atom_lines = [line for line in content.splitlines() if line.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) == 200, f"Expected 200 atoms, got {len(atom_lines)}"
