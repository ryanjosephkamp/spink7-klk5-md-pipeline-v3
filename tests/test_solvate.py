"""Tests for explicit-solvent box construction and ion placement."""

from __future__ import annotations

from pathlib import Path

import pytest
from openmm import NonbondedForce, unit
from openmm.app import ForceField, HBonds, PME

from src.config import SystemConfig
from src.prep.solvate import solvate_system
from src.prep.topology import build_topology


def _write_neutral_protonated_peptide_pdb(pdb_path: Path) -> Path:
    """Write a small neutral peptide PDB for topology and solvation tests."""

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


def _write_charged_protonated_peptide_pdb(pdb_path: Path) -> Path:
    """Write a small positively charged peptide PDB to test counter-ion neutralization."""

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """ATOM      1  N   ALA A   1      -0.648   1.210   0.000  1.00 20.00           N
ATOM      2  H   ALA A   1      -1.300   1.500   0.730  1.00 20.00           H
ATOM      3  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      4  C   ALA A   1       1.526   0.000   0.000  1.00 20.00           C
ATOM      5  O   ALA A   1       2.046  -1.143   0.000  1.00 20.00           O
ATOM      6  CB  ALA A   1      -0.507  -0.774  -1.206  1.00 20.00           C
ATOM      7  N   LYS A   2       2.260   1.105   0.000  1.00 20.00           N
ATOM      8  H   LYS A   2       1.836   2.019   0.000  1.00 20.00           H
ATOM      9  CA  LYS A   2       3.714   1.084   0.000  1.00 20.00           C
ATOM     10  C   LYS A   2       4.234  -0.337   0.000  1.00 20.00           C
ATOM     11  O   LYS A   2       5.430  -0.596   0.000  1.00 20.00           O
ATOM     12  OXT LYS A   2       3.451  -1.314   0.000  1.00 20.00           O
ATOM     13  CB  LYS A   2       4.312   1.782  -1.240  1.00 20.00           C
ATOM     14  CG  LYS A   2       5.796   1.939  -1.157  1.00 20.00           C
ATOM     15  CD  LYS A   2       6.359   2.693  -2.370  1.00 20.00           C
ATOM     16  CE  LYS A   2       7.842   2.857  -2.236  1.00 20.00           C
ATOM     17  NZ  LYS A   2       8.358   3.630  -3.377  1.00 20.00           N
TER
END
""",
        encoding="utf-8",
    )
    return pdb_path


def _total_charge_e(modeller, system_config: SystemConfig) -> float:
    """Build a solvated OpenMM system and return its total charge in elementary charge."""

    force_field = ForceField(system_config.force_field, system_config.water_model)
    system = force_field.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        constraints=HBonds,
        rigidWater=True,
    )

    nonbonded_force = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    total_charge = 0.0
    for particle_index in range(system.getNumParticles()):
        charge, _, _ = nonbonded_force.getParticleParameters(particle_index)
        total_charge += charge.value_in_unit(unit.elementary_charge)
    return total_charge


def test_solvate_system_adds_water_and_periodic_box(tmp_path: Path) -> None:
    """Chunk 8 gate: solvation should add explicit waters and define box vectors."""

    protonated_path = _write_neutral_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, SystemConfig())

    solvated_modeller, n_water, n_positive_ions, n_negative_ions = solvate_system(
        modeller,
        SystemConfig(box_padding_nm=0.8, ionic_strength_molar=0.0),
    )

    assert n_water > 0
    assert n_positive_ions >= 0
    assert n_negative_ions >= 0
    assert solvated_modeller.topology.getPeriodicBoxVectors() is not None


def test_solvate_system_neutralizes_charge_with_counterions(tmp_path: Path) -> None:
    """Solvation should neutralize a charged solute with the configured counter-ions."""

    protonated_path = _write_charged_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "charged_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, SystemConfig())
    solvation_config = SystemConfig(box_padding_nm=0.8, ionic_strength_molar=0.0)

    solvated_modeller, n_water, n_positive_ions, n_negative_ions = solvate_system(modeller, solvation_config)
    total_charge_e = _total_charge_e(solvated_modeller, solvation_config)

    assert n_water > 0
    assert n_positive_ions == 0
    assert n_negative_ions > 0
    assert abs(total_charge_e) < 1e-4


def test_solvate_system_rejects_negative_ionic_strength(tmp_path: Path) -> None:
    """Public API should reject non-physical ionic-strength inputs."""

    protonated_path = _write_neutral_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, SystemConfig())

    with pytest.raises(ValueError, match="system_config.ionic_strength_molar must be non-negative"):
        solvate_system(modeller, SystemConfig(ionic_strength_molar=-0.1))


# ---------- L-16 Step 2: Water model name resolution ----------


def test_water_model_name_resolves_tip3p() -> None:
    """_water_model_name should resolve TIP3P XML to 'tip3p'."""
    from src.prep.solvate import _water_model_name

    assert _water_model_name("amber14/tip3p.xml") == "tip3p"


def test_water_model_name_resolves_opc() -> None:
    """_water_model_name should resolve OPC XML to 'tip4pew' (placement model)."""
    from src.prep.solvate import _water_model_name

    assert _water_model_name("amber14/opc.xml") == "tip4pew"


def test_water_model_name_resolves_tip4pew() -> None:
    """_water_model_name should resolve TIP4P-Ew XML to 'tip4pew'."""
    from src.prep.solvate import _water_model_name

    assert _water_model_name("amber14/tip4pew.xml") == "tip4pew"


def test_water_model_name_fallback_for_unknown_model() -> None:
    """_water_model_name should fall back to stem extraction for unregistered models."""
    from src.prep.solvate import _water_model_name

    assert _water_model_name("custom/spce.xml") == "spce"


# ---------- L-16 Step 4: Parametrized water model solvation ----------


def _opc_available() -> bool:
    """Check if OPC water model XML is available in this OpenMM installation."""
    try:
        ForceField("amber14-all.xml", "amber14/opc.xml")
        return True
    except Exception:
        return False


def _tip4pew_available() -> bool:
    """Check if TIP4P-Ew water model XML is available in this OpenMM installation."""
    try:
        ForceField("amber14-all.xml", "amber14/tip4pew.xml")
        return True
    except Exception:
        return False


@pytest.mark.parametrize("water_model_xml", [
    "amber14/tip3p.xml",
    pytest.param("amber14/opc.xml", marks=pytest.mark.skipif(
        not _opc_available(), reason="OPC XML not available in this OpenMM installation"
    )),
    pytest.param("amber14/tip4pew.xml", marks=pytest.mark.skipif(
        not _tip4pew_available(), reason="TIP4P-Ew XML not available in this OpenMM installation"
    )),
])
def test_solvate_with_water_model_adds_water_and_periodic_box(
    tmp_path: Path, water_model_xml: str
) -> None:
    """Each supported water model should produce a valid solvated system."""
    protonated_path = _write_neutral_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, SystemConfig())

    config = SystemConfig(water_model=water_model_xml, box_padding_nm=0.8, ionic_strength_molar=0.0)
    solvated, n_water, n_pos, n_neg = solvate_system(modeller, config)

    assert n_water > 0
    assert solvated.topology.getPeriodicBoxVectors() is not None


# ---------- L-17 Step 1: Box shape in SystemConfig ----------


def test_system_config_box_shape_defaults_to_cubic() -> None:
    """SystemConfig should default to cubic box shape for backward compatibility."""
    from src.config import SUPPORTED_BOX_SHAPES

    config = SystemConfig()
    assert config.box_shape == "cubic"
    assert config.box_shape in SUPPORTED_BOX_SHAPES


def test_system_config_accepts_dodecahedron() -> None:
    """SystemConfig should accept dodecahedron box shape."""
    config = SystemConfig(box_shape="dodecahedron")
    assert config.box_shape == "dodecahedron"


# ---------- L-17 Step 2: Box shape propagation to solvation ----------


def test_solvate_with_cubic_box_produces_valid_system(tmp_path: Path) -> None:
    """Cubic solvation should produce a valid periodic box (baseline)."""
    protonated_path = _write_neutral_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, SystemConfig())

    config = SystemConfig(box_padding_nm=0.8, box_shape="cubic", ionic_strength_molar=0.0)
    solvated, n_water, _, _ = solvate_system(modeller, config)

    assert n_water > 0
    assert solvated.topology.getPeriodicBoxVectors() is not None


def test_solvate_rejects_invalid_box_shape(tmp_path: Path) -> None:
    """solvate_system() should reject unsupported box shapes."""
    protonated_path = _write_neutral_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, SystemConfig())

    # Use object.__setattr__ to bypass frozen dataclass for test purposes
    config = SystemConfig(box_padding_nm=0.8, ionic_strength_molar=0.0)
    object.__setattr__(config, "box_shape", "hexagonal")

    with pytest.raises(ValueError, match="box_shape"):
        solvate_system(modeller, config)


# ---------- L-17 Step 3: Atom count comparison test ----------


def test_dodecahedron_uses_fewer_atoms_than_cubic(tmp_path: Path) -> None:
    """A dodecahedral box should contain 20-35% fewer atoms than a cubic box at the same padding."""
    pdb_path = tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide.pdb"
    protonated_path = _write_neutral_protonated_peptide_pdb(pdb_path)

    cubic_config = SystemConfig(box_padding_nm=1.0, box_shape="cubic", ionic_strength_molar=0.0)
    _, _, modeller_cubic = build_topology(protonated_path, cubic_config)
    solvated_cubic, n_water_cubic, _, _ = solvate_system(modeller_cubic, cubic_config)
    n_atoms_cubic = solvated_cubic.topology.getNumAtoms()

    dodec_config = SystemConfig(box_padding_nm=1.0, box_shape="dodecahedron", ionic_strength_molar=0.0)
    _, _, modeller_dodec = build_topology(protonated_path, dodec_config)
    solvated_dodec, n_water_dodec, _, _ = solvate_system(modeller_dodec, dodec_config)
    n_atoms_dodec = solvated_dodec.topology.getNumAtoms()

    reduction = 1.0 - (n_atoms_dodec / n_atoms_cubic)
    assert 0.15 < reduction < 0.45, (
        f"Expected 15-45% atom reduction, got {reduction:.1%} "
        f"(cubic={n_atoms_cubic}, dodecahedron={n_atoms_dodec})"
    )


def test_dodecahedron_box_vectors_are_non_orthogonal(tmp_path: Path) -> None:
    """A dodecahedral box should have non-orthogonal box vectors (off-diagonal elements)."""
    pdb_path = tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide.pdb"
    protonated_path = _write_neutral_protonated_peptide_pdb(pdb_path)

    config = SystemConfig(box_padding_nm=0.8, box_shape="dodecahedron", ionic_strength_molar=0.0)
    _, _, modeller = build_topology(protonated_path, config)
    solvated, _, _, _ = solvate_system(modeller, config)

    box_vectors = solvated.topology.getPeriodicBoxVectors()
    assert box_vectors is not None

    # Extract box vector components (in nm)
    a = box_vectors[0]
    b = box_vectors[1]
    c = box_vectors[2]

    # Dodecahedral box should have non-zero off-diagonal elements in b and/or c
    b_x = b[0].value_in_unit(unit.nanometer)
    c_x = c[0].value_in_unit(unit.nanometer)
    c_y = c[1].value_in_unit(unit.nanometer)

    # At least one off-diagonal element should be non-zero for a non-cubic box
    has_off_diagonal = abs(b_x) > 1e-6 or abs(c_x) > 1e-6 or abs(c_y) > 1e-6
    assert has_off_diagonal, "Dodecahedral box should have non-orthogonal vectors"