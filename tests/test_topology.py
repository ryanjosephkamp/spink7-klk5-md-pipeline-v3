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
        tmp_path / "data" / "pdb" / "prepared" / "peptide.pdb"
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


def test_build_topology_default_uses_nocutoff(tmp_path: Path) -> None:
    """Default nonbonded method should be NoCutoff (pre-solvation, no box)."""

    protonated_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide.pdb"
    )
    _, system, _ = build_topology(protonated_path, SystemConfig())

    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if hasattr(force, "getNonbondedMethod"):
            assert force.getNonbondedMethod() == 0, (
                "Default nonbonded method should be NoCutoff (enum value 0)"
            )
            break


def test_build_topology_small_system_nocutoff_no_warning(tmp_path: Path) -> None:
    """NoCutoff on a small system should NOT emit a UserWarning."""

    import warnings

    protonated_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide.pdb"
    )
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        build_topology(protonated_path, SystemConfig())
        nocutoff_warnings = [x for x in w if "NoCutoff" in str(x.message)]
        assert len(nocutoff_warnings) == 0, "Small system should not trigger NoCutoff warning"


def test_build_topology_accepts_cutoff_nonperiodic(tmp_path: Path) -> None:
    """build_topology should accept CutoffNonPeriodic with nonbonded_cutoff_nm."""

    from openmm.app import CutoffNonPeriodic

    protonated_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide.pdb"
    )
    topology, system, _ = build_topology(
        protonated_path,
        SystemConfig(),
        nonbonded_method=CutoffNonPeriodic,
        nonbonded_cutoff_nm=0.9,
    )
    assert system.getNumParticles() == topology.getNumAtoms()

    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if hasattr(force, "getNonbondedMethod"):
            # CutoffNonPeriodic = enum value 1
            assert force.getNonbondedMethod() == 1
            break


def test_build_topology_system_nonbonded_method_is_inspectable(tmp_path: Path) -> None:
    """The returned system's NonbondedForce should report its nonbonded method."""

    protonated_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide.pdb"
    )
    _, system, _ = build_topology(protonated_path, SystemConfig())

    nonbonded_found = False
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if hasattr(force, "getNonbondedMethod"):
            nonbonded_found = True
            # NoCutoff = 0 in OpenMM's NonbondedForce enum
            assert force.getNonbondedMethod() == 0
            break
    assert nonbonded_found, "System should contain a NonbondedForce"


def test_build_topology_backward_compatible_signature(tmp_path: Path) -> None:
    """Calling build_topology with only (pdb_path, system_config) must still work."""

    protonated_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide.pdb"
    )
    # Two-argument call — must not raise
    topology, system, modeller = build_topology(protonated_path, SystemConfig())
    assert topology.getNumAtoms() > 0
    assert system.getNumParticles() == topology.getNumAtoms()


def test_build_topology_skip_hydrogens_preserves_atom_count(tmp_path: Path) -> None:
    """skip_hydrogens=True should not add or remove any atoms from the input PDB."""

    from openmm.app import PDBFile

    # Build a fully protonated structure via addHydrogens first
    bare_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide.pdb"
    )
    _, _, modeller_full = build_topology(bare_path, SystemConfig())

    # Save the fully protonated structure to a new PDB
    full_h_path = tmp_path / "data" / "pdb" / "prepared" / "full_h_protonated.pdb"
    with full_h_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller_full.topology, modeller_full.positions, handle)

    input_atom_count = modeller_full.topology.getNumAtoms()

    # Re-read with skip_hydrogens=True — atom count must be preserved exactly
    topology_skip, _, _ = build_topology(full_h_path, SystemConfig(), skip_hydrogens=True)
    assert topology_skip.getNumAtoms() == input_atom_count


def test_build_topology_auto_detects_pre_protonated_input(tmp_path: Path) -> None:
    """build_topology should auto-skip addHydrogens for *_protonated.pdb files with H atoms."""

    from openmm.app import PDBFile

    # Build a fully protonated structure first
    bare_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide_bare.pdb"
    )
    _, _, modeller_full = build_topology(bare_path, SystemConfig())

    # Save with _protonated naming convention (has H atoms)
    protonated_path = tmp_path / "data" / "pdb" / "prepared" / "peptide_protonated.pdb"
    with protonated_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller_full.topology, modeller_full.positions, handle)

    input_atom_count = modeller_full.topology.getNumAtoms()

    # Default skip_hydrogens=False, but auto-detection should kick in
    topology, _, _ = build_topology(protonated_path, SystemConfig())
    assert topology.getNumAtoms() == input_atom_count


def test_build_topology_does_not_auto_skip_for_non_protonated_name(tmp_path: Path) -> None:
    """build_topology should still add hydrogens when filename lacks _protonated suffix."""

    pdb_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide_cleaned.pdb"
    )

    input_atom_count = sum(
        1 for line in pdb_path.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    )

    topology, _, _ = build_topology(pdb_path, SystemConfig())
    # OpenMM should add hydrogens, so atom count should be >= input
    assert topology.getNumAtoms() >= input_atom_count


def _write_histidine_peptide_pdb(pdb_path: Path) -> Path:
    """Write a minimal HIS-containing dipeptide PDB for tautomer testing."""

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """ATOM      1  N   ALA A   1      -0.648   1.210   0.000  1.00 20.00           N
ATOM      2  H   ALA A   1      -1.300   1.500   0.730  1.00 20.00           H
ATOM      3  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      4  C   ALA A   1       1.526   0.000   0.000  1.00 20.00           C
ATOM      5  O   ALA A   1       2.046  -1.143   0.000  1.00 20.00           O
ATOM      6  CB  ALA A   1      -0.507  -0.774  -1.206  1.00 20.00           C
ATOM      7  N   HIS A   2       2.260   1.105   0.000  1.00 20.00           N
ATOM      8  CA  HIS A   2       3.714   1.084   0.000  1.00 20.00           C
ATOM      9  C   HIS A   2       4.234  -0.337   0.000  1.00 20.00           C
ATOM     10  O   HIS A   2       5.430  -0.596   0.000  1.00 20.00           O
ATOM     11  CB  HIS A   2       4.200   1.800  -1.260  1.00 20.00           C
ATOM     12  CG  HIS A   2       5.680   1.750  -1.400  1.00 20.00           C
ATOM     13  ND1 HIS A   2       6.400   0.600  -1.200  1.00 20.00           N
ATOM     14  CE1 HIS A   2       7.700   0.870  -1.350  1.00 20.00           C
ATOM     15  NE2 HIS A   2       7.860   2.130  -1.650  1.00 20.00           N
ATOM     16  CD2 HIS A   2       6.560   2.650  -1.710  1.00 20.00           C
ATOM     17  OXT HIS A   2       3.451  -1.314   0.000  1.00 20.00           O
TER
END
""",
        encoding="utf-8",
    )
    return pdb_path


def test_build_topology_preserves_histidine_tautomer_with_skip(tmp_path: Path) -> None:
    """skip_hydrogens=True must preserve the histidine tautomer from upstream protonation."""

    from openmm.app import PDBFile

    # Build a fully protonated structure containing HIS; OpenMM assigns a tautomer
    his_path = _write_histidine_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "histidine.pdb"
    )
    topo_orig, _, modeller_orig = build_topology(his_path, SystemConfig())

    # Determine which tautomer OpenMM chose
    orig_res_names = [res.name for res in topo_orig.residues()]
    his_variants = {"HIS", "HIE", "HID", "HIP"}
    his_name = next(n for n in orig_res_names if n in his_variants)

    # Save fully protonated structure with _protonated naming convention
    full_h_path = tmp_path / "data" / "pdb" / "prepared" / "histidine_protonated.pdb"
    with full_h_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller_orig.topology, modeller_orig.positions, handle)

    # Re-read with skip_hydrogens=True — tautomer must be preserved
    topo_skip, _, _ = build_topology(full_h_path, SystemConfig(), skip_hydrogens=True)
    skip_res_names = [res.name for res in topo_skip.residues()]
    assert his_name in skip_res_names, (
        f"Histidine tautomer {his_name} was not preserved; found: {skip_res_names}"
    )

    # Verify atom counts match exactly (no extra/missing hydrogens)
    assert topo_skip.getNumAtoms() == topo_orig.getNumAtoms()


def test_build_topology_warns_on_double_protonation(tmp_path: Path) -> None:
    """Calling addHydrogens on a pre-hydrogenated structure should emit a warning."""

    # Use a non-_protonated filename so auto-detection does NOT kick in,
    # but the PDB itself contains H atoms from the test helper
    pdb_path = _write_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "peptide_cleaned.pdb"
    )

    with pytest.warns(UserWarning, match="already contained.*hydrogen atoms"):
        build_topology(pdb_path, SystemConfig(), skip_hydrogens=False)