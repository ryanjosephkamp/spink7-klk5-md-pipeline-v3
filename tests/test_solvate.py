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
        tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide_protonated.pdb"
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
        tmp_path / "data" / "pdb" / "prepared" / "charged_peptide_protonated.pdb"
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
        tmp_path / "data" / "pdb" / "prepared" / "neutral_peptide_protonated.pdb"
    )
    _, _, modeller = build_topology(protonated_path, SystemConfig())

    with pytest.raises(ValueError, match="system_config.ionic_strength_molar must be non-negative"):
        solvate_system(modeller, SystemConfig(ionic_strength_molar=-0.1))