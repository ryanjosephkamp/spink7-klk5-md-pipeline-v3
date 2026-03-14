"""Shared pytest fixtures for the SPINK7-KLK5 MD pipeline test suite."""

from __future__ import annotations

from pathlib import Path

import openmm
import pytest
from openmm import LangevinMiddleIntegrator, unit
from openmm.app import ForceField, HBonds, PME, Simulation

from src.config import SystemConfig
from src.prep.solvate import solvate_system
from src.prep.topology import build_topology


def _write_alanine_dipeptide_pdb(pdb_path: Path) -> Path:
    """Write a compact alanine-dipeptide test structure to disk."""

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """ATOM      1  CH3 ACE A   1      -2.300   1.200   0.000  1.00 20.00           C
ATOM      2  C   ACE A   1      -1.000   0.800   0.000  1.00 20.00           C
ATOM      3  O   ACE A   1      -0.300   1.700   0.000  1.00 20.00           O
ATOM      4  N   ALA A   2      -0.650  -0.500   0.000  1.00 20.00           N
ATOM      5  CA  ALA A   2       0.700  -1.050   0.000  1.00 20.00           C
ATOM      6  C   ALA A   2       1.850  -0.050   0.000  1.00 20.00           C
ATOM      7  O   ALA A   2       2.980  -0.430   0.000  1.00 20.00           O
ATOM      8  CB  ALA A   2       0.820  -2.570   0.000  1.00 20.00           C
ATOM      9  N   NME A   3       1.560   1.210   0.000  1.00 20.00           N
ATOM     10  CH3 NME A   3       2.540   2.220   0.000  1.00 20.00           C
TER
END
""",
        encoding="utf-8",
    )
    return pdb_path


@pytest.fixture(scope="session")
def alanine_dipeptide_pdb(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Return a compact alanine-dipeptide PDB for fast CPU integration testing."""

    root = tmp_path_factory.mktemp("alanine_dipeptide")
    return _write_alanine_dipeptide_pdb(root / "alanine_dipeptide.pdb")


@pytest.fixture
def alanine_dipeptide_simulation(alanine_dipeptide_pdb: Path) -> Simulation:
    """Build a fully solvated, parameterized alanine-dipeptide simulation on CPU."""

    system_config = SystemConfig(box_padding_nm=1.2, ionic_strength_molar=0.0)
    _, _, modeller = build_topology(alanine_dipeptide_pdb, system_config)
    solvated_modeller, _, _, _ = solvate_system(modeller, system_config)

    force_field = ForceField(system_config.force_field, system_config.water_model)
    system = force_field.createSystem(
        solvated_modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )

    integrator = LangevinMiddleIntegrator(310.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(solvated_modeller.topology, system, integrator, platform)
    simulation.context.setPositions(solvated_modeller.positions)
    return simulation


@pytest.fixture
def tmp_output_dir(tmp_path: Path) -> Path:
    """Return a temporary output directory for streamed test artifacts."""

    output_dir = tmp_path / "test_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir