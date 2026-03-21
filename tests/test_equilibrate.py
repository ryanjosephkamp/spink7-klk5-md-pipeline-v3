"""Tests for NVT and NPT equilibration utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import openmm
import pytest
from openmm import LangevinMiddleIntegrator, unit
from openmm.app import ForceField, HBonds, PME, Simulation

from src.config import EquilibrationConfig, MinimizationConfig, SystemConfig
from src.physics.restraints import create_positional_restraints
from src.prep.solvate import solvate_system
from src.prep.topology import build_topology
from src.simulate.equilibrate import run_npt, run_nvt
from src.simulate.minimizer import minimize_energy


def _write_neutral_protonated_peptide_pdb(pdb_path: Path) -> Path:
    """Write a compact neutral peptide suitable for fast equilibration tests."""

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


def _make_solvated_simulation(tmp_path: Path) -> Simulation:
    """Create a minimized solvated peptide simulation with positional restraints."""

    system_config = SystemConfig(box_padding_nm=1.2, ionic_strength_molar=0.0)
    protonated_path = _write_neutral_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "equilibration_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, system_config)
    solvated_modeller, _, _, _ = solvate_system(modeller, system_config)

    force_field = ForceField(system_config.force_field, system_config.water_model)
    system = force_field.createSystem(
        solvated_modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )

    restrained_indices: list[int] = []
    reference_positions: list[list[float]] = []
    for atom_index, atom in enumerate(solvated_modeller.topology.atoms()):
        if atom.residue.name.upper() in {"HOH", "WAT", "NA", "CL"}:
            continue
        if atom.element is not None and atom.element.symbol != "H":
            restrained_indices.append(atom_index)
            position = solvated_modeller.positions[atom_index].value_in_unit(unit.nanometer)
            reference_positions.append([float(position[0]), float(position[1]), float(position[2])])

    create_positional_restraints(
        system,
        restrained_indices,
        reference_positions=np.asarray(reference_positions, dtype=float),
        force_constant_kj_mol_nm2=500.0,
    )

    integrator = LangevinMiddleIntegrator(310.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(solvated_modeller.topology, system, integrator, platform)
    simulation.context.setPositions(solvated_modeller.positions)
    minimize_energy(simulation, MinimizationConfig(max_iterations=1000, tolerance_kj_mol_nm=10.0))
    return simulation


def test_run_nvt_temperature_within_tolerance(tmp_path: Path) -> None:
    """IV-2: NVT equilibration should maintain the average temperature within 5 K of target."""

    simulation = _make_solvated_simulation(tmp_path)
    config = EquilibrationConfig(nvt_duration_ps=100.0, friction_per_ps=10.0, save_interval_ps=0.5)

    result = run_nvt(simulation, config, tmp_path / "nvt_output")

    assert result["trajectory_path"].exists()
    assert result["final_state_path"].exists()
    assert abs(result["avg_temperature_k"] - config.temperature_k) < 5.0
    assert result["temperature_std_k"] > 0.0
    assert simulation.currentStep > 0
    assert "t0_temperature" in result
    assert isinstance(result["t0_temperature"], int)
    assert result["t0_temperature"] >= 0


def test_run_npt_density_and_temperature_physical(tmp_path: Path) -> None:
    """IV-2 and IV-3: NPT equilibration should maintain physical temperature and density."""

    simulation = _make_solvated_simulation(tmp_path)
    config = EquilibrationConfig(npt_duration_ps=100.0, friction_per_ps=10.0, barostat_interval=5, save_interval_ps=1.0)

    result = run_npt(simulation, config, tmp_path / "npt_output")

    assert result["trajectory_path"].exists()
    assert result["final_state_path"].exists()
    assert abs(result["avg_temperature_k"] - config.temperature_k) < 5.0
    assert 0.95 < result["avg_density_g_cm3"] < 1.05
    assert result["box_vectors_nm"].shape == (3, 3)
    assert simulation.currentStep > 0
    assert "t0_temperature" in result
    assert "t0_density" in result
    assert isinstance(result["t0_temperature"], int)
    assert isinstance(result["t0_density"], int)


def test_run_nvt_rejects_non_positive_duration(tmp_path: Path) -> None:
    """Public API should reject non-physical NVT duration inputs."""

    simulation = _make_solvated_simulation(tmp_path)

    with pytest.raises(ValueError, match="config.nvt_duration_ps must be positive"):
        run_nvt(simulation, EquilibrationConfig(nvt_duration_ps=0.0), tmp_path / "invalid_nvt")


# ---------- L-18 Step 3: Configurable seeds in equilibration ----------


def test_run_nvt_with_explicit_seed_is_reproducible(tmp_path: Path) -> None:
    """NVT with an explicit seed should return that seed and produce valid results."""

    simulation = _make_solvated_simulation(tmp_path)
    config = EquilibrationConfig(
        nvt_duration_ps=10.0, friction_per_ps=10.0, save_interval_ps=0.5, random_seed=42
    )

    result = run_nvt(simulation, config, tmp_path / "nvt_1")

    assert result["random_seed"] == 42
    assert abs(result["avg_temperature_k"] - config.temperature_k) < 5.0


def test_run_nvt_with_none_seed_returns_auto_seed(tmp_path: Path) -> None:
    """NVT with random_seed=None should return the auto-generated seed."""

    simulation = _make_solvated_simulation(tmp_path)
    config = EquilibrationConfig(
        nvt_duration_ps=10.0, friction_per_ps=10.0, save_interval_ps=0.5, random_seed=None
    )

    result = run_nvt(simulation, config, tmp_path / "nvt_auto")

    assert isinstance(result["random_seed"], int)
    assert result["random_seed"] >= 0


# ---------- L-19 Step 1: Topology serialization during NPT ----------


def test_run_npt_writes_topology_pdb(tmp_path: Path) -> None:
    """NPT equilibration should produce a topology_reference.pdb file."""

    simulation = _make_solvated_simulation(tmp_path)
    config = EquilibrationConfig(
        npt_duration_ps=10.0, friction_per_ps=10.0, barostat_interval=5, save_interval_ps=0.5
    )

    result = run_npt(simulation, config, tmp_path / "npt_output")

    assert result["topology_path"].exists()
    assert result["topology_path"].name == "topology_reference.pdb"
    from openmm.app import PDBFile
    pdb = PDBFile(str(result["topology_path"]))
    n_chains = sum(1 for _ in pdb.topology.chains())
    assert n_chains >= 1