"""Tests for umbrella-sampling utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import openmm
import pytest
from openmm import XmlSerializer, unit
from openmm.app import Element, Simulation, Topology

from src.config import ProductionConfig, UmbrellaConfig
from src.physics.collective_variables import com_distance
from src.simulate.umbrella import generate_window_centers, run_umbrella_campaign, run_umbrella_window


def _make_umbrella_test_simulation() -> tuple[Simulation, list[int], list[int]]:
    """Create a small two-chain system suitable for umbrella-restraint tests."""

    topology = Topology()
    chain_a = topology.addChain()
    chain_b = topology.addChain()
    residue_a = topology.addResidue("AAA", chain_a)
    residue_b = topology.addResidue("BBB", chain_b)
    carbon = Element.getByAtomicNumber(6)
    atom_a0 = topology.addAtom("A0", carbon, residue_a)
    atom_a1 = topology.addAtom("A1", carbon, residue_a)
    atom_b0 = topology.addAtom("B0", carbon, residue_b)
    atom_b1 = topology.addAtom("B1", carbon, residue_b)
    topology.addBond(atom_a0, atom_a1)
    topology.addBond(atom_b0, atom_b1)

    system = openmm.System()
    for _ in range(4):
        system.addParticle(12.0)

    tether = openmm.CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    tether.addGlobalParameter("k", 300.0)
    tether.addPerParticleParameter("x0")
    tether.addPerParticleParameter("y0")
    tether.addPerParticleParameter("z0")
    tether.addParticle(0, [0.00, 0.00, 0.00])
    tether.addParticle(1, [0.12, 0.00, 0.00])
    tether.addParticle(2, [0.28, 0.00, 0.00])
    tether.addParticle(3, [0.40, 0.00, 0.00])
    system.addForce(tether)

    bond_force = openmm.HarmonicBondForce()
    bond_force.addBond(0, 1, 0.12, 800.0)
    bond_force.addBond(2, 3, 0.12, 800.0)
    system.addForce(bond_force)

    integrator = openmm.LangevinMiddleIntegrator(310.0 * unit.kelvin, 10.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(
        [[0.00, 0.00, 0.00], [0.12, 0.00, 0.00], [0.28, 0.00, 0.00], [0.40, 0.00, 0.00]] * unit.nanometer
    )
    simulation.context.setVelocitiesToTemperature(310.0 * unit.kelvin, 7)
    return simulation, [2, 3], [0, 1]


def test_generate_window_centers_respects_config_spacing() -> None:
    """Chunk 14 gate: window centers should span the configured xi range inclusively."""

    centers = generate_window_centers(UmbrellaConfig(xi_min_nm=0.20, xi_max_nm=0.30, window_spacing_nm=0.05))

    assert np.allclose(centers, np.asarray([0.20, 0.25, 0.30], dtype=float))


def test_run_umbrella_window_restrains_sampling_near_target_center(tmp_path: Path) -> None:
    """Umbrella sampling should keep the COM distance centered near the requested window target."""

    simulation, pull_group_1, pull_group_2 = _make_umbrella_test_simulation()
    initial_positions = np.asarray(
        simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer),
        dtype=float,
    )
    initial_distance_nm = com_distance(initial_positions, np.full(4, 12.0, dtype=float), np.asarray(pull_group_1), np.asarray(pull_group_2))
    target_center_nm = initial_distance_nm - 0.02
    config = UmbrellaConfig(
        xi_min_nm=0.20,
        xi_max_nm=0.30,
        window_spacing_nm=0.05,
        spring_constant_kj_mol_nm2=2500.0,
        per_window_duration_ns=0.01,
        save_interval_ps=0.5,
    )

    result = run_umbrella_window(simulation, target_center_nm, config, window_id=1, output_dir=tmp_path)

    assert result["trajectory_path"].exists()
    assert result["xi_timeseries_path"].exists()
    assert result["xi_timeseries"].ndim == 1
    assert result["xi_timeseries"].size == 20
    assert abs(result["mean_xi_nm"] - target_center_nm) < 0.03
    assert result["std_xi_nm"] > 0.0
    assert result["mean_xi_nm"] < initial_distance_nm
    assert np.allclose(np.load(result["xi_timeseries_path"]), result["xi_timeseries"])


def test_run_umbrella_campaign_enforces_iv8_overlap(tmp_path: Path) -> None:
    """IV-8: adjacent umbrella windows should produce at least 10% xi-histogram overlap."""

    simulation, pull_group_1, pull_group_2 = _make_umbrella_test_simulation()
    system_xml_path = tmp_path / "system.xml"
    state_xml_path = tmp_path / "state.xml"
    system_xml_path.write_text(XmlSerializer.serialize(simulation.system), encoding="utf-8")
    state_xml_path.write_text(
        XmlSerializer.serialize(simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)),
        encoding="utf-8",
    )

    config = UmbrellaConfig(
        xi_min_nm=0.22,
        xi_max_nm=0.28,
        window_spacing_nm=0.03,
        spring_constant_kj_mol_nm2=5000.0,
        per_window_duration_ns=0.02,
        save_interval_ps=0.5,
    )

    results = run_umbrella_campaign(
        state_xml_path,
        system_xml_path,
        config,
        pull_group_1=pull_group_1,
        pull_group_2=pull_group_2,
        output_dir=tmp_path / "campaign",
    )

    assert len(results) == 3
    assert [result["window_id"] for result in results] == [1, 2, 3]
    assert all(result["trajectory_path"].exists() for result in results)
    assert all(result["xi_timeseries_path"].exists() for result in results)
    assert all(result["xi_timeseries"].shape == (40,) for result in results)
    assert results[0]["mean_xi_nm"] <= results[1]["mean_xi_nm"] <= results[2]["mean_xi_nm"]


def test_run_umbrella_window_rejects_non_positive_center(tmp_path: Path) -> None:
    """Public API should reject non-physical umbrella window centers."""

    simulation, _, _ = _make_umbrella_test_simulation()

    with pytest.raises(ValueError, match="window_center_nm must be positive"):
        run_umbrella_window(simulation, 0.0, UmbrellaConfig(xi_min_nm=0.20, xi_max_nm=0.30), 1, tmp_path)