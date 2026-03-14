"""Tests for steered molecular dynamics utilities."""

from __future__ import annotations

import numpy as np
import openmm
import pytest
from openmm import unit
from openmm.app import Element, Simulation, Topology

from src.config import SMDConfig
from src.physics.collective_variables import com_distance
from src.simulate.smd import run_smd_replicate


def _make_smd_test_simulation() -> tuple[Simulation, list[int], list[int], np.ndarray, float]:
    """Build a tiny anchored two-group system for fast SMD validation."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("SMD", chain)
    carbon = Element.getByAtomicNumber(6)
    atom_0 = topology.addAtom("A0", carbon, residue)
    atom_1 = topology.addAtom("A1", carbon, residue)
    atom_2 = topology.addAtom("A2", carbon, residue)
    topology.addBond(atom_0, atom_1)
    topology.addBond(atom_1, atom_2)

    system = openmm.System()
    for _ in range(3):
        system.addParticle(12.0)

    tether = openmm.CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    tether.addGlobalParameter("k", 5000.0)
    tether.addPerParticleParameter("x0")
    tether.addPerParticleParameter("y0")
    tether.addPerParticleParameter("z0")
    tether.addParticle(0, [0.0, 0.0, 0.0])
    tether.addParticle(1, [0.2, 0.0, 0.0])
    system.addForce(tether)

    bond_force = openmm.HarmonicBondForce()
    bond_force.addBond(0, 1, 0.2, 2000.0)
    bond_force.addBond(1, 2, 0.2, 2000.0)
    system.addForce(bond_force)

    integrator = openmm.LangevinMiddleIntegrator(310.0 * unit.kelvin, 5.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.4, 0.0, 0.0]] * unit.nanometer
    )
    simulation.context.setVelocitiesToTemperature(310.0 * unit.kelvin, 11)

    initial_positions = np.asarray(
        simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer),
        dtype=float,
    )
    masses = np.asarray([12.0, 12.0, 12.0], dtype=float)
    initial_distance_nm = com_distance(initial_positions, masses, np.asarray([2]), np.asarray([0, 1]))
    return simulation, [2], [0, 1], np.asarray([1.0, 0.0, 0.0], dtype=float), initial_distance_nm


def test_run_smd_replicate_displaces_coordinate_and_logs_work(tmp_path) -> None:
    """SMD should increase the COM distance and record non-equilibrium work consistently."""

    simulation, pull_group_1, pull_group_2, pull_direction, initial_distance_nm = _make_smd_test_simulation()
    config = SMDConfig(
        spring_constant_kj_mol_nm2=4000.0,
        pulling_velocity_nm_per_ps=0.02,
        pull_distance_nm=0.10,
        save_interval_ps=0.25,
    )

    result = run_smd_replicate(
        simulation,
        config,
        replicate_id=1,
        output_dir=tmp_path,
        pull_group_1=pull_group_1,
        pull_group_2=pull_group_2,
        pull_direction=pull_direction,
    )

    assert result["trajectory_path"].exists()
    assert result["work_timeseries_path"].exists()
    assert result["force_timeseries_path"].exists()
    assert result["xi_timeseries_path"].exists()
    assert result["work_timeseries"].shape[1] == 2
    assert result["force_timeseries"].shape == result["work_timeseries"].shape
    assert result["xi_timeseries"].shape == result["work_timeseries"].shape
    assert np.all(np.diff(result["work_timeseries"][:, 0]) > 0.0)
    assert np.all(np.diff(result["xi_timeseries"][:, 0]) > 0.0)
    assert result["xi_timeseries"][-1, 1] > initial_distance_nm
    assert result["work_timeseries"][-1, 1] == pytest.approx(result["total_work_kj_mol"], rel=1e-9, abs=1e-9)
    assert np.any(np.abs(result["force_timeseries"][:, 1]) > 0.0)

    work_csv = np.loadtxt(result["work_timeseries_path"], delimiter=",", skiprows=1)
    force_csv = np.loadtxt(result["force_timeseries_path"], delimiter=",", skiprows=1)
    xi_csv = np.loadtxt(result["xi_timeseries_path"], delimiter=",", skiprows=1)
    assert work_csv.shape == result["work_timeseries"].shape
    assert force_csv.shape == result["force_timeseries"].shape
    assert xi_csv.shape == result["xi_timeseries"].shape
    assert np.allclose(work_csv, result["work_timeseries"])
    assert np.allclose(force_csv, result["force_timeseries"])
    assert np.allclose(xi_csv, result["xi_timeseries"])


def test_run_smd_replicate_rejects_misaligned_pull_direction(tmp_path) -> None:
    """Public API should reject pull directions that oppose the initial COM axis."""

    simulation, pull_group_1, pull_group_2, _, _ = _make_smd_test_simulation()

    with pytest.raises(ValueError, match="pull_direction must point from pull_group_2 toward pull_group_1"):
        run_smd_replicate(
            simulation,
            SMDConfig(pull_distance_nm=0.05, save_interval_ps=0.5),
            replicate_id=1,
            output_dir=tmp_path,
            pull_group_1=pull_group_1,
            pull_group_2=pull_group_2,
            pull_direction=np.asarray([-1.0, 0.0, 0.0], dtype=float),
        )