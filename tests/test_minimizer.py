"""Tests for energy-minimization utilities."""

from __future__ import annotations

import openmm
import pytest
from openmm import unit
from openmm.app import Element, Simulation, Topology

from src.config import MinimizationConfig
from src.simulate.minimizer import minimize_energy


def _make_stretched_bond_simulation() -> Simulation:
    """Create a tiny bonded system with a deliberately strained geometry."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MOL", chain)
    atom_a = topology.addAtom("A1", Element.getByAtomicNumber(6), residue)
    atom_b = topology.addAtom("A2", Element.getByAtomicNumber(6), residue)
    topology.addBond(atom_a, atom_b)

    system = openmm.System()
    system.addParticle(12.0)
    system.addParticle(12.0)

    bond_force = openmm.HarmonicBondForce()
    bond_force.addBond(0, 1, 0.10, 5000.0)
    system.addForce(bond_force)

    integrator = openmm.VerletIntegrator(0.001)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions([[0.0, 0.0, 0.0], [0.50, 0.0, 0.0]] * unit.nanometer)
    return simulation


def test_minimize_energy_reduces_potential_energy() -> None:
    """IV-1: post-minimization potential energy must be strictly lower than initial energy."""

    simulation = _make_stretched_bond_simulation()
    config = MinimizationConfig(max_iterations=200, tolerance_kj_mol_nm=0.1)

    result = minimize_energy(simulation, config)

    assert result["final_energy_kj_mol"] < result["initial_energy_kj_mol"]
    assert result["n_steps"] == config.max_iterations


def test_minimize_energy_updates_simulation_context_energy() -> None:
    """The simulation context should reflect the lowered post-minimization energy."""

    simulation = _make_stretched_bond_simulation()
    config = MinimizationConfig(max_iterations=200, tolerance_kj_mol_nm=0.1)

    initial_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    result = minimize_energy(simulation, config)
    final_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )

    assert final_energy < initial_energy
    assert final_energy == pytest.approx(result["final_energy_kj_mol"], abs=1e-9)


def test_minimize_energy_rejects_non_positive_tolerance() -> None:
    """Public API should reject non-physical minimization tolerances."""

    simulation = _make_stretched_bond_simulation()

    with pytest.raises(ValueError, match="config.tolerance_kj_mol_nm must be positive"):
        minimize_energy(simulation, MinimizationConfig(tolerance_kj_mol_nm=0.0))