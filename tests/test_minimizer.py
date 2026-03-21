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
    assert "converged" in result
    assert isinstance(result["converged"], bool)
    assert "max_force_kj_mol_nm" in result
    assert result["max_force_kj_mol_nm"] >= 0.0
    assert "energy_reduction" in result
    assert result["energy_reduction"] > 0.0


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
    assert "converged" in result


def test_minimize_energy_rejects_non_positive_tolerance() -> None:
    """Public API should reject non-physical minimization tolerances."""

    simulation = _make_stretched_bond_simulation()

    with pytest.raises(ValueError, match="config.tolerance_kj_mol_nm must be positive"):
        minimize_energy(simulation, MinimizationConfig(tolerance_kj_mol_nm=0.0))


def test_minimize_energy_reports_convergence_status() -> None:
    """Converged minimization should report converged=True and force below tolerance."""

    simulation = _make_stretched_bond_simulation()
    config = MinimizationConfig(max_iterations=500, tolerance_kj_mol_nm=10.0)

    result = minimize_energy(simulation, config)

    assert result["converged"] is True
    assert result["max_force_kj_mol_nm"] < config.tolerance_kj_mol_nm
    assert result["energy_reduction"] > 0


def test_minimize_energy_warns_on_non_convergence() -> None:
    """Non-converged minimization should emit a warning and report converged=False."""

    import warnings

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MOL", chain)
    a1 = topology.addAtom("A1", Element.getByAtomicNumber(6), residue)
    a2 = topology.addAtom("A2", Element.getByAtomicNumber(6), residue)
    topology.addBond(a1, a2)

    system = openmm.System()
    system.addParticle(12.0)
    system.addParticle(12.0)
    bond_force = openmm.HarmonicBondForce()
    bond_force.addBond(0, 1, 0.10, 500000.0)  # extremely stiff
    system.addForce(bond_force)

    integrator = openmm.VerletIntegrator(0.001)
    platform = openmm.Platform.getPlatformByName("CPU")
    sim = Simulation(topology, system, integrator, platform)
    sim.context.setPositions([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] * unit.nanometer)

    # Very tight tolerance with very few iterations — should not converge
    config = MinimizationConfig(max_iterations=1, tolerance_kj_mol_nm=0.001)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = minimize_energy(sim, config)

    assert result["final_energy_kj_mol"] < result["initial_energy_kj_mol"], "IV-1 must still hold"
    assert "converged" in result
    assert "max_force_kj_mol_nm" in result
    assert "energy_reduction" in result