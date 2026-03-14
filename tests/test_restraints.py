"""Tests for restraint-force construction utilities."""

from __future__ import annotations

import numpy as np
import openmm
import pytest
from openmm import unit

from src.physics.restraints import create_harmonic_distance_restraint, create_positional_restraints


def _make_system(n_particles: int = 3) -> openmm.System:
    """Create a minimal OpenMM system with identical particle masses."""

    system = openmm.System()
    for _ in range(n_particles):
        system.addParticle(12.0)
    return system


def _restraint_energy_kj_mol(system: openmm.System, positions_nm: np.ndarray) -> float:
    """Evaluate the potential energy of a restrained system at the given positions."""

    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator)
    try:
        context.setPositions(positions_nm * unit.nanometer)
        state = context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    finally:
        del context
        del integrator


def test_create_positional_restraints_adds_force_and_is_zero_at_reference_positions() -> None:
    """Chunk 9 gate: positional restraints should attach to the system and vanish at the reference."""

    system = _make_system(2)
    reference_positions = np.array([[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]], dtype=float)

    restraint = create_positional_restraints(
        system,
        atom_indices=[0, 1],
        reference_positions=reference_positions,
        force_constant_kj_mol_nm2=1000.0,
    )

    displaced_positions = np.array([[0.1, 0.0, 0.0], [0.2, 0.1, 0.0]], dtype=float)
    reference_energy = _restraint_energy_kj_mol(system, reference_positions)
    displaced_energy = _restraint_energy_kj_mol(system, displaced_positions)

    assert system.getNumForces() == 1
    assert isinstance(restraint, openmm.CustomExternalForce)
    assert restraint.getNumParticles() == 2
    assert reference_energy == pytest.approx(0.0, abs=1e-10)
    assert displaced_energy > 0.0


def test_create_harmonic_distance_restraint_adds_force_and_penalizes_distance_deviation() -> None:
    """Distance restraints should be zero at the target centroid separation and positive away from it."""

    system = _make_system(3)
    restraint = create_harmonic_distance_restraint(
        system,
        group_a_indices=[0],
        group_b_indices=[1, 2],
        target_distance_nm=1.0,
        force_constant_kj_mol_nm2=500.0,
    )

    target_positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    off_target_positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.4, 0.0, 0.0],
            [1.4, 0.0, 0.0],
        ],
        dtype=float,
    )

    target_energy = _restraint_energy_kj_mol(system, target_positions)
    off_target_energy = _restraint_energy_kj_mol(system, off_target_positions)

    assert system.getNumForces() == 1
    assert isinstance(restraint, openmm.CustomCentroidBondForce)
    assert restraint.getNumGroups() == 2
    assert restraint.getNumBonds() == 1
    assert target_energy == pytest.approx(0.0, abs=1e-10)
    assert off_target_energy > 0.0


def test_create_positional_restraints_rejects_shape_mismatch() -> None:
    """Public API should reject reference-position arrays with the wrong shape."""

    system = _make_system(2)

    with pytest.raises(ValueError, match=r"reference_positions must have shape \[N_restrained, 3\]"):
        create_positional_restraints(
            system,
            atom_indices=[0, 1],
            reference_positions=np.array([[0.0, 0.0, 0.0]], dtype=float),
            force_constant_kj_mol_nm2=1000.0,
        )


def test_create_harmonic_distance_restraint_rejects_invalid_group_indices() -> None:
    """Public API should reject empty or out-of-range particle groups."""

    system = _make_system(2)

    with pytest.raises(ValueError, match=r"group_a_indices must be within \[0, N_atoms\)"):
        create_harmonic_distance_restraint(
            system,
            group_a_indices=[2],
            group_b_indices=[1],
            target_distance_nm=1.0,
            force_constant_kj_mol_nm2=1000.0,
        )

    with pytest.raises(ValueError, match="group_b_indices must be non-empty"):
        create_harmonic_distance_restraint(
            system,
            group_a_indices=[0],
            group_b_indices=[],
            target_distance_nm=1.0,
            force_constant_kj_mol_nm2=1000.0,
        )