"""Restraint-force utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import openmm


def _validate_system(system: openmm.System) -> int:
    """Validate the target system and return its particle count."""

    n_particles = system.getNumParticles()
    if n_particles <= 0:
        raise ValueError("system must contain at least one particle")
    return n_particles


def _validate_atom_indices(atom_indices: Sequence[int], n_particles: int, parameter_name: str) -> list[int]:
    """Validate a particle-index selection against the system size."""

    indices = [int(index) for index in atom_indices]
    if not indices:
        raise ValueError(f"{parameter_name} must be non-empty")
    if any(index < 0 or index >= n_particles for index in indices):
        raise ValueError(f"{parameter_name} must be within [0, N_atoms)")
    return indices


def create_positional_restraints(
    system: openmm.System,
    atom_indices: list[int],
    reference_positions: np.ndarray,
    force_constant_kj_mol_nm2: float,
) -> openmm.CustomExternalForce:
    """Create and attach positional restraints for selected atoms.

    Invariants: None.

    Args:
        system: OpenMM system to augment with the restraint force.
        atom_indices: Particle indices to restrain.
        reference_positions: Reference coordinates in nm. Shape: [N_restrained, 3].
        force_constant_kj_mol_nm2: Harmonic force constant in kJ/mol/nm^2.

    Returns:
        openmm.CustomExternalForce: Force object added to the system.
    """

    n_particles = _validate_system(system)
    validated_indices = _validate_atom_indices(atom_indices, n_particles, "atom_indices")

    positions_nm = np.asarray(reference_positions, dtype=float)
    if positions_nm.shape != (len(validated_indices), 3):
        raise ValueError("reference_positions must have shape [N_restrained, 3]")
    if force_constant_kj_mol_nm2 <= 0.0:
        raise ValueError("force_constant_kj_mol_nm2 must be positive")

    restraint = openmm.CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("k", float(force_constant_kj_mol_nm2))
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom_index, position in zip(validated_indices, positions_nm, strict=True):
        restraint.addParticle(atom_index, [float(position[0]), float(position[1]), float(position[2])])

    system.addForce(restraint)
    return restraint


def create_harmonic_distance_restraint(
    system: openmm.System,
    group_a_indices: list[int],
    group_b_indices: list[int],
    target_distance_nm: float,
    force_constant_kj_mol_nm2: float,
) -> openmm.CustomCentroidBondForce:
    """Create and attach a harmonic centroid-distance restraint between two groups.

    Invariants: None.

    Args:
        system: OpenMM system to augment with the restraint force.
        group_a_indices: First particle group. Shape: [N_a].
        group_b_indices: Second particle group. Shape: [N_b].
        target_distance_nm: Target centroid separation in nm.
        force_constant_kj_mol_nm2: Harmonic force constant in kJ/mol/nm^2.

    Returns:
        openmm.CustomCentroidBondForce: Force object added to the system.
    """

    n_particles = _validate_system(system)
    validated_group_a = _validate_atom_indices(group_a_indices, n_particles, "group_a_indices")
    validated_group_b = _validate_atom_indices(group_b_indices, n_particles, "group_b_indices")

    if target_distance_nm <= 0.0:
        raise ValueError("target_distance_nm must be positive")
    if force_constant_kj_mol_nm2 <= 0.0:
        raise ValueError("force_constant_kj_mol_nm2 must be positive")

    restraint = openmm.CustomCentroidBondForce(2, "0.5*k*(distance(g1, g2)-r0)^2")
    restraint.addGlobalParameter("k", float(force_constant_kj_mol_nm2))
    restraint.addGlobalParameter("r0", float(target_distance_nm))

    group_a = restraint.addGroup(validated_group_a)
    group_b = restraint.addGroup(validated_group_b)
    restraint.addBond([group_a, group_b], [])

    system.addForce(restraint)
    return restraint