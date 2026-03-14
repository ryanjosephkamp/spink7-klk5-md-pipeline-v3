"""Collective variable utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import numpy as np


def _validate_group_indices(group_indices: np.ndarray, n_atoms: int) -> np.ndarray:
    """Validate atom-group indices at the public API boundary."""

    indices = np.asarray(group_indices, dtype=int)
    if indices.ndim != 1:
        raise ValueError("group_indices must be a one-dimensional array")
    if indices.size == 0:
        raise ValueError("group_indices must be non-empty")
    if np.any(indices < 0) or np.any(indices >= n_atoms):
        raise ValueError("group_indices must be within [0, N_atoms)")
    return indices


def _validate_positions_and_masses(positions: np.ndarray, masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Validate shared array inputs for collective variable calculations."""

    positions_array = np.asarray(positions, dtype=float)
    masses_array = np.asarray(masses, dtype=float)

    if positions_array.ndim != 2 or positions_array.shape[1] != 3:
        raise ValueError("positions must have shape [N_atoms, 3]")
    if masses_array.ndim != 1:
        raise ValueError("masses must have shape [N_atoms]")
    if positions_array.shape[0] != masses_array.shape[0]:
        raise ValueError("positions and masses must describe the same number of atoms")
    if np.any(masses_array <= 0.0):
        raise ValueError("masses must be strictly positive")

    return positions_array, masses_array


def com_vector(
    positions: np.ndarray,
    masses: np.ndarray,
    group_indices: np.ndarray,
) -> np.ndarray:
    """Compute the center-of-mass vector for a selected atom group.

    Invariants: None.

    Args:
        positions: Cartesian coordinates in nm. Shape: [N_atoms, 3].
        masses: Atomic masses in amu. Shape: [N_atoms].
        group_indices: Atom indices for the COM group. Shape: [N_group].

    Returns:
        np.ndarray: Center-of-mass vector in nm. Shape: [3].
    """

    positions_array, masses_array = _validate_positions_and_masses(positions, masses)
    indices = _validate_group_indices(group_indices, positions_array.shape[0])

    group_masses = masses_array[indices]
    group_positions = positions_array[indices]
    total_mass = np.sum(group_masses)

    return np.sum(group_positions * group_masses[:, np.newaxis], axis=0) / total_mass


def com_distance(
    positions: np.ndarray,
    masses: np.ndarray,
    group_a_indices: np.ndarray,
    group_b_indices: np.ndarray,
) -> float:
    """Compute the center-of-mass distance between two selected atom groups.

    Invariants: None.

    Args:
        positions: Cartesian coordinates in nm. Shape: [N_atoms, 3].
        masses: Atomic masses in amu. Shape: [N_atoms].
        group_a_indices: Atom indices for group A. Shape: [N_a].
        group_b_indices: Atom indices for group B. Shape: [N_b].

    Returns:
        float: Scalar COM distance in nm.
    """

    com_a = com_vector(positions, masses, group_a_indices)
    com_b = com_vector(positions, masses, group_b_indices)
    return float(np.linalg.norm(com_a - com_b))
