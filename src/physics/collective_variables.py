"""Collective variable utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field

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
    box_lengths: np.ndarray | None = None,
) -> np.ndarray:
    """Compute the center-of-mass vector for a selected atom group.

    When box_lengths is provided, uses the angular mean method [Bai & Breen,
    J. Comput. Phys. 199, 737 (2004)] to correctly handle atoms that wrap
    across periodic boundaries. When box_lengths is None, falls back to
    naive Cartesian averaging (suitable for non-periodic or unwrapped systems).

    Invariants: None.

    Args:
        positions: Cartesian coordinates in nm. Shape: [N_atoms, 3].
        masses: Atomic masses in amu. Shape: [N_atoms].
        group_indices: Atom indices for the COM group. Shape: [N_group].
        box_lengths: Periodic box edge lengths in nm. Shape: [3].
            If None, no periodic correction is applied.

    Returns:
        np.ndarray: Center-of-mass vector in nm. Shape: [3].
    """

    positions_array, masses_array = _validate_positions_and_masses(positions, masses)
    indices = _validate_group_indices(group_indices, positions_array.shape[0])

    group_masses = masses_array[indices]
    group_positions = positions_array[indices]
    total_mass = np.sum(group_masses)

    if box_lengths is not None:
        box = np.asarray(box_lengths, dtype=float)
        if box.shape != (3,):
            raise ValueError("box_lengths must have shape [3]")
        if np.any(box <= 0.0):
            raise ValueError("box_lengths must be strictly positive")

        # Angular mean method (Bai & Breen, 2004)
        theta = 2.0 * np.pi * group_positions / box  # [N_group, 3]
        weights = group_masses[:, np.newaxis]          # [N_group, 1]
        sin_sum = np.sum(weights * np.sin(theta), axis=0)  # [3]
        cos_sum = np.sum(weights * np.cos(theta), axis=0)  # [3]
        theta_bar = np.arctan2(sin_sum, cos_sum)            # [3]
        com = (box / (2.0 * np.pi)) * (theta_bar % (2.0 * np.pi))
        return com

    # Fallback: naive arithmetic mean (non-periodic systems)
    return np.sum(group_positions * group_masses[:, np.newaxis], axis=0) / total_mass


def com_distance(
    positions: np.ndarray,
    masses: np.ndarray,
    group_a_indices: np.ndarray,
    group_b_indices: np.ndarray,
    box_lengths: np.ndarray | None = None,
) -> float:
    """Compute the center-of-mass distance between two selected atom groups.

    When box_lengths is provided, both COM vectors are computed using the
    PBC-aware angular mean method, and the minimum image convention is applied
    to the inter-COM displacement vector.

    Invariants: None.

    Args:
        positions: Cartesian coordinates in nm. Shape: [N_atoms, 3].
        masses: Atomic masses in amu. Shape: [N_atoms].
        group_a_indices: Atom indices for group A. Shape: [N_a].
        group_b_indices: Atom indices for group B. Shape: [N_b].
        box_lengths: Periodic box edge lengths in nm. Shape: [3].
            If None, no periodic correction is applied.

    Returns:
        float: Scalar COM distance in nm.
    """

    com_a = com_vector(positions, masses, group_a_indices, box_lengths=box_lengths)
    com_b = com_vector(positions, masses, group_b_indices, box_lengths=box_lengths)
    delta = com_a - com_b

    if box_lengths is not None:
        box = np.asarray(box_lengths, dtype=float)
        # Minimum image convention: wrap displacement to [-L/2, +L/2]
        delta = delta - box * np.round(delta / box)

    return float(np.linalg.norm(delta))


def com_angle(
    positions: np.ndarray,
    masses: np.ndarray,
    group_a_indices: np.ndarray,
    group_b_indices: np.ndarray,
    reference_direction: np.ndarray,
    box_lengths: np.ndarray | None = None,
) -> float:
    """Compute the angular deviation of the inter-COM axis from a reference direction.

    .. math::
        \\theta = \\arccos\\!\\left(\\text{clamp}\\!\\left(
        \\frac{(\\mathbf{R}_{COM}^A - \\mathbf{R}_{COM}^B) \\cdot \\hat{\\mathbf{n}}_{ref}}
        {\\|\\mathbf{R}_{COM}^A - \\mathbf{R}_{COM}^B\\|}, -1, 1\\right)\\right)

    Args:
        positions: Cartesian coordinates in nm. Shape: [N_atoms, 3].
        masses: Atomic masses in amu. Shape: [N_atoms].
        group_a_indices: Atom indices for group A. Shape: [N_a].
        group_b_indices: Atom indices for group B. Shape: [N_b].
        reference_direction: Reference direction vector. Shape: [3]. Must have non-zero norm.
        box_lengths: Periodic box edge lengths in nm. Shape: [3].
            If None, no periodic correction is applied.

    Returns:
        float: Angle theta in radians, in [0, pi].
    """
    positions_array, masses_array = _validate_positions_and_masses(positions, masses)
    n_atoms = positions_array.shape[0]
    _validate_group_indices(np.asarray(group_a_indices, dtype=int), n_atoms)
    _validate_group_indices(np.asarray(group_b_indices, dtype=int), n_atoms)

    ref = np.asarray(reference_direction, dtype=float)
    if ref.shape != (3,):
        raise ValueError("reference_direction must have shape [3]")
    if not np.all(np.isfinite(ref)):
        raise ValueError("reference_direction must contain only finite values")
    ref_norm = float(np.linalg.norm(ref))
    if ref_norm <= 0.0:
        raise ValueError("reference_direction must have non-zero norm")
    ref_hat = ref / ref_norm

    com_a = com_vector(positions_array, masses_array, np.asarray(group_a_indices, dtype=int), box_lengths=box_lengths)
    com_b = com_vector(positions_array, masses_array, np.asarray(group_b_indices, dtype=int), box_lengths=box_lengths)
    delta = com_a - com_b

    if box_lengths is not None:
        box = np.asarray(box_lengths, dtype=float)
        delta = delta - box * np.round(delta / box)

    delta_norm = float(np.linalg.norm(delta))
    if delta_norm <= 0.0:
        return 0.0

    cos_theta = float(np.dot(delta / delta_norm, ref_hat))
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return float(np.arccos(cos_theta))


def contact_fraction(
    positions: np.ndarray,
    contact_pairs: np.ndarray,
    contact_reference_distances_nm: np.ndarray,
    switch_exponent: int = 12,
    distance_factor: float = 1.5,
) -> float:
    """Compute the differentiable interface contact fraction Q.

    Uses a soft switching function to count native contacts:

    .. math::
        Q = \\frac{1}{N_{pairs}} \\sum_{k} \\frac{1}{1 + (r_k / r_0^{(k)})^n}

    where :math:`r_0^{(k)} = \\text{distance\\_factor} \\times r_{ref}^{(k)}`.

    Args:
        positions: Cartesian coordinates in nm. Shape: [N_atoms, 3].
        contact_pairs: Native contact atom index pairs. Shape: [N_pairs, 2].
        contact_reference_distances_nm: Reference distances in nm. Shape: [N_pairs].
        switch_exponent: Exponent n for the switching function (default 12).
        distance_factor: Multiplier for reference distances to set switching midpoint (default 1.5).

    Returns:
        float: Contact fraction Q in [0, 1].
    """
    positions_array = np.asarray(positions, dtype=float)
    if positions_array.ndim != 2 or positions_array.shape[1] != 3:
        raise ValueError("positions must have shape [N_atoms, 3]")

    pairs = np.asarray(contact_pairs, dtype=int)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError("contact_pairs must have shape [N_pairs, 2]")
    if pairs.shape[0] == 0:
        raise ValueError("contact_pairs must be non-empty")

    n_atoms = positions_array.shape[0]
    if np.any(pairs < 0) or np.any(pairs >= n_atoms):
        raise ValueError("contact_pairs indices must be within [0, N_atoms)")

    ref_dists = np.asarray(contact_reference_distances_nm, dtype=float)
    if ref_dists.ndim != 1 or ref_dists.shape[0] != pairs.shape[0]:
        raise ValueError("contact_reference_distances_nm must have shape [N_pairs]")
    if np.any(ref_dists <= 0.0):
        raise ValueError("contact_reference_distances_nm must be strictly positive")

    # Compute pairwise distances
    r_ij = np.linalg.norm(positions_array[pairs[:, 0]] - positions_array[pairs[:, 1]], axis=1)
    r_0 = distance_factor * ref_dists
    ratio = r_ij / r_0
    switching = 1.0 / (1.0 + ratio ** switch_exponent)
    return float(np.mean(switching))


@dataclass(frozen=True)
class CollectiveVariableSpec:
    """Specification for a collective variable to compute during simulation."""

    name: str  # "com_distance", "com_angle", "contact_fraction"
    reference_direction: np.ndarray | None = field(default=None, repr=False)
    contact_pairs: np.ndarray | None = field(default=None, repr=False)
    contact_reference_distances_nm: np.ndarray | None = field(default=None, repr=False)
    switch_exponent: int = 12
    distance_factor: float = 1.5


def compute_cv(
    spec: CollectiveVariableSpec,
    positions: np.ndarray,
    masses: np.ndarray,
    group_a_indices: np.ndarray,
    group_b_indices: np.ndarray,
    box_lengths: np.ndarray | None = None,
) -> float:
    """Dispatch to the appropriate CV function based on the spec name.

    Args:
        spec: Collective variable specification.
        positions: Cartesian coordinates in nm. Shape: [N_atoms, 3].
        masses: Atomic masses in amu. Shape: [N_atoms].
        group_a_indices: Atom indices for group A.
        group_b_indices: Atom indices for group B.
        box_lengths: Periodic box edge lengths in nm. Shape: [3].

    Returns:
        float: Computed CV value.
    """
    if spec.name == "com_distance":
        return com_distance(positions, masses, group_a_indices, group_b_indices, box_lengths=box_lengths)
    if spec.name == "com_angle":
        if spec.reference_direction is None:
            raise ValueError("CollectiveVariableSpec for 'com_angle' requires reference_direction")
        return com_angle(positions, masses, group_a_indices, group_b_indices, spec.reference_direction, box_lengths=box_lengths)
    if spec.name == "contact_fraction":
        if spec.contact_pairs is None or spec.contact_reference_distances_nm is None:
            raise ValueError("CollectiveVariableSpec for 'contact_fraction' requires contact_pairs and contact_reference_distances_nm")
        return contact_fraction(positions, spec.contact_pairs, spec.contact_reference_distances_nm, spec.switch_exponent, spec.distance_factor)
    raise ValueError(f"Unknown collective variable: {spec.name}")


def default_com_distance_spec() -> CollectiveVariableSpec:
    """Return the default COM-distance CV specification."""
    return CollectiveVariableSpec(name="com_distance")
