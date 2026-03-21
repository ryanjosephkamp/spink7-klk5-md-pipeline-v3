"""Structural observables for molecular dynamics trajectories."""

from __future__ import annotations

import mdtraj as md
import numpy as np

from src.analyze.trajectory import align_trajectory, unwrap_trajectory


def _selection_indices(trajectory: md.Trajectory, atom_selection: str) -> np.ndarray:
    """Resolve and validate an MDTraj atom selection."""

    indices = np.asarray(trajectory.topology.select(atom_selection), dtype=int)
    if indices.ndim != 1 or indices.size == 0:
        raise ValueError("atom_selection must match at least one atom")
    return indices


def compute_rmsd(
    trajectory: md.Trajectory,
    reference: md.Trajectory,
    atom_selection: str = "backbone",
    unwrap: bool = True,
) -> np.ndarray:
    """Compute per-frame RMSD relative to a reference structure.

    Invariants: None.

    Args:
        trajectory: MD trajectory. Shape: [N_frames, N_atoms, 3].
        reference: Reference structure. Shape: [N_ref_frames, N_atoms, 3].
        atom_selection: MDTraj atom selection string.
        unwrap: If True, apply PBC unwrapping before computation.

    Returns:
        np.ndarray: RMSD values in nm. Shape: [N_frames].
    """

    if unwrap:
        trajectory = unwrap_trajectory(trajectory)
        reference = unwrap_trajectory(reference)

    if reference.n_frames < 1:
        raise ValueError("reference must contain at least one frame")

    atom_indices = _selection_indices(trajectory, atom_selection)
    reference_indices = _selection_indices(reference, atom_selection)
    if atom_indices.size != reference_indices.size:
        raise ValueError("atom_selection must produce the same number of atoms in trajectory and reference")

    aligned = align_trajectory(trajectory, reference, atom_selection=atom_selection)
    aligned_coordinates = aligned.xyz[:, atom_indices, :]
    reference_coordinates = reference.xyz[0, reference_indices, :]
    squared_displacements = np.sum((aligned_coordinates - reference_coordinates[np.newaxis, :, :]) ** 2, axis=2)
    rmsd_nm = np.sqrt(np.mean(squared_displacements, axis=1, dtype=float))
    return np.asarray(rmsd_nm, dtype=float)


def compute_rmsf(
    trajectory: md.Trajectory,
    atom_selection: str = "name CA",
    unwrap: bool = True,
) -> np.ndarray:
    """Compute per-atom RMSF over an aligned trajectory.

    Invariants: None.

    Args:
        trajectory: MD trajectory. Shape: [N_frames, N_atoms, 3].
        atom_selection: MDTraj atom selection string.
        unwrap: If True, apply PBC unwrapping before computation.

    Returns:
        np.ndarray: RMSF values in nm. Shape: [N_selected_atoms].
    """

    if unwrap:
        trajectory = unwrap_trajectory(trajectory)

    atom_indices = _selection_indices(trajectory, atom_selection)
    aligned = align_trajectory(trajectory, trajectory[0], atom_selection=atom_selection)
    rmsf_nm = md.rmsf(aligned, aligned[0], frame=0, atom_indices=atom_indices)
    return np.asarray(rmsf_nm, dtype=float)


def compute_radius_of_gyration(
    trajectory: md.Trajectory,
    unwrap: bool = True,
) -> np.ndarray:
    """Compute the radius of gyration for each trajectory frame.

    Invariants: None.

    Args:
        trajectory: MD trajectory. Shape: [N_frames, N_atoms, 3].
        unwrap: If True, apply PBC unwrapping before computation.

    Returns:
        np.ndarray: Radius-of-gyration values in nm. Shape: [N_frames].
    """

    if unwrap:
        trajectory = unwrap_trajectory(trajectory)

    rg_nm = md.compute_rg(trajectory)
    return np.asarray(rg_nm, dtype=float)


def compute_sasa(
    trajectory: md.Trajectory,
    probe_radius_nm: float = 0.14,
    unwrap: bool = True,
) -> np.ndarray:
    """Compute per-frame solvent accessible surface area.

    Invariants: None.

    Args:
        trajectory: MD trajectory. Shape: [N_frames, N_atoms, 3].
        probe_radius_nm: Probe radius in nm.
        unwrap: If True, apply PBC unwrapping before computation.

    Returns:
        np.ndarray: Per-frame SASA values in nm^2. Shape: [N_frames].
    """

    if unwrap:
        trajectory = unwrap_trajectory(trajectory)

    if probe_radius_nm <= 0.0:
        raise ValueError("probe_radius_nm must be positive")

    atom_sasa_nm2 = md.shrake_rupley(trajectory, probe_radius=probe_radius_nm, mode="atom")
    total_sasa_nm2 = np.zeros(trajectory.n_frames, dtype=float)
    total_sasa_nm2[:] = np.sum(atom_sasa_nm2, axis=1, dtype=float)
    return total_sasa_nm2