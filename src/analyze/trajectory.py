"""Trajectory I/O, frame selection, and alignment utilities."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import mdtraj as md
import numpy as np


def _validate_trajectory_path(trajectory_path: Path) -> Path:
    """Validate an on-disk trajectory path at the public API boundary."""

    path = Path(trajectory_path)
    if not path.exists():
        raise FileNotFoundError(f"trajectory file not found: {path}")
    if path.suffix.lower() not in {".dcd", ".xtc", ".h5", ".nc", ".trr"}:
        raise ValueError("trajectory_path must have a supported trajectory extension")
    return path


def _validate_topology_path(topology_path: Path) -> Path:
    """Validate a topology file path for trajectory loading."""

    path = Path(topology_path)
    if not path.exists():
        raise FileNotFoundError(f"topology file not found: {path}")
    if path.suffix.lower() not in {".pdb", ".cif", ".gro", ".prmtop", ".parm7"}:
        raise ValueError("topology_path must have a supported topology extension")
    return path


def _validate_stride(stride: int) -> int:
    """Validate trajectory stride parameters."""

    validated_stride = int(stride)
    if validated_stride <= 0:
        raise ValueError("stride must be positive")
    return validated_stride


def _validate_chunk_size(chunk_size: int) -> int:
    """Validate chunk sizes for iterloaded trajectories."""

    validated_chunk_size = int(chunk_size)
    if validated_chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    return validated_chunk_size


def load_trajectory(
    trajectory_path: Path,
    topology_path: Path,
    stride: int = 1,
) -> md.Trajectory:
    """Load a trajectory from disk into an MDTraj trajectory object.

    Invariants: None.

    Args:
        trajectory_path: On-disk trajectory file path.
        topology_path: On-disk topology file path.
        stride: Frame stride for loading.

    Returns:
        md.Trajectory: Loaded coordinates. Shape: [N_frames, N_atoms, 3].
    """

    trajectory_file = _validate_trajectory_path(trajectory_path)
    topology_file = _validate_topology_path(topology_path)
    validated_stride = _validate_stride(stride)
    return md.load(str(trajectory_file), top=str(topology_file), stride=validated_stride)


def iterload_trajectory(
    trajectory_path: Path,
    topology_path: Path,
    chunk_size: int = 100,
    stride: int = 1,
) -> Iterator[md.Trajectory]:
    """Iteratively stream a trajectory from disk in bounded chunks.

    Invariants: None.

    Args:
        trajectory_path: On-disk trajectory file path.
        topology_path: On-disk topology file path.
        chunk_size: Number of frames per chunk.
        stride: Frame stride for loading.

    Yields:
        md.Trajectory: Chunked trajectory slice. Shape: [N_chunk_frames, N_atoms, 3].
    """

    trajectory_file = _validate_trajectory_path(trajectory_path)
    topology_file = _validate_topology_path(topology_path)
    validated_chunk_size = _validate_chunk_size(chunk_size)
    validated_stride = _validate_stride(stride)

    return md.iterload(
        str(trajectory_file),
        top=str(topology_file),
        chunk=validated_chunk_size,
        stride=validated_stride,
    )


def select_frames(trajectory: md.Trajectory, frame_indices: np.ndarray) -> md.Trajectory:
    """Select an ordered subset of frames from a trajectory.

    Invariants: None.

    Args:
        trajectory: MDTraj trajectory. Shape: [N_frames, N_atoms, 3].
        frame_indices: Selected frame indices. Shape: [N_selected].

    Returns:
        md.Trajectory: Selected frame subset. Shape: [N_selected, N_atoms, 3].
    """

    indices = np.asarray(frame_indices, dtype=int)
    if indices.ndim != 1:
        raise ValueError("frame_indices must be a one-dimensional array")
    if indices.size == 0:
        raise ValueError("frame_indices must be non-empty")
    if np.any(indices < 0) or np.any(indices >= trajectory.n_frames):
        raise ValueError("frame_indices must be within [0, N_frames)")
    return trajectory.slice(indices, copy=True)


def align_trajectory(
    trajectory: md.Trajectory,
    reference: md.Trajectory,
    atom_selection: str = "backbone",
) -> md.Trajectory:
    """Align a trajectory to a reference structure using an MDTraj atom selection.

    Invariants: None.

    Args:
        trajectory: MDTraj trajectory. Shape: [N_frames, N_atoms, 3].
        reference: MDTraj reference trajectory. Shape: [N_ref_frames, N_atoms, 3].
        atom_selection: MDTraj atom-selection query used for superposition.

    Returns:
        md.Trajectory: Aligned trajectory. Shape: [N_frames, N_atoms, 3].
    """

    if reference.n_frames < 1:
        raise ValueError("reference must contain at least one frame")

    atom_indices = trajectory.topology.select(atom_selection)
    reference_indices = reference.topology.select(atom_selection)
    if atom_indices.size == 0 or reference_indices.size == 0:
        raise ValueError("atom_selection must match at least one atom in both trajectories")
    if atom_indices.size != reference_indices.size:
        raise ValueError("atom_selection must produce the same number of atoms in trajectory and reference")

    aligned = trajectory[:]
    aligned.superpose(reference, frame=0, atom_indices=atom_indices, ref_atom_indices=reference_indices)
    return aligned