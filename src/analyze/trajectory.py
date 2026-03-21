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
    unwrap: bool = True,
) -> md.Trajectory:
    """Load a trajectory from disk into an MDTraj trajectory object.

    Invariants: None.

    Args:
        trajectory_path: On-disk trajectory file path.
        topology_path: On-disk topology file path.
        stride: Frame stride for loading.
        unwrap: If True, apply PBC unwrapping via image_molecules.

    Returns:
        md.Trajectory: Loaded coordinates. Shape: [N_frames, N_atoms, 3].
    """

    trajectory_file = _validate_trajectory_path(trajectory_path)
    topology_file = _validate_topology_path(topology_path)
    validated_stride = _validate_stride(stride)
    trajectory = md.load(str(trajectory_file), top=str(topology_file), stride=validated_stride)
    if unwrap:
        trajectory = unwrap_trajectory(trajectory)
    return trajectory


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


def unwrap_trajectory(trajectory: md.Trajectory) -> md.Trajectory:
    """Re-image molecules to undo periodic boundary wrapping.

    For trajectories with periodic box information, each molecule is
    translated by integer multiples of the box vectors so that it is
    contiguous (no bond crosses a box face).  Trajectories without
    unitcell information are returned unchanged.

    Args:
        trajectory: MDTraj trajectory. Shape: [N_frames, N_atoms, 3].

    Returns:
        md.Trajectory: Unwrapped trajectory. Shape: [N_frames, N_atoms, 3].
    """
    if trajectory.unitcell_lengths is None:
        return trajectory

    unwrapped = trajectory[:]
    molecules = unwrapped.topology.find_molecules()
    unwrapped.image_molecules(anchor_molecules=molecules, inplace=True)
    return unwrapped


def compute_water_diffusion_coefficient(
    positions_nm: np.ndarray,
    time_ps: np.ndarray,
    water_oxygen_indices: list[int],
    start_fraction: float = 0.2,
    end_fraction: float = 0.8,
) -> float:
    """Compute the self-diffusion coefficient of water from oxygen MSD.

    Uses the Einstein relation: D = (1/6) * d(MSD)/dt, where MSD is the
    mean-square displacement of water oxygen atoms averaged over all molecules.
    The slope is extracted from a linear fit to the MSD in the region
    [start_fraction, end_fraction] of the trajectory to avoid ballistic
    and poor-statistics regimes.

    Args:
        positions_nm: Trajectory positions in nm. Shape: [N_frames, N_atoms, 3].
        time_ps: Frame timestamps in ps. Shape: [N_frames].
        water_oxygen_indices: Atom indices of water oxygen atoms.
        start_fraction: Fraction of trajectory at which to begin the linear fit.
        end_fraction: Fraction of trajectory at which to end the linear fit.

    Returns:
        Self-diffusion coefficient in units of 10^{-5} cm^2/s.

    Raises:
        ValueError: If inputs are insufficient or invalid.
    """
    positions_nm = np.asarray(positions_nm)
    time_ps = np.asarray(time_ps)

    if positions_nm.ndim != 3 or positions_nm.shape[2] != 3:
        raise ValueError("positions_nm must have shape [N_frames, N_atoms, 3]")
    if time_ps.ndim != 1:
        raise ValueError("time_ps must be a 1-D array")
    if positions_nm.shape[0] != time_ps.shape[0]:
        raise ValueError("positions_nm and time_ps must have the same number of frames")
    if positions_nm.shape[0] < 100:
        raise ValueError("At least 100 frames are required for reliable MSD fitting")
    if len(water_oxygen_indices) < 10:
        raise ValueError("At least 10 water oxygen atoms are required")
    if not (0.0 <= start_fraction < end_fraction <= 1.0):
        raise ValueError("start_fraction must be less than end_fraction, both in [0, 1]")

    # Extract water oxygen positions: shape [N_frames, N_water, 3]
    oxy_pos = positions_nm[:, water_oxygen_indices, :]

    # Displacement from initial frame: r_i(t) - r_i(0)
    displacement = oxy_pos - oxy_pos[0:1, :, :]  # broadcast over frames

    # MSD(t) = <|r_i(t) - r_i(0)|^2> averaged over water molecules
    msd = np.mean(np.sum(displacement ** 2, axis=2), axis=1)  # shape [N_frames]

    # Select linear-fit region
    n_frames = len(time_ps)
    i_start = int(start_fraction * n_frames)
    i_end = int(end_fraction * n_frames)
    t_fit = time_ps[i_start:i_end]
    msd_fit = msd[i_start:i_end]

    # Linear fit: MSD = slope * t + intercept
    slope, _ = np.polyfit(t_fit, msd_fit, 1)

    # D = slope / 6, in nm^2/ps
    # Convert to 10^{-5} cm^2/s: 1 nm^2/ps = 1e-18 m^2 / 1e-12 s = 1e-6 m^2/s = 1e-1 * 10^{-5} cm^2/s
    diffusion_coefficient = (slope / 6.0) * 1e-1

    return diffusion_coefficient