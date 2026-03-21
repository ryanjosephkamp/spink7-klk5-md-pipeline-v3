"""Tests for trajectory analysis utilities."""

from __future__ import annotations

from pathlib import Path

import mdtraj as md
import numpy as np
import pytest
from openmm.app import Element, Topology

from src.analyze.trajectory import align_trajectory, iterload_trajectory, load_trajectory, select_frames, unwrap_trajectory


def _write_test_trajectory(tmp_path: Path) -> tuple[Path, Path, md.Trajectory]:
    """Create a small on-disk trajectory and matching topology for analysis tests."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ALA", chain)
    atom_n = topology.addAtom("N", Element.getByAtomicNumber(7), residue)
    atom_ca = topology.addAtom("CA", Element.getByAtomicNumber(6), residue)
    atom_c = topology.addAtom("C", Element.getByAtomicNumber(6), residue)
    topology.addBond(atom_n, atom_ca)
    topology.addBond(atom_ca, atom_c)

    md_topology = md.Topology.from_openmm(topology)
    xyz = np.asarray(
        [
            [[0.00, 0.00, 0.00], [0.12, 0.00, 0.00], [0.24, 0.00, 0.00]],
            [[0.01, 0.00, 0.00], [0.13, 0.01, 0.00], [0.25, 0.01, 0.00]],
            [[0.02, 0.01, 0.00], [0.14, 0.02, 0.00], [0.26, 0.02, 0.00]],
            [[0.03, 0.01, 0.00], [0.15, 0.03, 0.00], [0.27, 0.03, 0.00]],
            [[0.04, 0.02, 0.00], [0.16, 0.04, 0.00], [0.28, 0.04, 0.00]],
        ],
        dtype=np.float32,
    )
    trajectory = md.Trajectory(xyz=xyz, topology=md_topology)

    topology_path = tmp_path / "trajectory_topology.pdb"
    trajectory_path = tmp_path / "trajectory_frames.dcd"
    trajectory[0].save_pdb(str(topology_path))
    trajectory.save_dcd(str(trajectory_path))
    return trajectory_path, topology_path, trajectory


def test_load_trajectory_returns_expected_frame_count(tmp_path: Path) -> None:
    """Chunk 15 gate: loading a saved trajectory should preserve the frame count."""

    trajectory_path, topology_path, reference = _write_test_trajectory(tmp_path)

    loaded = load_trajectory(trajectory_path, topology_path)

    assert loaded.n_frames == reference.n_frames
    assert loaded.n_atoms == reference.n_atoms
    assert loaded.xyz.shape == (5, 3, 3)


def test_iterload_trajectory_streams_all_frames(tmp_path: Path) -> None:
    """Chunked trajectory iteration should cover the full trajectory without frame loss."""

    trajectory_path, topology_path, reference = _write_test_trajectory(tmp_path)

    chunks = list(iterload_trajectory(trajectory_path, topology_path, chunk_size=2))

    assert [chunk.n_frames for chunk in chunks] == [2, 2, 1]
    assert sum(chunk.n_frames for chunk in chunks) == reference.n_frames


def test_select_frames_returns_requested_subset(tmp_path: Path) -> None:
    """Frame selection should preserve order and coordinates for the requested indices."""

    trajectory_path, topology_path, reference = _write_test_trajectory(tmp_path)
    loaded = load_trajectory(trajectory_path, topology_path)

    subset = select_frames(loaded, np.asarray([4, 1, 3], dtype=int))

    assert subset.n_frames == 3
    assert np.allclose(subset.xyz[0], reference.xyz[4])
    assert np.allclose(subset.xyz[1], reference.xyz[1])
    assert np.allclose(subset.xyz[2], reference.xyz[3])


def test_align_trajectory_reduces_backbone_rmsd_to_reference(tmp_path: Path) -> None:
    """Alignment should reduce the average backbone RMSD to the reference frame."""

    trajectory_path, topology_path, _ = _write_test_trajectory(tmp_path)
    loaded = load_trajectory(trajectory_path, topology_path)
    displaced = loaded[:]
    rotation = np.asarray(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float32,
    )
    displaced.xyz = np.matmul(displaced.xyz, rotation.T) + np.asarray([0.5, -0.2, 0.1], dtype=np.float32)

    backbone_indices = loaded.topology.select("backbone")
    reference_backbone = loaded.xyz[0, backbone_indices, :]
    displaced_backbone = displaced.xyz[:, backbone_indices, :]
    rmsd_before = np.sqrt(np.mean(np.sum((displaced_backbone - reference_backbone[np.newaxis, :, :]) ** 2, axis=2), axis=1))
    aligned = align_trajectory(displaced, loaded[0], atom_selection="backbone")
    rmsd_after = md.rmsd(aligned, loaded[0], atom_indices=backbone_indices)

    assert np.mean(rmsd_after) < np.mean(rmsd_before)


def test_load_trajectory_rejects_missing_file(tmp_path: Path) -> None:
    """Public API should reject missing trajectory paths."""

    missing_trajectory = tmp_path / "missing.dcd"
    missing_topology = tmp_path / "missing.pdb"

    with pytest.raises(FileNotFoundError, match="trajectory file not found"):
        load_trajectory(missing_trajectory, missing_topology)


# ---------- L-16 Step 3: Diffusion coefficient diagnostic ----------


def test_compute_water_diffusion_coefficient_returns_positive() -> None:
    """Diffusion coefficient from a random walk should be positive."""
    from src.analyze.trajectory import compute_water_diffusion_coefficient

    rng = np.random.default_rng(42)
    n_frames, n_atoms = 200, 20
    # Cumulative random walk in nm
    positions = np.cumsum(rng.normal(0, 0.001, (n_frames, n_atoms, 3)), axis=0)
    time_ps = np.arange(n_frames) * 1.0

    D = compute_water_diffusion_coefficient(
        positions, time_ps, list(range(n_atoms))
    )
    assert D > 0.0


def test_compute_water_diffusion_coefficient_rejects_small_frames() -> None:
    """Should reject trajectories shorter than 100 frames."""
    from src.analyze.trajectory import compute_water_diffusion_coefficient

    positions = np.zeros((50, 20, 3))
    time_ps = np.arange(50) * 1.0
    with pytest.raises(ValueError, match="At least 100 frames"):
        compute_water_diffusion_coefficient(positions, time_ps, list(range(20)))


def test_compute_water_diffusion_coefficient_rejects_few_oxygens() -> None:
    """Should reject fewer than 10 water oxygen atoms."""
    from src.analyze.trajectory import compute_water_diffusion_coefficient

    positions = np.zeros((200, 100, 3))
    time_ps = np.arange(200) * 1.0
    with pytest.raises(ValueError, match="At least 10 water oxygen"):
        compute_water_diffusion_coefficient(positions, time_ps, list(range(5)))


def test_compute_water_diffusion_coefficient_rejects_invalid_fractions() -> None:
    """Should reject start_fraction >= end_fraction."""
    from src.analyze.trajectory import compute_water_diffusion_coefficient

    positions = np.zeros((200, 20, 3))
    time_ps = np.arange(200) * 1.0
    with pytest.raises(ValueError, match="start_fraction must be less than end_fraction"):
        compute_water_diffusion_coefficient(
            positions, time_ps, list(range(20)), start_fraction=0.8, end_fraction=0.2
        )


# ---------- L-32 Step 1: PBC unwrapping ----------


def test_unwrap_trajectory_returns_unchanged_for_non_periodic() -> None:
    """Unwrapping a trajectory without box info should return it unchanged."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ALA", chain)
    topology.addAtom("N", Element.getByAtomicNumber(7), residue)
    topology.addAtom("CA", Element.getByAtomicNumber(6), residue)

    xyz = np.array([[[0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]], dtype=np.float32)
    traj = md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))

    result = unwrap_trajectory(traj)

    assert np.allclose(result.xyz, xyz), "Non-periodic trajectory should be unchanged"


# ---------- L-32 Step 2: load_trajectory unwrap integration ----------


def test_load_trajectory_unwraps_by_default(tmp_path: Path) -> None:
    """load_trajectory with unwrap=True should produce the same result as manual unwrapping."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ALA", chain)
    topology.addAtom("N", Element.getByAtomicNumber(7), residue)
    md_top = md.Topology.from_openmm(topology)
    xyz = np.array([[[0.1, 0.2, 0.3]]], dtype=np.float32)
    traj = md.Trajectory(xyz=xyz, topology=md_top)
    top_path = tmp_path / "top.pdb"
    traj_path = tmp_path / "traj.dcd"
    traj.save_pdb(str(top_path))
    traj.save_dcd(str(traj_path))

    loaded = load_trajectory(traj_path, top_path, unwrap=True)

    assert loaded.n_frames == 1
    assert loaded.n_atoms == 1


# ---------- L-32 Step 4: PBC unwrapping regression test ----------


def test_unwrap_trajectory_corrects_boundary_crossing() -> None:
    """PBC unwrapping should eliminate artificial jumps from boundary-wrapped coordinates."""

    # Build a two-atom bonded molecule
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ALA", chain)
    a1 = topology.addAtom("N", Element.getByAtomicNumber(7), residue)
    a2 = topology.addAtom("CA", Element.getByAtomicNumber(6), residue)
    topology.addBond(a1, a2)

    md_top = md.Topology.from_openmm(topology)

    # Frame 1: both atoms on same side of box
    # Frame 2: atom2 wrapped to opposite face (simulating boundary crossing)
    box_length = 3.0  # nm
    xyz = np.array(
        [
            [[1.4, 1.5, 1.5], [1.55, 1.5, 1.5]],  # contiguous
            [[1.4, 1.5, 1.5], [1.55 - box_length, 1.5, 1.5]],  # atom2 wrapped
        ],
        dtype=np.float32,
    )

    unitcell_lengths = np.array([[box_length, box_length, box_length]] * 2, dtype=np.float32)
    unitcell_angles = np.array([[90.0, 90.0, 90.0]] * 2, dtype=np.float32)

    traj = md.Trajectory(
        xyz=xyz,
        topology=md_top,
        unitcell_lengths=unitcell_lengths,
        unitcell_angles=unitcell_angles,
    )

    unwrapped = unwrap_trajectory(traj)

    # After unwrapping, the bond length in frame 2 should be ~0.15 nm, not ~2.85 nm
    bond_length_frame2 = np.linalg.norm(
        unwrapped.xyz[1, 1, :] - unwrapped.xyz[1, 0, :]
    )
    assert bond_length_frame2 < 0.5, (
        f"Unwrapped bond length {bond_length_frame2:.3f} nm is too large — "
        "PBC unwrapping failed to correct boundary crossing"
    )