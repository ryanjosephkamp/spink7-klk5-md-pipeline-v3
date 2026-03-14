"""Tests for trajectory analysis utilities."""

from __future__ import annotations

from pathlib import Path

import mdtraj as md
import numpy as np
import pytest
from openmm.app import Element, Topology

from src.analyze.trajectory import align_trajectory, iterload_trajectory, load_trajectory, select_frames


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