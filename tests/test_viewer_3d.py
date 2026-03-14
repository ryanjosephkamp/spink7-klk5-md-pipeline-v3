"""Tests for py3Dmol visualization helpers."""

from __future__ import annotations

from pathlib import Path

import mdtraj as md
import numpy as np
import py3Dmol
from openmm.app import Element, Topology

from src.visualization.viewer_3d import render_complex, render_trajectory_frame


def _write_test_complex_pdb(tmp_path: Path) -> Path:
    """Write a two-chain PDB fixture with HIS/ASP/SER residues for styling tests."""

    topology = Topology()

    chain_a = topology.addChain("A")
    residue_his = topology.addResidue("HIS", chain_a, id="1")
    residue_ser = topology.addResidue("SER", chain_a, id="2")
    atom_a_n = topology.addAtom("N", Element.getByAtomicNumber(7), residue_his)
    atom_a_ca = topology.addAtom("CA", Element.getByAtomicNumber(6), residue_his)
    atom_a_cb = topology.addAtom("CB", Element.getByAtomicNumber(6), residue_ser)
    topology.addBond(atom_a_n, atom_a_ca)
    topology.addBond(atom_a_ca, atom_a_cb)

    chain_b = topology.addChain("B")
    residue_asp = topology.addResidue("ASP", chain_b, id="3")
    residue_gly = topology.addResidue("GLY", chain_b, id="4")
    atom_b_n = topology.addAtom("N", Element.getByAtomicNumber(7), residue_asp)
    atom_b_ca = topology.addAtom("CA", Element.getByAtomicNumber(6), residue_asp)
    atom_b_c = topology.addAtom("C", Element.getByAtomicNumber(6), residue_gly)
    topology.addBond(atom_b_n, atom_b_ca)
    topology.addBond(atom_b_ca, atom_b_c)

    xyz = np.asarray(
        [[
            [0.00, 0.00, 0.00],
            [0.15, 0.00, 0.00],
            [0.30, 0.05, 0.00],
            [0.55, 0.00, 0.00],
            [0.70, 0.00, 0.00],
            [0.85, 0.05, 0.00],
        ]],
        dtype=np.float32,
    )
    trajectory = md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))
    pdb_path = tmp_path / "complex_fixture.pdb"
    trajectory.save_pdb(str(pdb_path))
    return pdb_path


def _make_test_trajectory() -> md.Trajectory:
    """Build a small multi-frame trajectory for visualization tests."""

    topology = Topology()
    chain_a = topology.addChain("A")
    residue_a = topology.addResidue("ALA", chain_a, id="1")
    atom_a_n = topology.addAtom("N", Element.getByAtomicNumber(7), residue_a)
    atom_a_ca = topology.addAtom("CA", Element.getByAtomicNumber(6), residue_a)
    topology.addBond(atom_a_n, atom_a_ca)

    chain_b = topology.addChain("B")
    residue_b = topology.addResidue("SER", chain_b, id="2")
    atom_b_n = topology.addAtom("N", Element.getByAtomicNumber(7), residue_b)
    atom_b_ca = topology.addAtom("CA", Element.getByAtomicNumber(6), residue_b)
    topology.addBond(atom_b_n, atom_b_ca)

    xyz = np.asarray(
        [
            [[0.00, 0.00, 0.00], [0.14, 0.00, 0.00], [0.50, 0.00, 0.00], [0.64, 0.00, 0.00]],
            [[0.02, 0.00, 0.00], [0.16, 0.01, 0.00], [0.52, 0.00, 0.00], [0.66, 0.01, 0.00]],
        ],
        dtype=np.float32,
    )
    return md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))


def test_render_complex_returns_py3dmol_view(tmp_path: Path) -> None:
    """Chunk 21 gate: complex rendering should return a valid py3Dmol view object."""

    pdb_path = _write_test_complex_pdb(tmp_path)

    view = render_complex(pdb_path, highlight_interface_residues=[1, 3], style="cartoon")

    assert isinstance(view, type(py3Dmol.view()))


def test_render_trajectory_frame_returns_py3dmol_view_for_valid_frame() -> None:
    """Trajectory-frame rendering should return a valid py3Dmol view object for a legal frame index."""

    trajectory = _make_test_trajectory()

    view = render_trajectory_frame(trajectory, frame_index=1, color_by="chain")

    assert isinstance(view, type(py3Dmol.view()))


def test_render_trajectory_frame_rejects_invalid_frame_index() -> None:
    """Public API should reject out-of-range trajectory frame indices."""

    trajectory = _make_test_trajectory()

    try:
        render_trajectory_frame(trajectory, frame_index=2, color_by="chain")
    except ValueError as exc:
        assert "frame_index" in str(exc)
    else:
        raise AssertionError("render_trajectory_frame should reject out-of-range frame indices")