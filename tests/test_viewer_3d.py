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


def test_apply_base_chain_style_does_not_raise_for_many_chains(tmp_path: Path) -> None:
    """Palette cycling should handle more chains than palette entries without error."""

    from src.visualization.viewer_3d import _apply_base_chain_style

    chain_ids = [chr(ord("A") + i) for i in range(20)]
    lines = []
    for idx, cid in enumerate(chain_ids):
        x = float(idx) * 0.5
        lines.append(
            f"ATOM  {idx + 1:5d}  CA  ALA {cid}{idx + 1:4d}    {x:8.3f}   0.000   0.000  1.00  0.00           C"
        )
    lines.append("END")
    pdb_text = "\n".join(lines)

    view = py3Dmol.view(width=400, height=300)
    view.addModel(pdb_text, "pdb")

    # Should not raise IndexError
    _apply_base_chain_style(view, "cartoon", pdb_text)


def test_render_complex_accepts_custom_chain_colors(tmp_path: Path) -> None:
    """render_complex should accept a custom chain_colors parameter."""

    topology = Topology()
    chain_a = topology.addChain("A")
    res_a = topology.addResidue("ALA", chain_a)
    topology.addAtom("CA", Element.getByAtomicNumber(6), res_a)
    chain_b = topology.addChain("B")
    res_b = topology.addResidue("GLY", chain_b)
    topology.addAtom("CA", Element.getByAtomicNumber(6), res_b)

    xyz = np.array([[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]], dtype=np.float32)
    traj = md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))
    pdb_path = tmp_path / "test.pdb"
    traj.save_pdb(str(pdb_path))

    view = render_complex(pdb_path, chain_colors=["gold", "silver"])
    assert isinstance(view, type(py3Dmol.view()))


def test_render_complex_handles_many_chains(tmp_path: Path) -> None:
    """Complexes with more chains than palette entries should render without error."""

    topology = Topology()
    atoms = []
    for i in range(15):
        chain = topology.addChain(chr(ord("A") + i))
        residue = topology.addResidue("ALA", chain)
        atom = topology.addAtom("CA", Element.getByAtomicNumber(6), residue)
        atoms.append(atom)

    xyz = np.zeros((1, 15, 3), dtype=np.float32)
    for i in range(15):
        xyz[0, i, 0] = float(i) * 0.3

    traj = md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))
    pdb_path = tmp_path / "many_chains.pdb"
    traj.save_pdb(str(pdb_path))

    view = render_complex(pdb_path, style="cartoon")
    assert isinstance(view, type(py3Dmol.view()))