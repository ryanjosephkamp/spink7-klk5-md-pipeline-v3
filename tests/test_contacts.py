"""Tests for interface-contact and hydrogen-bond analysis."""

from __future__ import annotations

import mdtraj as md
import numpy as np
import pytest
from openmm.app import Element, Topology

from src.analyze.contacts import compute_hbonds, compute_interface_contacts


def _make_contact_test_trajectory() -> tuple[md.Trajectory, np.ndarray, np.ndarray]:
    """Build a two-chain trajectory with controlled interface contacts."""

    topology = Topology()
    chain_a = topology.addChain()
    residue_a = topology.addResidue("ALA", chain_a)
    atom_a0 = topology.addAtom("CA", Element.getByAtomicNumber(6), residue_a)
    atom_a1 = topology.addAtom("CB", Element.getByAtomicNumber(6), residue_a)
    topology.addBond(atom_a0, atom_a1)

    chain_b = topology.addChain()
    residue_b = topology.addResidue("GLY", chain_b)
    atom_b0 = topology.addAtom("O", Element.getByAtomicNumber(8), residue_b)
    atom_b1 = topology.addAtom("N", Element.getByAtomicNumber(7), residue_b)
    topology.addBond(atom_b0, atom_b1)

    xyz = np.asarray(
        [
            [[0.00, 0.00, 0.00], [0.20, 0.00, 0.00], [0.30, 0.00, 0.00], [0.95, 0.00, 0.00]],
            [[0.00, 0.00, 0.00], [0.20, 0.00, 0.00], [0.70, 0.00, 0.00], [0.95, 0.00, 0.00]],
            [[0.00, 0.00, 0.00], [0.20, 0.00, 0.00], [0.08, 0.00, 0.00], [0.28, 0.00, 0.00]],
        ],
        dtype=np.float32,
    )

    trajectory = md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))
    chain_a_indices = np.asarray([0, 1], dtype=int)
    chain_b_indices = np.asarray([2, 3], dtype=int)
    return trajectory, chain_a_indices, chain_b_indices


def _make_hbond_test_trajectory() -> tuple[md.Trajectory, np.ndarray, np.ndarray]:
    """Build a minimal two-chain trajectory with one cross-interface hydrogen bond."""

    topology = Topology()
    chain_a = topology.addChain()
    residue_a = topology.addResidue("SER", chain_a)
    atom_n = topology.addAtom("N", Element.getByAtomicNumber(7), residue_a)
    atom_h = topology.addAtom("H", Element.getByAtomicNumber(1), residue_a)
    topology.addBond(atom_n, atom_h)

    chain_b = topology.addChain()
    residue_b = topology.addResidue("ASP", chain_b)
    atom_o = topology.addAtom("O", Element.getByAtomicNumber(8), residue_b)

    xyz = np.asarray(
        [
            [[0.00, 0.00, 0.00], [0.10, 0.00, 0.00], [0.22, 0.00, 0.00]],
            [[0.00, 0.00, 0.00], [0.10, 0.00, 0.00], [0.23, 0.00, 0.00]],
            [[0.00, 0.00, 0.00], [0.10, 0.00, 0.00], [0.45, 0.00, 0.00]],
            [[0.00, 0.00, 0.00], [0.10, 0.00, 0.00], [0.50, 0.00, 0.00]],
        ],
        dtype=np.float32,
    )

    trajectory = md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))
    chain_a_indices = np.asarray([0, 1], dtype=int)
    chain_b_indices = np.asarray([2], dtype=int)
    return trajectory, chain_a_indices, chain_b_indices


def test_compute_interface_contacts_returns_schema_and_bounded_frequencies() -> None:
    """Chunk 17 gate: contact analysis must match the ratified data contract."""

    trajectory, chain_a_indices, chain_b_indices = _make_contact_test_trajectory()

    results = compute_interface_contacts(trajectory, chain_a_indices, chain_b_indices, cutoff_nm=0.15)

    assert results["contact_map"].shape == (trajectory.n_frames, chain_a_indices.size, chain_b_indices.size)
    assert results["contact_map"].dtype == np.bool_
    assert results["n_contacts_per_frame"].shape == (trajectory.n_frames,)
    assert results["contact_frequency"].shape == (1, 1)
    assert np.all(results["contact_frequency"] >= 0.0)
    assert np.all(results["contact_frequency"] <= 1.0)
    assert np.array_equal(results["n_contacts_per_frame"], np.asarray([1, 0, 3], dtype=int))
    assert results["contact_frequency"][0, 0] == pytest.approx(2.0 / 3.0)


def test_compute_interface_contacts_rejects_overlapping_groups() -> None:
    """Public API should reject overlapping interface groups."""

    trajectory, _, chain_b_indices = _make_contact_test_trajectory()

    with pytest.raises(ValueError, match="must be disjoint"):
        compute_interface_contacts(trajectory, np.asarray([0, 1], dtype=int), chain_b_indices=np.asarray([1, 2], dtype=int))


def test_compute_hbonds_identifies_cross_interface_donor_acceptor_triplets() -> None:
    """Hydrogen-bond analysis should report donor-hydrogen-acceptor triplets across the interface."""

    trajectory, chain_a_indices, chain_b_indices = _make_hbond_test_trajectory()

    results = compute_hbonds(trajectory, chain_a_indices, chain_b_indices)

    assert results["hbond_triplets"].shape == (1, 3)
    assert results["hbond_frequency"].shape == (1,)
    assert results["hbond_triplets"][0, 0] in chain_a_indices
    assert results["hbond_triplets"][0, 1] in chain_a_indices
    assert results["hbond_triplets"][0, 2] in chain_b_indices
    assert np.all(results["hbond_frequency"] >= 0.0)
    assert np.all(results["hbond_frequency"] <= 1.0)
    assert results["hbond_frequency"][0] == pytest.approx(0.5, abs=1e-8)


def test_compute_hbonds_returns_empty_arrays_without_cross_interface_candidates() -> None:
    """Hydrogen-bond analysis should return empty arrays when no valid donor-acceptor pairs exist."""

    trajectory, chain_a_indices, chain_b_indices = _make_contact_test_trajectory()

    results = compute_hbonds(trajectory, chain_a_indices, chain_b_indices)

    assert results["hbond_triplets"].shape == (0, 3)
    assert results["hbond_frequency"].shape == (0,)