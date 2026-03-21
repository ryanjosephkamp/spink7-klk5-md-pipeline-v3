"""Tests for structural observables."""

from __future__ import annotations

from pathlib import Path

import mdtraj as md
import numpy as np
import pytest
from openmm.app import Element, Topology

from src.analyze.structural import compute_radius_of_gyration, compute_rmsd, compute_rmsf, compute_sasa


def _make_structural_test_trajectory() -> md.Trajectory:
    """Create a compact peptide-like trajectory for structural analysis tests."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ALA", chain)
    atom_n = topology.addAtom("N", Element.getByAtomicNumber(7), residue)
    atom_ca = topology.addAtom("CA", Element.getByAtomicNumber(6), residue)
    atom_c = topology.addAtom("C", Element.getByAtomicNumber(6), residue)
    atom_o = topology.addAtom("O", Element.getByAtomicNumber(8), residue)
    topology.addBond(atom_n, atom_ca)
    topology.addBond(atom_ca, atom_c)
    topology.addBond(atom_c, atom_o)

    xyz = np.asarray(
        [
            [[0.00, 0.00, 0.00], [0.12, 0.00, 0.00], [0.24, 0.00, 0.00], [0.32, 0.03, 0.00]],
            [[0.00, 0.01, 0.00], [0.12, 0.01, 0.01], [0.24, 0.02, 0.00], [0.32, 0.04, 0.01]],
            [[0.01, 0.01, 0.00], [0.13, 0.02, 0.01], [0.25, 0.03, 0.01], [0.33, 0.05, 0.01]],
            [[0.01, 0.02, 0.01], [0.13, 0.03, 0.02], [0.25, 0.04, 0.01], [0.34, 0.06, 0.02]],
        ],
        dtype=np.float32,
    )
    return md.Trajectory(xyz=xyz, topology=md.Topology.from_openmm(topology))


def test_compute_rmsd_returns_zero_for_reference_against_itself() -> None:
    """Chunk 16 gate: RMSD should be zero when a trajectory is compared to itself."""

    trajectory = _make_structural_test_trajectory()

    rmsd_nm = compute_rmsd(trajectory, trajectory[0], atom_selection="backbone")

    assert rmsd_nm.shape == (trajectory.n_frames,)
    assert rmsd_nm[0] == pytest.approx(0.0, abs=1e-8)
    assert np.all(rmsd_nm >= 0.0)


def test_compute_rmsf_returns_selected_atom_count() -> None:
    """RMSF should return one value per selected atom."""

    trajectory = _make_structural_test_trajectory()

    rmsf_nm = compute_rmsf(trajectory, atom_selection="name CA")

    assert rmsf_nm.shape == (1,)
    assert rmsf_nm[0] >= 0.0


def test_compute_radius_of_gyration_is_positive_per_frame() -> None:
    """Radius of gyration should be positive for every frame."""

    trajectory = _make_structural_test_trajectory()

    rg_nm = compute_radius_of_gyration(trajectory)

    assert rg_nm.shape == (trajectory.n_frames,)
    assert np.all(rg_nm > 0.0)


def test_compute_sasa_returns_positive_framewise_values() -> None:
    """SASA should return strictly positive total surface area per frame."""

    trajectory = _make_structural_test_trajectory()

    sasa_nm2 = compute_sasa(trajectory, probe_radius_nm=0.14)

    assert sasa_nm2.shape == (trajectory.n_frames,)
    assert np.all(sasa_nm2 > 0.0)


def test_compute_sasa_rejects_non_positive_probe_radius() -> None:
    """Public API should reject non-physical SASA probe radii."""

    trajectory = _make_structural_test_trajectory()

    with pytest.raises(ValueError, match="probe_radius_nm must be positive"):
        compute_sasa(trajectory, probe_radius_nm=0.0)


# ---------- L-32 Step 3: Structural functions accept unwrap parameter ----------


def test_compute_rmsd_accepts_unwrap_parameter() -> None:
    """compute_rmsd should accept and honor the unwrap parameter."""

    trajectory = _make_structural_test_trajectory()

    rmsd_unwrap = compute_rmsd(trajectory, trajectory[0], unwrap=True)
    rmsd_no_unwrap = compute_rmsd(trajectory, trajectory[0], unwrap=False)

    assert rmsd_unwrap.shape == (trajectory.n_frames,)
    assert rmsd_no_unwrap.shape == (trajectory.n_frames,)
    assert np.allclose(rmsd_unwrap, rmsd_no_unwrap, atol=1e-6)