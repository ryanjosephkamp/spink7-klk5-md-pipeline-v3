"""Interface contact and hydrogen-bond analysis utilities."""

from __future__ import annotations

import mdtraj as md
import numpy as np


def _validate_atom_index_array(name: str, atom_indices: np.ndarray, n_atoms: int) -> np.ndarray:
    """Validate a one-dimensional atom-index array at the public API boundary."""

    validated = np.asarray(atom_indices, dtype=int)
    assert validated.ndim == 1, f"{name} must have shape [N]"
    if validated.ndim != 1 or validated.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional array")
    if np.any(validated < 0) or np.any(validated >= n_atoms):
        raise ValueError(f"{name} must contain valid atom indices")
    if np.unique(validated).size != validated.size:
        raise ValueError(f"{name} must not contain duplicate atom indices")
    return validated


def _validate_disjoint_groups(chain_a_indices: np.ndarray, chain_b_indices: np.ndarray) -> None:
    """Reject overlapping interface groups."""

    if np.intersect1d(chain_a_indices, chain_b_indices).size != 0:
        raise ValueError("chain_a_indices and chain_b_indices must be disjoint")


def _residue_inverse_indices(topology: md.Topology, atom_indices: np.ndarray) -> tuple[np.ndarray, int]:
    """Map atom indices to a compact residue-index domain for aggregation."""

    residue_ids = np.asarray([topology.atom(int(atom_index)).residue.index for atom_index in atom_indices], dtype=int)
    _, inverse = np.unique(residue_ids, return_inverse=True)
    return inverse.astype(int, copy=False), int(inverse.max(initial=-1) + 1)


def _collect_donor_hydrogen_pairs(topology: md.Topology, atom_mask: np.ndarray) -> list[tuple[int, int]]:
    """Identify donor-heavy-atom / hydrogen pairs inside a selected group."""

    donor_pairs: list[tuple[int, int]] = []
    donor_elements = {"N", "O", "S"}
    for atom_a, atom_b in topology.bonds:
        atom_a_symbol = atom_a.element.symbol if atom_a.element is not None else ""
        atom_b_symbol = atom_b.element.symbol if atom_b.element is not None else ""
        if atom_a_symbol == "H" and atom_b_symbol in donor_elements and atom_mask[atom_a.index] and atom_mask[atom_b.index]:
            donor_pairs.append((atom_b.index, atom_a.index))
        elif atom_b_symbol == "H" and atom_a_symbol in donor_elements and atom_mask[atom_a.index] and atom_mask[atom_b.index]:
            donor_pairs.append((atom_a.index, atom_b.index))
    return donor_pairs


def _collect_acceptors(topology: md.Topology, atom_mask: np.ndarray, donor_heavy_atoms: set[int]) -> np.ndarray:
    """Identify acceptor heavy atoms inside a selected group."""

    acceptor_elements = {"N", "O", "S"}
    acceptors = [
        atom.index
        for atom in topology.atoms
        if atom_mask[atom.index]
        and atom.element is not None
        and atom.element.symbol in acceptor_elements
        and atom.index not in donor_heavy_atoms
    ]
    return np.asarray(acceptors, dtype=int)


def compute_interface_contacts(
    trajectory: md.Trajectory,
    chain_a_indices: np.ndarray,
    chain_b_indices: np.ndarray,
    cutoff_nm: float = 0.45,
) -> dict[str, np.ndarray]:
    """Compute atom-level interface contacts and residue-pair contact frequencies.

    Args:
        trajectory: MD trajectory. Shape: [N_frames, N_atoms, 3].
        chain_a_indices: Atom indices for interface group A. Shape: [N_a].
        chain_b_indices: Atom indices for interface group B. Shape: [N_b].
        cutoff_nm: Distance cutoff for a contact in nm.

    Returns:
        dict[str, np.ndarray]:
            "contact_map": bool array. Shape: [N_frames, N_a, N_b].
            "n_contacts_per_frame": int array. Shape: [N_frames].
            "contact_frequency": float array. Shape: [N_residues_a, N_residues_b].
    """

    assert trajectory.xyz.ndim == 3, "trajectory coordinates must have shape [N_frames, N_atoms, 3]"
    assert trajectory.xyz.shape[2] == 3, "trajectory coordinates must have shape [N_frames, N_atoms, 3]"

    group_a = _validate_atom_index_array("chain_a_indices", chain_a_indices, trajectory.n_atoms)
    group_b = _validate_atom_index_array("chain_b_indices", chain_b_indices, trajectory.n_atoms)
    _validate_disjoint_groups(group_a, group_b)

    if cutoff_nm <= 0.0:
        raise ValueError("cutoff_nm must be positive")

    # Pairwise distance matrix via broadcasting: [frames, N_a, N_b, 3].
    positions_a = trajectory.xyz[:, group_a, :]
    positions_b = trajectory.xyz[:, group_b, :]
    deltas = positions_a[:, :, np.newaxis, :] - positions_b[:, np.newaxis, :, :]
    squared_distances = np.sum(deltas * deltas, axis=3, dtype=float)
    contact_map = squared_distances <= float(cutoff_nm * cutoff_nm)
    n_contacts_per_frame = np.count_nonzero(contact_map, axis=(1, 2)).astype(int, copy=False)

    # Map atom-level contacts to residue-pair domain for frequency analysis.
    residue_a_inverse, n_residues_a = _residue_inverse_indices(trajectory.topology, group_a)
    residue_b_inverse, n_residues_b = _residue_inverse_indices(trajectory.topology, group_b)
    # Flatten residue pair (i, j) to a single index: i * n_residues_b + j.
    residue_pair_ids = (residue_a_inverse[:, np.newaxis] * n_residues_b + residue_b_inverse[np.newaxis, :]).reshape(-1)

    residue_contact_map = np.zeros((trajectory.n_frames, n_residues_a * n_residues_b), dtype=bool)
    contact_frames, contact_pair_indices = np.nonzero(contact_map.reshape(trajectory.n_frames, -1))
    if contact_frames.size != 0:
        residue_contact_map[contact_frames, residue_pair_ids[contact_pair_indices]] = True
    contact_frequency = np.mean(residue_contact_map, axis=0, dtype=float).reshape(n_residues_a, n_residues_b)

    return {
        "contact_map": contact_map,
        "n_contacts_per_frame": n_contacts_per_frame,
        "contact_frequency": np.asarray(contact_frequency, dtype=float),
    }


def compute_hbonds(
    trajectory: md.Trajectory,
    chain_a_indices: np.ndarray,
    chain_b_indices: np.ndarray,
) -> dict[str, np.ndarray]:
    """Compute cross-interface hydrogen-bond triplets and their per-trajectory frequencies.

    A hydrogen bond is counted when the H···A distance is at most 0.25 nm and the
    D-H···A angle is at least 120 degrees.

    Args:
        trajectory: MD trajectory. Shape: [N_frames, N_atoms, 3].
        chain_a_indices: Atom indices for interface group A. Shape: [N_a].
        chain_b_indices: Atom indices for interface group B. Shape: [N_b].

    Returns:
        dict[str, np.ndarray]:
            "hbond_triplets": int array of donor, hydrogen, acceptor indices. Shape: [N_hbonds, 3].
            "hbond_frequency": float array in [0, 1]. Shape: [N_hbonds].
    """

    assert trajectory.xyz.ndim == 3, "trajectory coordinates must have shape [N_frames, N_atoms, 3]"
    assert trajectory.xyz.shape[2] == 3, "trajectory coordinates must have shape [N_frames, N_atoms, 3]"

    group_a = _validate_atom_index_array("chain_a_indices", chain_a_indices, trajectory.n_atoms)
    group_b = _validate_atom_index_array("chain_b_indices", chain_b_indices, trajectory.n_atoms)
    _validate_disjoint_groups(group_a, group_b)

    atom_mask_a = np.zeros(trajectory.n_atoms, dtype=bool)
    atom_mask_b = np.zeros(trajectory.n_atoms, dtype=bool)
    atom_mask_a[group_a] = True
    atom_mask_b[group_b] = True

    donor_pairs_a = _collect_donor_hydrogen_pairs(trajectory.topology, atom_mask_a)
    donor_pairs_b = _collect_donor_hydrogen_pairs(trajectory.topology, atom_mask_b)
    donor_heavy_atoms = {donor for donor, _ in donor_pairs_a + donor_pairs_b}
    acceptors_a = _collect_acceptors(trajectory.topology, atom_mask_a, donor_heavy_atoms)
    acceptors_b = _collect_acceptors(trajectory.topology, atom_mask_b, donor_heavy_atoms)

    # Build cross-interface D-H···A candidate triplets: donors from one
    # chain paired with acceptors from the other chain.
    candidate_triplets: list[tuple[int, int, int]] = []
    for donor_index, hydrogen_index in donor_pairs_a:
        for acceptor_index in acceptors_b:
            if donor_index != acceptor_index:
                candidate_triplets.append((donor_index, hydrogen_index, int(acceptor_index)))
    for donor_index, hydrogen_index in donor_pairs_b:
        for acceptor_index in acceptors_a:
            if donor_index != acceptor_index:
                candidate_triplets.append((donor_index, hydrogen_index, int(acceptor_index)))

    if not candidate_triplets:
        return {
            "hbond_triplets": np.zeros((0, 3), dtype=int),
            "hbond_frequency": np.zeros((0,), dtype=float),
        }

    triplets = np.asarray(candidate_triplets, dtype=int)
    donor_xyz = trajectory.xyz[:, triplets[:, 0], :]
    hydrogen_xyz = trajectory.xyz[:, triplets[:, 1], :]
    acceptor_xyz = trajectory.xyz[:, triplets[:, 2], :]

    # H-bond geometry: compute H···A distance and D-H···A angle.
    hydrogen_acceptor = acceptor_xyz - hydrogen_xyz
    donor_hydrogen = donor_xyz - hydrogen_xyz
    hydrogen_acceptor_distances = np.linalg.norm(hydrogen_acceptor, axis=2)

    # cos(theta) = (D-H) . (H-A) / (|D-H| * |H-A|); clip for arccos stability.
    donor_hydrogen_norm = np.linalg.norm(donor_hydrogen, axis=2)
    denominator = np.clip(donor_hydrogen_norm * hydrogen_acceptor_distances, a_min=np.finfo(float).eps, a_max=None)
    cos_theta = np.sum(donor_hydrogen * hydrogen_acceptor, axis=2, dtype=float) / denominator
    angles_deg = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

    # Apply geometric criteria: H···A <= 0.25 nm and D-H···A >= 120 degrees.
    present = (hydrogen_acceptor_distances <= 0.25) & (angles_deg >= 120.0)
    frequencies = np.mean(present, axis=0, dtype=float)
    keep_mask = frequencies > 0.0

    return {
        "hbond_triplets": triplets[keep_mask],
        "hbond_frequency": np.asarray(frequencies[keep_mask], dtype=float),
    }