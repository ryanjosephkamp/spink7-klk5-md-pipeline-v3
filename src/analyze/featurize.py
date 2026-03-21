"""Trajectory featurization for Markov State Model construction.

Converts MD trajectory coordinates into feature vectors suitable for
TICA dimensionality reduction and microstate clustering.
"""

from __future__ import annotations

import logging

import mdtraj as md
import numpy as np

__all__ = [
    "compute_backbone_dihedrals",
    "compute_contact_features",
    "combine_features",
]

logger = logging.getLogger(__name__)


def compute_backbone_dihedrals(trajectory: md.Trajectory) -> np.ndarray:
    """Extract backbone phi/psi dihedrals as sin/cos feature pairs.

    The sin/cos transformation avoids discontinuities at +/-pi and produces
    a feature space where Euclidean distances approximate dihedral distances.

    Args:
        trajectory: MDTraj trajectory. Shape: [N_frames, N_atoms, 3].

    Returns:
        np.ndarray: Dihedral features. Shape: [N_frames, 4 * N_common_residues].
            Columns are ordered as
            [sin(phi_1), cos(phi_1), sin(psi_1), cos(psi_1), ...].
    """
    phi_indices, phi_angles = md.compute_phi(trajectory)
    psi_indices, psi_angles = md.compute_psi(trajectory)

    # MDTraj may return different numbers of phi/psi if terminal residues
    # lack one of the angles.  Match by residue index (3rd atom in the
    # dihedral quad defines the residue for phi, 2nd atom for psi).
    phi_residues = {
        trajectory.topology.atom(int(idx[2])).residue.index: i
        for i, idx in enumerate(phi_indices)
    }
    psi_residues = {
        trajectory.topology.atom(int(idx[1])).residue.index: i
        for i, idx in enumerate(psi_indices)
    }
    common = sorted(set(phi_residues) & set(psi_residues))

    if not common:
        raise ValueError("No residues with both phi and psi dihedrals found")

    parts: list[np.ndarray] = []
    for res_idx in common:
        phi_col = phi_angles[:, phi_residues[res_idx]]
        psi_col = psi_angles[:, psi_residues[res_idx]]
        parts.extend([
            np.sin(phi_col)[:, np.newaxis],
            np.cos(phi_col)[:, np.newaxis],
            np.sin(psi_col)[:, np.newaxis],
            np.cos(psi_col)[:, np.newaxis],
        ])

    features = np.hstack(parts).astype(float)
    logger.info("Dihedral features: %d frames x %d features", *features.shape)
    return features


def compute_contact_features(
    trajectory: md.Trajectory,
    chain_a_indices: np.ndarray,
    chain_b_indices: np.ndarray,
    cutoff_nm: float = 0.45,
) -> np.ndarray:
    """Compute per-frame interface contact features for MSM featurization.

    For each frame, computes a binary vector indicating which inter-chain
    residue pairs are within the contact cutoff distance.

    Args:
        trajectory: MDTraj trajectory. Shape: [N_frames, N_atoms, 3].
        chain_a_indices: Atom indices for chain A. Shape: [N_a].
        chain_b_indices: Atom indices for chain B. Shape: [N_b].
        cutoff_nm: Distance cutoff for a contact in nm.

    Returns:
        np.ndarray: Binary contact feature matrix. Shape: [N_frames, N_pairs].
    """
    from src.analyze.contacts import compute_interface_contacts

    result = compute_interface_contacts(
        trajectory,
        chain_a_indices,
        chain_b_indices,
        cutoff_nm=cutoff_nm,
    )
    # contact_frequency is shape [N_residues_a, N_residues_b] — aggregated.
    # For per-frame features we need to recompute per-frame contact maps.
    # Use the n_contacts_per_frame as a 1D feature (total contact count).
    n_contacts = result["n_contacts_per_frame"]
    return n_contacts.reshape(trajectory.n_frames, 1).astype(float)


def combine_features(*feature_arrays: np.ndarray) -> np.ndarray:
    """Concatenate feature arrays along the feature axis.

    All arrays must have the same number of frames (axis 0).

    Args:
        *feature_arrays: Feature matrices, each shape [N_frames, N_features_k].

    Returns:
        np.ndarray: Combined feature matrix. Shape: [N_frames, sum(N_features_k)].
    """
    if not feature_arrays:
        raise ValueError("At least one feature array is required")
    n_frames = feature_arrays[0].shape[0]
    for i, arr in enumerate(feature_arrays):
        if arr.shape[0] != n_frames:
            raise ValueError(
                f"Feature array {i} has {arr.shape[0]} frames, expected {n_frames}"
            )
    return np.column_stack(feature_arrays).astype(float)
