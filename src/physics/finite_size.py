"""Finite-size electrostatic correction utilities for the SPINK7-KLK5 MD pipeline.

Implements analytical corrections for artifacts arising from periodic boundary
conditions with Particle Mesh Ewald electrostatics, following Hünenberger and
McCammon (1999) and Rocklin et al. (2013).
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

from src.config import FiniteSizeCorrectionConfig

logger = logging.getLogger(__name__)

KCAL_PER_KJ: float = 1.0 / 4.184


def hunenberger_mccammon_correction(
    net_charge_e: float,
    box_length_nm: float,
    config: FiniteSizeCorrectionConfig | None = None,
) -> float:
    """Compute the Hünenberger-McCammon finite-size correction.

    Implements the leading-order correction for the interaction of a net
    charge with its periodic images in PME electrostatics:

        ΔG_HM = -(Q² · ξ_EW · k_e) / (2 · ε_s · L)

    Args:
        net_charge_e: Net charge of the solute in elementary charges.
        box_length_nm: Cubic box edge length in nm.
        config: Correction parameters. Uses defaults if None.

    Returns:
        Correction energy in kJ/mol.

    Raises:
        ValueError: If box_length_nm is not positive.
    """
    if box_length_nm <= 0.0:
        raise ValueError(f"box_length_nm must be positive, got {box_length_nm}")

    if config is None:
        config = FiniteSizeCorrectionConfig()

    return -(
        net_charge_e**2
        * config.ewald_self_energy_constant
        * config.coulomb_constant_kj_nm_per_mol_e2
    ) / (2.0 * config.solvent_dielectric * box_length_nm)


def finite_size_correction(
    net_charge_e: float,
    box_vectors_nm: np.ndarray,
    config: FiniteSizeCorrectionConfig | None = None,
) -> dict[str, float]:
    """Compute finite-size electrostatic correction for arbitrary box geometries.

    For non-cubic boxes, the effective box length is the geometric mean
    of the three box edge lengths.

    Args:
        net_charge_e: Net charge of the solute in elementary charges.
        box_vectors_nm: Box vectors as a (3, 3) array in nm.
        config: Correction parameters.

    Returns:
        Dictionary with 'correction_kj_mol', 'effective_box_length_nm',
        and 'correction_kcal_mol'.
    """
    box_vectors_nm = np.asarray(box_vectors_nm, dtype=float)
    if box_vectors_nm.shape != (3, 3):
        raise ValueError(f"box_vectors_nm must have shape (3, 3), got {box_vectors_nm.shape}")

    edge_lengths = np.array([
        np.linalg.norm(box_vectors_nm[0]),
        np.linalg.norm(box_vectors_nm[1]),
        np.linalg.norm(box_vectors_nm[2]),
    ])
    effective_length = float(np.cbrt(np.prod(edge_lengths)))

    correction_kj = hunenberger_mccammon_correction(net_charge_e, effective_length, config)

    return {
        "correction_kj_mol": correction_kj,
        "effective_box_length_nm": effective_length,
        "correction_kcal_mol": correction_kj * KCAL_PER_KJ,
    }


def extrapolate_to_infinite_box(
    box_lengths_nm: np.ndarray,
    energies_kj_mol: np.ndarray,
) -> dict[str, float]:
    """Extrapolate energies to the infinite-box-size limit using 1/L fitting.

    Fits the model: E(L) = E_inf + a/L + b/L^3

    Args:
        box_lengths_nm: Array of box edge lengths in nm.
        energies_kj_mol: Corresponding energies in kJ/mol.

    Returns:
        Dictionary with 'energy_infinite_kj_mol', 'slope_1_over_L',
        'r_squared', and 'residuals_kj_mol'.

    Raises:
        ValueError: If fewer than 3 data points are provided.
    """
    box_lengths_nm = np.asarray(box_lengths_nm, dtype=float)
    energies_kj_mol = np.asarray(energies_kj_mol, dtype=float)

    if len(box_lengths_nm) < 3:
        raise ValueError("At least 3 data points are required for extrapolation")
    if len(box_lengths_nm) != len(energies_kj_mol):
        raise ValueError("box_lengths_nm and energies_kj_mol must have the same length")

    inv_l = 1.0 / box_lengths_nm
    inv_l3 = inv_l**3

    # Design matrix: [1, 1/L, 1/L^3]
    design = np.column_stack([np.ones_like(inv_l), inv_l, inv_l3])
    coeffs, residuals_arr, _, _ = np.linalg.lstsq(design, energies_kj_mol, rcond=None)

    e_inf = float(coeffs[0])
    slope = float(coeffs[1])
    fitted = design @ coeffs
    residuals = energies_kj_mol - fitted

    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((energies_kj_mol - np.mean(energies_kj_mol))**2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0

    return {
        "energy_infinite_kj_mol": e_inf,
        "slope_1_over_L": slope,
        "r_squared": r_squared,
        "residuals_kj_mol": residuals.tolist(),
    }


def compute_solute_net_charge(
    simulation: Any,
    solute_indices: list[int],
) -> float:
    """Compute the net charge of a subset of atoms in elementary charges.

    Args:
        simulation: OpenMM Simulation object.
        solute_indices: Atom indices of the solute.

    Returns:
        Net charge in elementary charges.
    """
    from openmm import NonbondedForce

    system = simulation.system
    nb_force = None
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nb_force = force
            break

    if nb_force is None:
        logger.warning("No NonbondedForce found; assuming net charge = 0.0")
        return 0.0

    total_charge = 0.0
    for idx in solute_indices:
        charge, _, _ = nb_force.getParticleParameters(idx)
        total_charge += charge.value_in_unit(charge.unit)

    return float(total_charge)


def run_box_size_convergence_study(
    energies_by_box_size: dict[float, float],
) -> dict[str, Any]:
    """Analyze electrostatic energy convergence with box size.

    Args:
        energies_by_box_size: Mapping from box edge length (nm) to computed
            electrostatic energy (kJ/mol). Must contain at least 3 entries.

    Returns:
        Dictionary with:
        - 'energy_infinite_kj_mol': Extrapolated infinite-box energy
        - 'correction_at_default_padding_kj_mol': Difference from infinite-box at 1.2 nm padding
        - 'r_squared': Quality of 1/L fit
        - 'scaling_exponent': Fitted dominant exponent (should be ~1.0)
    """
    if len(energies_by_box_size) < 3:
        raise ValueError("At least 3 box sizes are required for convergence study")

    box_lengths = np.array(sorted(energies_by_box_size.keys()))
    energies = np.array([energies_by_box_size[L] for L in box_lengths])

    result = extrapolate_to_infinite_box(box_lengths, energies)

    # Estimate the scaling exponent from log-log regression of (E - E_inf) vs L
    e_inf = result["energy_infinite_kj_mol"]
    deviations = np.abs(energies - e_inf)
    mask = deviations > 0.0
    if np.sum(mask) >= 2:
        log_l = np.log(box_lengths[mask])
        log_dev = np.log(deviations[mask])
        poly = np.polyfit(log_l, log_dev, 1)
        scaling_exponent = float(-poly[0])
    else:
        scaling_exponent = 1.0

    # Compute correction at default padding (smallest box in the study)
    smallest_box_energy = energies[0]
    correction_at_default = smallest_box_energy - e_inf

    return {
        "energy_infinite_kj_mol": e_inf,
        "correction_at_default_padding_kj_mol": float(correction_at_default),
        "r_squared": result["r_squared"],
        "scaling_exponent": scaling_exponent,
    }
