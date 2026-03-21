"""Alchemical free energy analysis using BAR and MBAR estimators.

Computes alchemical free energy differences from lambda-windowed
simulations and assembles the thermodynamic cycle for DDG calculation.
"""

from __future__ import annotations

import logging

import numpy as np

from src.config import BOLTZMANN_KJ, KCAL_TO_KJ

logger = logging.getLogger(__name__)


def compute_delta_g_bar(
    energy_forward: np.ndarray,
    energy_reverse: np.ndarray,
    temperature_k: float,
) -> float:
    """Compute free energy difference via Bennett Acceptance Ratio.

    Args:
        energy_forward: Delta-U values (kJ/mol) from lambda_i configurations
            evaluated at lambda_{i+1}. Shape: [N_samples].
        energy_reverse: Delta-U values (kJ/mol) from lambda_{i+1}
            configurations evaluated at lambda_i. Shape: [N_samples].
        temperature_k: Temperature in Kelvin.

    Returns:
        float: Free energy difference in kJ/mol.

    Raises:
        ValueError: If inputs are empty or temperature is non-positive.
    """
    if energy_forward.size == 0 or energy_reverse.size == 0:
        raise ValueError("Energy arrays must be non-empty")
    if temperature_k <= 0.0:
        raise ValueError("temperature_k must be positive")

    from pymbar import BAR as pymbar_BAR

    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)
    result = pymbar_BAR(energy_forward * beta, energy_reverse * beta)
    return float(result["Delta_f"]) / beta


def compute_delta_g_mbar(
    energy_matrix: np.ndarray,
    n_samples_per_state: np.ndarray,
    temperature_k: float,
) -> dict[str, float | np.ndarray]:
    """Compute alchemical free energies via MBAR.

    Args:
        energy_matrix: Reduced potential energy matrix (dimensionless).
            Shape: [K, N_total] where K is the number of lambda states
            and N_total is the total number of samples across all states.
        n_samples_per_state: Number of samples from each state. Shape: [K].
        temperature_k: Temperature in Kelvin.

    Returns:
        dict with keys:
            - delta_g_kj_mol: total alchemical free energy (kJ/mol)
            - delta_g_kcal_mol: total alchemical free energy (kcal/mol)
            - delta_g_std_kj_mol: uncertainty (kJ/mol)
            - delta_g_std_kcal_mol: uncertainty (kcal/mol)
            - free_energy_profile_kj_mol: per-window free energies (kJ/mol)

    Raises:
        ValueError: If inputs have inconsistent shapes or temperature
            is non-positive.
    """
    if temperature_k <= 0.0:
        raise ValueError("temperature_k must be positive")
    if energy_matrix.ndim != 2:
        raise ValueError("energy_matrix must be 2-D")
    if n_samples_per_state.ndim != 1:
        raise ValueError("n_samples_per_state must be 1-D")
    if energy_matrix.shape[0] != len(n_samples_per_state):
        raise ValueError(
            "energy_matrix rows must match length of n_samples_per_state"
        )
    if energy_matrix.shape[1] != int(n_samples_per_state.sum()):
        raise ValueError(
            "energy_matrix columns must equal total samples"
        )

    from pymbar import MBAR

    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)
    mbar = MBAR(energy_matrix, n_samples_per_state)
    free_energies = mbar.compute_free_energy_differences()

    delta_g_kj = float(free_energies["Delta_f"][0, -1]) / beta
    delta_g_std_kj = float(free_energies["dDelta_f"][0, -1]) / beta

    profile = np.array(free_energies["Delta_f"][0, :], dtype=float) / beta

    return {
        "delta_g_kj_mol": delta_g_kj,
        "delta_g_kcal_mol": delta_g_kj / KCAL_TO_KJ,
        "delta_g_std_kj_mol": delta_g_std_kj,
        "delta_g_std_kcal_mol": delta_g_std_kj / KCAL_TO_KJ,
        "free_energy_profile_kj_mol": profile,
    }


def compute_delta_delta_g(
    delta_g_complex: dict[str, float],
    delta_g_free: dict[str, float],
) -> dict[str, float]:
    """Compute the relative binding free energy DDG.

    DDG = DG_alch_complex - DG_alch_free

    A positive DDG indicates that the mutation destabilizes binding.

    Args:
        delta_g_complex: MBAR result dict for the complex leg.
            Must contain "delta_g_kj_mol" and "delta_g_std_kj_mol".
        delta_g_free: MBAR result dict for the free-solution leg.
            Must contain "delta_g_kj_mol" and "delta_g_std_kj_mol".

    Returns:
        dict with keys:
            - delta_delta_g_kj_mol: DDG in kJ/mol
            - delta_delta_g_kcal_mol: DDG in kcal/mol
            - delta_delta_g_std_kj_mol: propagated uncertainty (kJ/mol)
            - delta_delta_g_std_kcal_mol: propagated uncertainty (kcal/mol)
    """
    ddg_kj = delta_g_complex["delta_g_kj_mol"] - delta_g_free["delta_g_kj_mol"]
    ddg_std_kj = np.sqrt(
        delta_g_complex["delta_g_std_kj_mol"] ** 2
        + delta_g_free["delta_g_std_kj_mol"] ** 2
    )
    return {
        "delta_delta_g_kj_mol": ddg_kj,
        "delta_delta_g_kcal_mol": ddg_kj / KCAL_TO_KJ,
        "delta_delta_g_std_kj_mol": ddg_std_kj,
        "delta_delta_g_std_kcal_mol": ddg_std_kj / KCAL_TO_KJ,
    }
