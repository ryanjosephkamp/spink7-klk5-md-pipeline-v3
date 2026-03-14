"""Jarzynski free-energy estimators for nonequilibrium SMD work data."""

from __future__ import annotations

import numpy as np

from src.config import BOLTZMANN_KJ, KCAL_TO_KJ


def _validate_work_values(work_values: np.ndarray) -> np.ndarray:
    """Validate a one-dimensional work array at the public API boundary."""

    values = np.asarray(work_values, dtype=float)
    assert values.ndim == 1, "work_values must have shape [N_replicates]"
    if values.ndim != 1 or values.size == 0:
        raise ValueError("work_values must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("work_values must contain only finite values")
    return values


def _validate_temperature(temperature_k: float) -> float:
    """Validate the thermodynamic temperature input."""

    temperature = float(temperature_k)
    if temperature <= 0.0:
        raise ValueError("temperature_k must be positive")
    return temperature


def _beta(temperature_k: float) -> float:
    """Return the inverse thermal energy in mol/kJ."""

    return 1.0 / (BOLTZMANN_KJ * temperature_k)


def _jarzynski_exact_kj_mol(work_values: np.ndarray, temperature_k: float) -> float:
    """Compute the numerically stable exact Jarzynski estimator in kJ/mol."""

    # Numerically stable Jarzynski: shift by max(-beta * W) before
    # exponentiation to avoid overflow in the exponential average.
    beta = _beta(temperature_k)
    scaled_work = -beta * work_values
    maximum = np.max(scaled_work)
    log_average = maximum + np.log(np.mean(np.exp(scaled_work - maximum), dtype=float))
    return -(1.0 / beta) * log_average


def jarzynski_free_energy(
    work_values: np.ndarray,
    temperature_k: float,
) -> dict[str, float]:
    """Estimate free-energy differences from nonequilibrium work values.

    Args:
        work_values: Total work values in kJ/mol. Shape: [N_replicates].
        temperature_k: Temperature in Kelvin.

    Returns:
        dict[str, float]: Exact Jarzynski and second-order cumulant estimates.
    """

    values = _validate_work_values(work_values)
    temperature = _validate_temperature(temperature_k)
    beta = _beta(temperature)

    mean_work_kj_mol = float(np.mean(values, dtype=float))
    work_variance_kj_mol2 = float(np.var(values, ddof=0, dtype=float))
    delta_g_kj_mol = float(_jarzynski_exact_kj_mol(values, temperature))
    delta_g_cumulant2_kj_mol = float(mean_work_kj_mol - 0.5 * beta * work_variance_kj_mol2)

    return {
        "delta_g_kj_mol": delta_g_kj_mol,
        "delta_g_kcal_mol": delta_g_kj_mol / KCAL_TO_KJ,
        "mean_work_kj_mol": mean_work_kj_mol,
        "work_variance_kj_mol2": work_variance_kj_mol2,
        "delta_g_cumulant2_kj_mol": delta_g_cumulant2_kj_mol,
        "delta_g_cumulant2_kcal_mol": delta_g_cumulant2_kj_mol / KCAL_TO_KJ,
    }


def evaluate_convergence(
    work_values: np.ndarray,
    temperature_k: float,
    n_subsets: int = 10,
) -> dict[str, np.ndarray]:
    """Evaluate Jarzynski estimator convergence as a function of trajectory count.

    Args:
        work_values: Total work values in kJ/mol. Shape: [N_replicates].
        temperature_k: Temperature in Kelvin.
        n_subsets: Number of prefix subset sizes to evaluate.

    Returns:
        dict[str, np.ndarray]: Prefix sizes, exact Jarzynski estimates, and bootstrap SE.
    """

    values = _validate_work_values(work_values)
    temperature = _validate_temperature(temperature_k)

    subset_count = int(n_subsets)
    if subset_count <= 0:
        raise ValueError("n_subsets must be positive")

    # Evaluate convergence by computing the Jarzynski estimate for
    # growing prefix subsets and estimating bootstrap standard errors.
    subset_sizes = np.unique(np.linspace(1, values.size, subset_count, dtype=int))
    delta_g_vs_n = np.empty(subset_sizes.size, dtype=float)
    std_vs_n = np.empty(subset_sizes.size, dtype=float)
    rng = np.random.default_rng(0)
    n_bootstrap = 200

    for subset_index, subset_size in enumerate(subset_sizes):
        subset_values = values[:subset_size]
        delta_g_vs_n[subset_index] = _jarzynski_exact_kj_mol(subset_values, temperature)

        bootstrap_estimates = np.empty(n_bootstrap, dtype=float)
        for bootstrap_index in range(n_bootstrap):
            sampled_indices = rng.integers(0, subset_size, size=subset_size)
            bootstrap_estimates[bootstrap_index] = _jarzynski_exact_kj_mol(subset_values[sampled_indices], temperature)
        std_vs_n[subset_index] = np.std(bootstrap_estimates, ddof=1, dtype=float)

    return {
        "subset_sizes": subset_sizes,
        "delta_g_vs_n": delta_g_vs_n,
        "std_vs_n": std_vs_n,
    }