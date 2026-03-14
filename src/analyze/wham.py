"""Weighted histogram analysis method (WHAM) solvers."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from src import PhysicalValidityError
from src.config import BOLTZMANN_KJ, KCAL_TO_KJ, WHAMConfig


def _validate_timeseries_list(xi_timeseries_list: list[np.ndarray]) -> list[np.ndarray]:
    """Validate the public WHAM timeseries input contract."""

    if len(xi_timeseries_list) == 0:
        raise ValueError("xi_timeseries_list must contain at least one window")

    validated: list[np.ndarray] = []
    for window_index, samples in enumerate(xi_timeseries_list):
        array = np.asarray(samples, dtype=float)
        assert array.ndim == 1, "each xi timeseries must have shape [N_samples]"
        if array.ndim != 1 or array.size == 0:
            raise ValueError(f"xi_timeseries_list[{window_index}] must be a non-empty one-dimensional array")
        validated.append(array)
    return validated


def _validate_window_parameters(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate that window metadata matches the number of timeseries windows."""

    centers = np.asarray(window_centers, dtype=float)
    springs = np.asarray(spring_constants, dtype=float)
    n_windows = len(xi_timeseries_list)

    assert centers.ndim == 1, "window_centers must have shape [M]"
    assert springs.ndim == 1, "spring_constants must have shape [M]"
    if centers.ndim != 1 or centers.size != n_windows:
        raise ValueError("window_centers must be a one-dimensional array with one value per window")
    if springs.ndim != 1 or springs.size != n_windows:
        raise ValueError("spring_constants must be a one-dimensional array with one value per window")
    if np.any(springs <= 0.0):
        raise ValueError("spring_constants must be strictly positive")

    return centers, springs


def _validate_wham_inputs(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: WHAMConfig,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray, float]:
    """Validate all public solver inputs."""

    validated_timeseries = _validate_timeseries_list(xi_timeseries_list)
    centers, springs = _validate_window_parameters(validated_timeseries, window_centers, spring_constants)

    temperature = float(temperature_k)
    if temperature <= 0.0:
        raise ValueError("temperature_k must be positive")
    if config.histogram_bins < 2:
        raise ValueError("config.histogram_bins must be at least 2")
    if config.tolerance <= 0.0:
        raise ValueError("config.tolerance must be positive")
    if config.max_iterations <= 0:
        raise ValueError("config.max_iterations must be positive")
    if config.n_bootstrap <= 0:
        raise ValueError("config.n_bootstrap must be positive")

    return validated_timeseries, centers, springs, temperature


def _common_bin_edges(xi_timeseries_list: list[np.ndarray], n_bins: int) -> np.ndarray:
    """Build a shared histogram grid for all windows."""

    minima = np.asarray([np.min(samples) for samples in xi_timeseries_list], dtype=float)
    maxima = np.asarray([np.max(samples) for samples in xi_timeseries_list], dtype=float)
    lower = float(np.min(minima))
    upper = float(np.max(maxima))
    if np.isclose(lower, upper):
        lower -= 1e-3
        upper += 1e-3
    padding = max(1e-6, 0.02 * (upper - lower))
    return np.linspace(lower - padding, upper + padding, n_bins + 1, dtype=float)


def _histogram_counts(xi_timeseries_list: list[np.ndarray], bin_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Histogram all windows onto a shared bin grid."""

    n_windows = len(xi_timeseries_list)
    n_bins = bin_edges.size - 1
    counts = np.zeros((n_windows, n_bins), dtype=float)
    n_samples = np.zeros(n_windows, dtype=float)
    for window_index, samples in enumerate(xi_timeseries_list):
        counts[window_index], _ = np.histogram(samples, bins=bin_edges)
        n_samples[window_index] = float(samples.size)
    return counts, n_samples


def _logsumexp(values: np.ndarray, axis: int) -> np.ndarray:
    """Compute a numerically stable log-sum-exp reduction.

    Subtracts the per-slice maximum before exponentiation to prevent
    overflow, then adds it back in log-space.
    """

    maximum = np.max(values, axis=axis, keepdims=True)
    stabilized = values - maximum
    summed = np.sum(np.exp(stabilized), axis=axis, keepdims=True)
    return np.squeeze(maximum + np.log(summed), axis=axis)


def _compute_probability_density(
    counts: np.ndarray,
    n_samples: np.ndarray,
    free_energies: np.ndarray,
    beta: float,
    bias_energies: np.ndarray,
    bin_width: float,
) -> np.ndarray:
    """Evaluate the WHAM unbiased probability density on the shared grid."""

    numerator = np.sum(counts, axis=0, dtype=float)
    log_denominator = _logsumexp(
        np.log(n_samples)[:, np.newaxis] + beta * (free_energies[:, np.newaxis] - bias_energies),
        axis=0,
    )
    probability_density = np.zeros_like(numerator, dtype=float)
    nonzero_mask = numerator > 0.0
    probability_density[nonzero_mask] = numerator[nonzero_mask] / np.exp(log_denominator[nonzero_mask])

    normalization = np.sum(probability_density, dtype=float) * bin_width
    if normalization <= 0.0:
        raise PhysicalValidityError("WHAM produced a non-normalizable probability density")
    probability_density /= normalization
    return probability_density


def _update_free_energies(
    probability_density: np.ndarray,
    beta: float,
    bias_energies: np.ndarray,
    bin_width: float,
) -> np.ndarray:
    """Update WHAM free-energy offsets from the current unbiased density."""

    weighted_probability = probability_density * bin_width
    log_weighted_probability = np.full_like(weighted_probability, -np.inf, dtype=float)
    positive_mask = weighted_probability > 0.0
    log_weighted_probability[positive_mask] = np.log(weighted_probability[positive_mask])
    log_integrals = _logsumexp(log_weighted_probability[np.newaxis, :] - beta * bias_energies, axis=1)
    return -(1.0 / beta) * log_integrals


def _shift_pmf_to_dissociated_state(pmf_kj_mol: np.ndarray) -> np.ndarray:
    """Set the final finite PMF bin to zero as the dissociated-state reference."""

    shifted = np.asarray(pmf_kj_mol, dtype=float).copy()
    finite_indices = np.flatnonzero(np.isfinite(shifted))
    if finite_indices.size == 0:
        raise PhysicalValidityError("WHAM produced no finite PMF bins")
    shifted -= shifted[finite_indices[-1]]
    return shifted


def _solve_wham_core(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: WHAMConfig,
) -> dict[str, np.ndarray | int | bool]:
    """Internal WHAM solver implementation shared by the bootstrap routine."""

    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)
    bin_edges = _common_bin_edges(xi_timeseries_list, config.histogram_bins)
    xi_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = float(bin_edges[1] - bin_edges[0])
    counts, n_samples = _histogram_counts(xi_timeseries_list, bin_edges)
    # Precompute the bias energy matrix: U_bias[i, k] = 0.5 * k_i * (xi_k - xi_0_i)^2.
    bias_energies = 0.5 * spring_constants[:, np.newaxis] * (xi_bins[np.newaxis, :] - window_centers[:, np.newaxis]) ** 2

    # Iterative WHAM self-consistency loop: alternate between updating the
    # unbiased probability density and the per-window free-energy offsets
    # until the maximum absolute change falls below the tolerance.
    free_energies = np.zeros(window_centers.size, dtype=float)
    n_iterations = 0
    converged = False
    for iteration_index in range(1, config.max_iterations + 1):
        probability_density = _compute_probability_density(counts, n_samples, free_energies, beta, bias_energies, bin_width)
        updated_free_energies = _update_free_energies(probability_density, beta, bias_energies, bin_width)
        delta = np.max(np.abs(updated_free_energies - free_energies))
        free_energies = updated_free_energies
        n_iterations = iteration_index
        if delta < config.tolerance:
            converged = True
            break

    if not converged:
        raise PhysicalValidityError(
            f"WHAM failed IV-9 convergence after {config.max_iterations} iterations; max|df| remained above tolerance"
        )

    probability_density = _compute_probability_density(counts, n_samples, free_energies, beta, bias_energies, bin_width)
    pmf_kj_mol = np.full(config.histogram_bins, np.inf, dtype=float)
    positive_mask = probability_density > 0.0
    pmf_kj_mol[positive_mask] = -(1.0 / beta) * np.log(probability_density[positive_mask])
    pmf_kj_mol = _shift_pmf_to_dissociated_state(pmf_kj_mol)

    return {
        "xi_bins": xi_bins,
        "pmf_kj_mol": pmf_kj_mol,
        "pmf_kcal_mol": pmf_kj_mol / KCAL_TO_KJ,
        "free_energies_f_i": free_energies,
        "n_iterations": n_iterations,
        "converged": True,
    }


def solve_wham(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: WHAMConfig,
) -> dict[str, np.ndarray | int | bool]:
    """Solve the WHAM equations for an unbiased PMF along a shared reaction coordinate."""

    validated_timeseries, centers, springs, temperature = _validate_wham_inputs(
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        config,
    )
    return _solve_wham_core(validated_timeseries, centers, springs, temperature, config)


def _bootstrap_resample_window(samples: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Block-bootstrap a single window timeseries with bounded memory use."""

    n_samples = samples.size
    block_size = max(1, int(round(np.sqrt(n_samples))))
    if block_size >= n_samples:
        return samples.copy()

    n_blocks = int(np.ceil(n_samples / block_size))
    max_start = n_samples - block_size
    start_indices = rng.integers(0, max_start + 1, size=n_blocks)
    resampled = np.empty(n_samples, dtype=float)
    write_offset = 0
    for start_index in start_indices:
        end_offset = min(write_offset + block_size, n_samples)
        resampled[write_offset:end_offset] = samples[start_index : start_index + (end_offset - write_offset)]
        write_offset = end_offset
        if write_offset >= n_samples:
            break
    return resampled


def bootstrap_pmf_uncertainty(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: WHAMConfig,
) -> dict[str, np.ndarray]:
    """Estimate WHAM PMF uncertainty from bootstrap-resampled window trajectories."""

    validated_timeseries, centers, springs, temperature = _validate_wham_inputs(
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        config,
    )

    base_result = _solve_wham_core(validated_timeseries, centers, springs, temperature, config)
    # Re-solve WHAM on each bootstrap-resampled dataset to estimate
    # the uncertainty (standard deviation) of the PMF at each bin.
    pmf_bootstrap_samples = np.empty((config.n_bootstrap, config.histogram_bins), dtype=float)
    rng = np.random.default_rng(0)

    bootstrap_config = replace(config)
    for bootstrap_index in range(config.n_bootstrap):
        resampled_timeseries = [
            _bootstrap_resample_window(samples, rng)
            for samples in validated_timeseries
        ]
        bootstrap_result = _solve_wham_core(resampled_timeseries, centers, springs, temperature, bootstrap_config)
        pmf_bootstrap_samples[bootstrap_index, :] = bootstrap_result["pmf_kj_mol"]

    finite_mask = np.isfinite(pmf_bootstrap_samples)
    finite_counts = np.sum(finite_mask, axis=0)
    pmf_mean = np.zeros(config.histogram_bins, dtype=float)
    pmf_std = np.zeros(config.histogram_bins, dtype=float)
    for bin_index in range(config.histogram_bins):
        finite_values = pmf_bootstrap_samples[finite_mask[:, bin_index], bin_index]
        if finite_values.size != 0:
            pmf_mean[bin_index] = np.mean(finite_values, dtype=float)
        if finite_values.size > 1:
            pmf_std[bin_index] = np.std(finite_values, ddof=1, dtype=float)

    return {
        "pmf_mean": pmf_mean,
        "pmf_std": pmf_std,
        "pmf_bootstrap_samples": pmf_bootstrap_samples,
        "xi_bins": base_result["xi_bins"],
    }