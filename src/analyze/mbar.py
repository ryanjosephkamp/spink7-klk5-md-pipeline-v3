"""Multistate Bennett Acceptance Ratio (MBAR) solver for PMF estimation.

Implements MBAR as an alternative to WHAM for potential of mean force
reconstruction from umbrella sampling data.  Backed by the ``pymbar``
library (>= 4.0).
"""

from __future__ import annotations

import logging

import numpy as np
import pymbar
from scipy.special import logsumexp as _scipy_logsumexp

from src import PhysicalValidityError
from src.config import BOLTZMANN_KJ, KCAL_TO_KJ, MBARConfig

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def _validate_mbar_inputs(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: MBARConfig,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray, float]:
    """Validate all public MBAR solver inputs."""

    if len(xi_timeseries_list) == 0:
        raise ValueError("xi_timeseries_list must contain at least one window")

    validated: list[np.ndarray] = []
    for i, samples in enumerate(xi_timeseries_list):
        arr = np.asarray(samples, dtype=float)
        if arr.ndim != 1 or arr.size == 0:
            raise ValueError(
                f"xi_timeseries_list[{i}] must be a non-empty 1-D array"
            )
        validated.append(arr)

    centers = np.asarray(window_centers, dtype=float)
    springs = np.asarray(spring_constants, dtype=float)
    n_windows = len(validated)

    if centers.ndim != 1 or centers.size != n_windows:
        raise ValueError(
            "window_centers must be a 1-D array with one value per window"
        )
    if springs.ndim != 1 or springs.size != n_windows:
        raise ValueError(
            "spring_constants must be a 1-D array with one value per window"
        )
    if np.any(springs <= 0.0):
        raise ValueError("spring_constants must be strictly positive")

    temperature = float(temperature_k)
    if temperature <= 0.0:
        raise ValueError("temperature_k must be positive")
    if config.n_pmf_bins < 2:
        raise ValueError("config.n_pmf_bins must be at least 2")

    return validated, centers, springs, temperature


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_reduced_potential_matrix(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    beta: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the (K x N) reduced potential matrix u_kn.

    Returns
    -------
    u_kn : ndarray, shape (K, N)
    N_k  : ndarray, shape (K,)
    xi_all : ndarray, shape (N,)
    """

    xi_all = np.concatenate(xi_timeseries_list)
    N_k = np.array([ts.size for ts in xi_timeseries_list], dtype=int)
    K = len(xi_timeseries_list)
    N = xi_all.size

    u_kn = np.zeros((K, N), dtype=float)
    for k in range(K):
        u_kn[k] = (
            0.5 * beta * float(spring_constants[k])
            * (xi_all - float(window_centers[k])) ** 2
        )

    return u_kn, N_k, xi_all


def _shift_pmf_to_dissociated_state(pmf_kj_mol: np.ndarray) -> np.ndarray:
    """Set the last finite PMF bin to zero as the dissociated-state reference."""

    shifted = pmf_kj_mol.copy()
    finite = np.flatnonzero(np.isfinite(shifted))
    if finite.size == 0:
        raise PhysicalValidityError("MBAR produced no finite PMF bins")
    shifted -= shifted[finite[-1]]
    return shifted


def _bootstrap_resample_window(
    samples: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Block-bootstrap a single window timeseries with bounded memory use."""

    n = samples.size
    block_size = max(1, int(round(np.sqrt(n))))
    if block_size >= n:
        return samples.copy()

    n_blocks = int(np.ceil(n / block_size))
    max_start = n - block_size
    starts = rng.integers(0, max_start + 1, size=n_blocks)
    resampled = np.empty(n, dtype=float)
    offset = 0
    for start in starts:
        end = min(offset + block_size, n)
        resampled[offset:end] = samples[start : start + (end - offset)]
        offset = end
        if offset >= n:
            break
    return resampled


# ---------------------------------------------------------------------------
# Core solver
# ---------------------------------------------------------------------------

def _solve_mbar_core(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: MBARConfig,
) -> dict[str, np.ndarray | bool]:
    """Internal MBAR solver shared by the public API and bootstrap routine."""

    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)

    u_kn, N_k, xi_all = _build_reduced_potential_matrix(
        xi_timeseries_list, window_centers, spring_constants, beta,
    )
    N_total = xi_all.size

    # --- Solve MBAR ---
    try:
        mbar = pymbar.MBAR(
            u_kn,
            N_k,
            solver_protocol=config.solver_protocol,
            relative_tolerance=config.relative_tolerance,
            maximum_iterations=config.maximum_iterations,
        )
    except Exception as exc:
        raise PhysicalValidityError(
            f"MBAR failed IV-9 convergence: {exc}"
        ) from exc

    f_k_kj_mol = mbar.f_k / beta  # dimensionless -> kJ/mol

    # --- PMF bin grid ---
    xi_range = float(xi_all.max() - xi_all.min())
    padding = max(1e-6, 0.02 * xi_range)
    bin_edges = np.linspace(
        xi_all.min() - padding,
        xi_all.max() + padding,
        config.n_pmf_bins + 1,
    )
    xi_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    bin_indices = np.clip(
        np.digitize(xi_all, bin_edges) - 1, 0, config.n_pmf_bins - 1,
    )

    # Identify populated bins (bins with at least one sample).
    bin_counts = np.bincount(bin_indices, minlength=config.n_pmf_bins)
    populated = np.flatnonzero(bin_counts > 0)

    if populated.size == 0:
        raise PhysicalValidityError("MBAR: no bins contain samples")

    # --- PMF via MBAR unbiased weights (bin-free reweighting) ---
    # Compute log weights for the unbiased state (u = 0):
    #   log w_n = -logsumexp_k(log N_k + f_k - u_kn[k, n])
    f_k = mbar.f_k  # dimensionless free energies, shape (K,)
    log_nk = np.log(N_k.astype(float))
    log_denom = log_nk[:, np.newaxis] + f_k[:, np.newaxis] - u_kn  # (K, N)
    log_w_n = -_scipy_logsumexp(log_denom, axis=0)  # (N,)
    log_w_n -= _scipy_logsumexp(log_w_n)  # normalize
    w_n = np.exp(log_w_n)

    # Weighted histogram -> probability per bin -> PMF
    weighted_counts = np.zeros(config.n_pmf_bins, dtype=float)
    for j in populated:
        weighted_counts[j] = float(np.sum(w_n[bin_indices == j]))

    pmf_kj_mol = np.full(config.n_pmf_bins, np.inf, dtype=float)
    nonzero = weighted_counts > 0.0
    pmf_kj_mol[nonzero] = -(1.0 / beta) * np.log(weighted_counts[nonzero])

    # --- Per-bin uncertainty via compute_perturbed_free_energies ---
    n_pop = populated.size
    u_ln = np.full((n_pop, N_total), 1e10, dtype=float)
    for idx, b in enumerate(populated):
        u_ln[idx, bin_indices == b] = 0.0

    perturbed = mbar.compute_perturbed_free_energies(
        u_ln, compute_uncertainty=True,
    )
    d_delta_f = perturbed["dDelta_f"][0]  # shape (n_pop,)

    pmf_uncertainty = np.full(config.n_pmf_bins, np.inf, dtype=float)
    pmf_uncertainty[populated] = d_delta_f / beta

    pmf_kj_mol = _shift_pmf_to_dissociated_state(pmf_kj_mol)

    return {
        "xi_bins": xi_bins,
        "pmf_kj_mol": pmf_kj_mol,
        "pmf_kcal_mol": pmf_kj_mol / KCAL_TO_KJ,
        "pmf_uncertainty_kj_mol": pmf_uncertainty,
        "free_energies_f_k": f_k_kj_mol,
        "converged": True,
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def solve_mbar(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: MBARConfig,
) -> dict[str, np.ndarray | bool]:
    """Solve the MBAR equations for an unbiased PMF along a shared reaction coordinate.

    Accepts the same data format as ``solve_wham`` for direct cross-validation.
    """

    validated, centers, springs, temperature = _validate_mbar_inputs(
        xi_timeseries_list, window_centers, spring_constants,
        temperature_k, config,
    )
    return _solve_mbar_core(validated, centers, springs, temperature, config)


def bootstrap_mbar_uncertainty(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    spring_constants: np.ndarray,
    temperature_k: float,
    config: MBARConfig,
) -> dict[str, np.ndarray]:
    """Estimate MBAR PMF uncertainty from bootstrap-resampled window trajectories."""

    validated, centers, springs, temperature = _validate_mbar_inputs(
        xi_timeseries_list, window_centers, spring_constants,
        temperature_k, config,
    )

    base = _solve_mbar_core(validated, centers, springs, temperature, config)

    pmf_samples = np.empty(
        (config.n_bootstrap, config.n_pmf_bins), dtype=float,
    )
    rng = np.random.default_rng(0)

    for i in range(config.n_bootstrap):
        resampled = [_bootstrap_resample_window(s, rng) for s in validated]
        try:
            boot = _solve_mbar_core(
                resampled, centers, springs, temperature, config,
            )
            pmf_samples[i] = boot["pmf_kj_mol"]
        except (PhysicalValidityError, Exception):
            pmf_samples[i] = np.nan

    finite_mask = np.isfinite(pmf_samples)
    pmf_mean = np.zeros(config.n_pmf_bins, dtype=float)
    pmf_std = np.zeros(config.n_pmf_bins, dtype=float)
    for j in range(config.n_pmf_bins):
        vals = pmf_samples[finite_mask[:, j], j]
        if vals.size > 0:
            pmf_mean[j] = float(np.mean(vals))
        if vals.size > 1:
            pmf_std[j] = float(np.std(vals, ddof=1))

    return {
        "pmf_mean": pmf_mean,
        "pmf_std": pmf_std,
        "pmf_bootstrap_samples": pmf_samples,
        "xi_bins": base["xi_bins"],
    }
