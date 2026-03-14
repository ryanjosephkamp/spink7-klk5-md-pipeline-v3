"""Tests for the WHAM solver and bootstrap uncertainty estimation."""

from __future__ import annotations

import numpy as np
import pytest

from src.analyze.wham import bootstrap_pmf_uncertainty, solve_wham
from src.config import BOLTZMANN_KJ, WHAMConfig


def _make_harmonic_wham_dataset(
    *,
    n_windows: int = 9,
    samples_per_window: int = 4000,
    temperature_k: float = 310.0,
    unbiased_center_nm: float = 0.0,
    unbiased_spring_kj_mol_nm2: float = 40.0,
    umbrella_spring_kj_mol_nm2: float = 300.0,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray, float, float, float]:
    """Generate synthetic umbrella samples from an analytically solvable harmonic PMF."""

    rng = np.random.default_rng(12345)
    window_centers = np.linspace(-0.6, 0.6, n_windows, dtype=float)
    spring_constants = np.full(n_windows, umbrella_spring_kj_mol_nm2, dtype=float)
    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)
    total_spring = unbiased_spring_kj_mol_nm2 + umbrella_spring_kj_mol_nm2
    variance_nm2 = 1.0 / (beta * total_spring)
    sigma_nm = float(np.sqrt(variance_nm2))

    xi_timeseries_list: list[np.ndarray] = []
    for window_center in window_centers:
        biased_mean = (
            unbiased_spring_kj_mol_nm2 * unbiased_center_nm + umbrella_spring_kj_mol_nm2 * window_center
        ) / total_spring
        xi_timeseries_list.append(rng.normal(loc=biased_mean, scale=sigma_nm, size=samples_per_window))

    return (
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        unbiased_center_nm,
        unbiased_spring_kj_mol_nm2,
    )


def test_solve_wham_recovers_harmonic_pmf_from_synthetic_biased_samples() -> None:
    """Chunk 18 gate: WHAM should recover a known harmonic free-energy profile from biased data."""

    (
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        unbiased_center_nm,
        unbiased_spring_kj_mol_nm2,
    ) = _make_harmonic_wham_dataset()
    config = WHAMConfig(tolerance=1e-6, max_iterations=25_000, n_bootstrap=16, histogram_bins=80)

    result = solve_wham(xi_timeseries_list, window_centers, spring_constants, temperature_k, config)

    assert result["converged"] is True
    assert result["xi_bins"].shape == (config.histogram_bins,)
    assert result["pmf_kj_mol"].shape == (config.histogram_bins,)
    assert result["pmf_kcal_mol"].shape == (config.histogram_bins,)
    assert result["free_energies_f_i"].shape == (window_centers.size,)
    assert result["n_iterations"] > 0

    xi_bins = result["xi_bins"]
    analytic_pmf = 0.5 * unbiased_spring_kj_mol_nm2 * (xi_bins - unbiased_center_nm) ** 2
    populated_mask = np.isfinite(result["pmf_kj_mol"]) & (np.abs(xi_bins) <= 0.5)
    shifted_numeric = result["pmf_kj_mol"][populated_mask] - np.min(result["pmf_kj_mol"][populated_mask])
    shifted_analytic = analytic_pmf[populated_mask] - np.min(analytic_pmf[populated_mask])
    rmse_kj_mol = np.sqrt(np.mean((shifted_numeric - shifted_analytic) ** 2, dtype=float))

    assert np.count_nonzero(populated_mask) >= 20
    assert rmse_kj_mol < 0.75


def test_bootstrap_pmf_uncertainty_returns_expected_shapes_and_nonnegative_spread() -> None:
    """Bootstrap uncertainty estimation should return bounded arrays with non-negative standard deviations."""

    xi_timeseries_list, window_centers, spring_constants, temperature_k, _, _ = _make_harmonic_wham_dataset(
        samples_per_window=1500,
    )
    config = WHAMConfig(tolerance=1e-6, max_iterations=20_000, n_bootstrap=12, histogram_bins=60)

    result = bootstrap_pmf_uncertainty(xi_timeseries_list, window_centers, spring_constants, temperature_k, config)

    assert result["xi_bins"].shape == (config.histogram_bins,)
    assert result["pmf_mean"].shape == (config.histogram_bins,)
    assert result["pmf_std"].shape == (config.histogram_bins,)
    assert result["pmf_bootstrap_samples"].shape == (config.n_bootstrap, config.histogram_bins)
    assert np.all(result["pmf_std"] >= 0.0)
    assert np.any(result["pmf_std"] > 0.0)


def test_solve_wham_rejects_mismatched_window_metadata() -> None:
    """Public API should reject window metadata with a length mismatch."""

    xi_timeseries_list, window_centers, spring_constants, temperature_k, _, _ = _make_harmonic_wham_dataset()
    config = WHAMConfig()

    with pytest.raises(ValueError, match="window_centers"):
        solve_wham(
            xi_timeseries_list,
            window_centers[:-1],
            spring_constants,
            temperature_k,
            config,
        )