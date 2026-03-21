"""Tests for the MBAR solver and bootstrap uncertainty estimation."""

from __future__ import annotations

import numpy as np
import pytest

from tests.test_wham import _make_harmonic_wham_dataset


# ---------- L-11 Step 2: solve_mbar returns expected keys and shapes ----------


def test_solve_mbar_returns_expected_keys_and_shapes() -> None:
    """MBAR solver should return a dict with the expected keys and array shapes."""

    from src.analyze.mbar import solve_mbar
    from src.config import MBARConfig

    xi_timeseries_list, window_centers, spring_constants, temperature_k, _, _ = (
        _make_harmonic_wham_dataset(samples_per_window=2000)
    )
    config = MBARConfig(n_pmf_bins=80, n_bootstrap=10)
    result = solve_mbar(
        xi_timeseries_list, window_centers, spring_constants, temperature_k, config,
    )

    assert "xi_bins" in result
    assert "pmf_kj_mol" in result
    assert "pmf_kcal_mol" in result
    assert "pmf_uncertainty_kj_mol" in result
    assert "free_energies_f_k" in result
    assert result["xi_bins"].shape == (config.n_pmf_bins,)
    assert result["pmf_kj_mol"].shape == (config.n_pmf_bins,)
    assert result["converged"] is True


# ---------- L-11 Step 3: Cross-validation tests ----------


def test_mbar_recovers_harmonic_pmf() -> None:
    """MBAR should recover a known harmonic PMF within 0.5 kJ/mol RMSE on the populated region."""

    from src.analyze.mbar import solve_mbar
    from src.config import MBARConfig

    (
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        center,
        spring_kj,
    ) = _make_harmonic_wham_dataset(samples_per_window=4000)

    config = MBARConfig(n_pmf_bins=80)
    result = solve_mbar(
        xi_timeseries_list, window_centers, spring_constants, temperature_k, config,
    )

    xi = result["xi_bins"]
    analytic = 0.5 * spring_kj * (xi - center) ** 2
    mask = np.isfinite(result["pmf_kj_mol"]) & (np.abs(xi) <= 0.5)
    shifted_mbar = result["pmf_kj_mol"][mask] - np.min(result["pmf_kj_mol"][mask])
    shifted_analytic = analytic[mask] - np.min(analytic[mask])
    rmse = float(np.sqrt(np.mean((shifted_mbar - shifted_analytic) ** 2)))

    assert rmse < 0.5, f"MBAR RMSE {rmse:.3f} kJ/mol exceeds 0.5 kJ/mol threshold"


def test_mbar_agrees_with_wham_on_harmonic_pmf() -> None:
    """MBAR and WHAM should produce PMFs within 0.75 kJ/mol of each other on synthetic data."""

    from src.analyze.mbar import solve_mbar
    from src.analyze.wham import solve_wham
    from src.config import MBARConfig, WHAMConfig

    (
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        _,
        _,
    ) = _make_harmonic_wham_dataset(samples_per_window=4000)

    n_bins = 80
    wham_result = solve_wham(
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        WHAMConfig(
            histogram_bins=n_bins,
            tolerance=1e-6,
            max_iterations=25_000,
            n_bootstrap=10,
        ),
    )
    mbar_result = solve_mbar(
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        MBARConfig(n_pmf_bins=n_bins),
    )

    wham_finite = np.isfinite(wham_result["pmf_kj_mol"])
    mbar_finite = np.isfinite(mbar_result["pmf_kj_mol"])
    shared = wham_finite & mbar_finite

    wham_pmf = wham_result["pmf_kj_mol"][shared] - np.min(wham_result["pmf_kj_mol"][shared])
    mbar_pmf = mbar_result["pmf_kj_mol"][shared] - np.min(mbar_result["pmf_kj_mol"][shared])

    max_diff = float(np.max(np.abs(wham_pmf - mbar_pmf)))
    assert max_diff < 0.75, (
        f"Max WHAM-MBAR difference {max_diff:.3f} kJ/mol exceeds 0.75 kJ/mol"
    )


def test_mbar_uncertainties_bounded_by_wham() -> None:
    """MBAR should produce uncertainties no larger than WHAM bootstrap uncertainties."""

    from src.analyze.mbar import bootstrap_mbar_uncertainty
    from src.analyze.wham import bootstrap_pmf_uncertainty
    from src.config import MBARConfig, WHAMConfig

    (
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        _,
        _,
    ) = _make_harmonic_wham_dataset(samples_per_window=3000)

    n_bins = 60
    n_boot = 50
    wham_unc = bootstrap_pmf_uncertainty(
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        WHAMConfig(
            histogram_bins=n_bins,
            n_bootstrap=n_boot,
            tolerance=1e-6,
            max_iterations=20_000,
        ),
    )
    mbar_unc = bootstrap_mbar_uncertainty(
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        MBARConfig(n_pmf_bins=n_bins, n_bootstrap=n_boot),
    )

    populated = (wham_unc["pmf_std"] > 0) & (mbar_unc["pmf_std"] > 0)
    ratio = mbar_unc["pmf_std"][populated] / wham_unc["pmf_std"][populated]
    median_ratio = float(np.median(ratio))
    assert median_ratio <= 1.05, (
        f"Median MBAR/WHAM uncertainty ratio {median_ratio:.3f} exceeds 1.05"
    )


def test_mbar_insensitive_to_bin_count() -> None:
    """MBAR PMF should be insensitive to the number of reporting bins."""

    from src.analyze.mbar import solve_mbar
    from src.config import MBARConfig

    (
        xi_timeseries_list,
        window_centers,
        spring_constants,
        temperature_k,
        _,
        _,
    ) = _make_harmonic_wham_dataset(samples_per_window=4000)

    results: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for n_bins in [50, 100, 200]:
        r = solve_mbar(
            xi_timeseries_list,
            window_centers,
            spring_constants,
            temperature_k,
            MBARConfig(n_pmf_bins=n_bins),
        )
        mask = np.isfinite(r["pmf_kj_mol"]) & (np.abs(r["xi_bins"]) <= 0.5)
        results[n_bins] = (
            r["xi_bins"][mask],
            r["pmf_kj_mol"][mask] - np.min(r["pmf_kj_mol"][mask]),
        )

    xi_50, pmf_50 = results[50]
    xi_200, pmf_200 = results[200]
    pmf_200_interp = np.interp(xi_50, xi_200, pmf_200)
    max_diff = float(np.max(np.abs(pmf_50 - pmf_200_interp)))
    # Threshold 0.5 kJ/mol accounts for finite-sample noise at the 200-bin
    # resolution (~180 samples/bin → σ ≈ 0.19 kJ/mol/bin).
    assert max_diff < 0.5, (
        f"MBAR bin-sensitivity: max diff = {max_diff:.4f} kJ/mol exceeds 0.5 kJ/mol"
    )
