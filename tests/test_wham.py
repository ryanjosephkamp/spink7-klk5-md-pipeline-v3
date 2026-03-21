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


# ---------- L-20 Step 3: assert → ValueError in wham.py ----------


def test_wham_raises_valueerror_for_2d_window_centers() -> None:
    """WHAM validation must raise ValueError for non-1D window centers."""

    with pytest.raises(ValueError, match="one-dimensional"):
        solve_wham(
            [np.asarray([1.0, 2.0])],
            window_centers=np.asarray([[1.5]]),  # 2D instead of 1D
            spring_constants=np.asarray([1000.0]),
            temperature_k=310.0,
            config=WHAMConfig(),
        )


# ---------- L-10 Step 5: Pre-WHAM coverage validation ----------


def test_wham_warns_on_coverage_gaps() -> None:
    """WHAM solver should warn when input histograms have zero-count bins."""
    import warnings

    rng = np.random.default_rng(42)
    # Two well-separated windows with a gap.
    xi_list = [
        rng.normal(loc=1.5, scale=0.02, size=200),
        rng.normal(loc=1.55, scale=0.02, size=200),
        rng.normal(loc=2.5, scale=0.02, size=200),
    ]
    centers = np.array([1.5, 1.55, 2.5])
    springs = np.full(3, 1000.0)
    config = WHAMConfig(histogram_bins=100, n_bootstrap=5, max_iterations=5000)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            solve_wham(xi_list, centers, springs, 310.0, config)
        except Exception:
            pass  # May raise PhysicalValidityError due to large gap.
        warning_msgs = [str(x.message) for x in w]
        assert any("zero counts" in msg for msg in warning_msgs), (
            f"Expected coverage warning, got: {warning_msgs}"
        )


# ---------- L-12 Step 1: Integrated autocorrelation time estimator ----------


def test_integrated_autocorrelation_time_recovers_known_tau():
    """The ACF estimator should recover tau_int for a synthetic AR(1) process with known correlation."""
    from src.analyze.wham import _integrated_autocorrelation_time

    rng = np.random.default_rng(42)
    n = 100_000
    phi = 0.95  # AR(1) coefficient
    # Analytic tau_int for AR(1): tau_int = (1 + phi) / (2*(1 - phi))
    analytic_tau = (1.0 + phi) / (2.0 * (1.0 - phi))  # = 19.5

    # Generate AR(1) process: x_t = phi * x_{t-1} + eps_t
    noise = rng.normal(0, 1, size=n)
    x = np.empty(n)
    x[0] = noise[0]
    for t in range(1, n):
        x[t] = phi * x[t - 1] + noise[t]

    tau_est = _integrated_autocorrelation_time(x)
    # Accept within 25% of analytic value (IPS has finite-sample bias)
    assert abs(tau_est - analytic_tau) / analytic_tau < 0.25, (
        f"Estimated tau_int={tau_est:.2f}, expected ~{analytic_tau:.1f}"
    )


def test_integrated_autocorrelation_time_uncorrelated_data():
    """For uncorrelated data, tau_int should be approximately 0.5."""
    from src.analyze.wham import _integrated_autocorrelation_time

    rng = np.random.default_rng(99)
    x = rng.normal(size=10_000)
    tau_est = _integrated_autocorrelation_time(x)
    assert tau_est < 2.0, f"tau_int={tau_est:.2f} for uncorrelated data; expected ~0.5"


# ---------- L-12 Step 2: Calibrated block sizes ----------


def test_bootstrap_calibrated_block_sizes_match_autocorrelation():
    """Block sizes should be approximately 2*tau_int for each window."""
    from src.analyze.wham import bootstrap_pmf_uncertainty
    from src.config import WHAMConfig

    # Create correlated synthetic data with known tau_int ~ 19.5
    rng = np.random.default_rng(7)
    phi = 0.95  # AR(1) with tau_int ~ 19.5
    n = 5000
    n_windows = 5

    xi_timeseries_list = []
    for _ in range(n_windows):
        noise = rng.normal(0, 0.01, size=n)
        x = np.empty(n)
        x[0] = noise[0]
        for t in range(1, n):
            x[t] = phi * x[t - 1] + noise[t]
        xi_timeseries_list.append(x)

    window_centers = np.linspace(-0.5, 0.5, n_windows)
    spring_constants = np.full(n_windows, 300.0)
    config = WHAMConfig(n_bootstrap=5, histogram_bins=30, tolerance=1e-5, max_iterations=10_000)

    result = bootstrap_pmf_uncertainty(
        xi_timeseries_list, window_centers, spring_constants, 310.0, config
    )

    assert "tau_int_per_window" in result
    assert "n_eff_per_window" in result
    assert "block_sizes_per_window" in result
    assert result["tau_int_per_window"].shape == (n_windows,)
    assert result["n_eff_per_window"].shape == (n_windows,)

    # block_size should be >= 2*tau_int for each window
    for bs, tau in zip(result["block_sizes_per_window"], result["tau_int_per_window"]):
        assert bs >= int(np.ceil(2.0 * tau)), (
            f"Block size {bs} < ceil(2*tau_int)={int(np.ceil(2.0 * tau))}"
        )


# ---------- L-12 Step 3: Autocorrelation-aware blocking validation ----------


def test_autocorrelation_aware_bootstrap_correct_uncertainty():
    """Bootstrap with block size = 2*tau_int should produce uncertainty within 60% of analytic value."""
    from src.analyze.wham import _bootstrap_resample_window, _integrated_autocorrelation_time

    rng_gen = np.random.default_rng(123)
    n = 20_000
    phi = 0.98  # AR(1) with tau_int = (1+phi)/(2*(1-phi)) = 49.5
    analytic_tau = (1.0 + phi) / (2.0 * (1.0 - phi))  # 49.5
    sigma_noise = 1.0
    # Analytic variance of AR(1) stationary distribution: sigma^2 / (1 - phi^2)
    analytic_var = sigma_noise**2 / (1.0 - phi**2)
    # Analytic variance of the sample mean: 2 * tau_int * var / N
    analytic_mean_var = 2.0 * analytic_tau * analytic_var / n

    # Generate AR(1) reference
    noise = rng_gen.normal(0, sigma_noise, size=n)
    x = np.empty(n)
    x[0] = noise[0]
    for t in range(1, n):
        x[t] = phi * x[t - 1] + noise[t]

    tau_est = _integrated_autocorrelation_time(x)
    correct_block = max(1, int(np.ceil(2.0 * tau_est)))

    # Bootstrap with correct block size
    n_boot = 500
    boot_rng = np.random.default_rng(456)
    means_correct = np.array([
        np.mean(_bootstrap_resample_window(x, boot_rng, block_size=correct_block))
        for _ in range(n_boot)
    ])
    bootstrap_var_correct = np.var(means_correct, ddof=1)

    ratio = bootstrap_var_correct / analytic_mean_var
    assert 0.4 < ratio < 2.0, (
        f"Correct-block bootstrap var ratio {ratio:.2f} outside [0.4, 2.0]"
    )


def test_naive_block_underestimates_uncertainty():
    """Bootstrap with sqrt(N) block (naive) should underestimate uncertainty for highly correlated data."""
    from src.analyze.wham import _bootstrap_resample_window, _integrated_autocorrelation_time

    rng_gen = np.random.default_rng(789)
    n = 20_000
    # phi = 0.999 gives tau_int ~ (1.999)/(0.002) = 999.5
    # sqrt(20000) ~ 141, much smaller than 2*tau_int ~ 1999
    phi = 0.999
    sigma_noise = 1.0

    noise = rng_gen.normal(0, sigma_noise, size=n)
    x = np.empty(n)
    x[0] = noise[0]
    for t in range(1, n):
        x[t] = phi * x[t - 1] + noise[t]

    tau_est = _integrated_autocorrelation_time(x)
    correct_block = max(1, int(np.ceil(2.0 * tau_est)))
    naive_block = max(1, int(round(np.sqrt(n))))  # ~ 141

    n_boot = 500
    boot_rng_naive = np.random.default_rng(111)
    boot_rng_correct = np.random.default_rng(222)

    means_naive = np.array([
        np.mean(_bootstrap_resample_window(x, boot_rng_naive, block_size=naive_block))
        for _ in range(n_boot)
    ])
    means_correct = np.array([
        np.mean(_bootstrap_resample_window(x, boot_rng_correct, block_size=correct_block))
        for _ in range(n_boot)
    ])

    var_naive = np.var(means_naive, ddof=1)
    var_correct = np.var(means_correct, ddof=1)

    # Naive should underestimate: var_naive < var_correct
    ratio = var_naive / var_correct
    assert ratio < 0.70, (
        f"Naive/correct variance ratio {ratio:.2f}; expected < 0.70 demonstrating underestimation"
    )