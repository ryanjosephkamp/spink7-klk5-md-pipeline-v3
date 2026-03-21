"""Statistical convergence diagnostics for simulation timeseries."""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)


def _validate_timeseries(timeseries: np.ndarray) -> np.ndarray:
    """Validate a one-dimensional finite timeseries at the public API boundary."""

    values = np.asarray(timeseries, dtype=float)
    if values.ndim != 1:
        raise ValueError(f"timeseries must be one-dimensional, got shape {values.shape}")
    if values.size == 0:
        raise ValueError("timeseries must be non-empty")
    if not np.all(np.isfinite(values)):
        raise ValueError("timeseries must contain only finite values")
    return values


def block_average(
    timeseries: np.ndarray,
    n_blocks: int = 10,
) -> dict[str, float]:
    """Estimate the mean and standard error by block averaging a timeseries."""

    values = _validate_timeseries(timeseries)
    block_count = int(n_blocks)
    if block_count <= 1:
        raise ValueError("n_blocks must be greater than 1")

    block_size = values.size // block_count
    if block_size == 0:
        raise ValueError("n_blocks must not exceed the number of samples")

    trimmed = values[: block_count * block_size]
    reshaped = trimmed.reshape(block_count, block_size)
    block_means = np.mean(reshaped, axis=1, dtype=float)
    mean = float(np.mean(block_means, dtype=float))
    sem = float(np.std(block_means, ddof=1, dtype=float) / np.sqrt(block_count))

    return {
        "mean": mean,
        "sem": sem,
        "block_size": float(block_size),
    }


def autocorrelation_time(
    timeseries: np.ndarray,
) -> float:
    """Estimate the integrated autocorrelation time in samples using FFT autocovariance."""

    values = _validate_timeseries(timeseries)
    n_samples = values.size
    if n_samples == 1:
        return 0.5

    centered = values - np.mean(values, dtype=float)
    variance = float(np.var(centered, ddof=0, dtype=float))
    if variance <= np.finfo(float).eps:
        return 0.5

    # FFT-based autocovariance: pad to the next power of two, compute
    # the power spectrum, and inverse-transform to obtain the biased
    # autocovariance function normalized by overlap count per lag.
    fft_size = 1 << (2 * n_samples - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=fft_size)
    autocovariance = np.fft.irfft(spectrum * np.conjugate(spectrum), n=fft_size)[:n_samples]
    normalization = np.arange(n_samples, 0, -1, dtype=float)
    autocovariance /= normalization
    autocorrelation = autocovariance / autocovariance[0]

    # Integrate the autocorrelation up to the first non-positive lag
    # to estimate the integrated autocorrelation time.
    positive_lags = autocorrelation[1:]
    negative_mask = positive_lags <= 0.0
    cutoff = int(np.argmax(negative_mask)) if np.any(negative_mask) else positive_lags.size
    truncated = positive_lags[:cutoff]
    tau_int = 0.5 + float(np.sum(truncated, dtype=float))
    return max(0.5, tau_int)


def effective_sample_size(
    timeseries: np.ndarray,
) -> int:
    """Estimate the effective sample size implied by the integrated autocorrelation time."""

    values = _validate_timeseries(timeseries)
    tau_int = autocorrelation_time(values)
    n_eff = int(max(1.0, np.floor(values.size / (2.0 * tau_int))))
    return min(values.size, n_eff)


def compare_fes_profiles(
    fes_a_grid: np.ndarray,
    fes_a_values: np.ndarray,
    fes_b_grid: np.ndarray,
    fes_b_values: np.ndarray,
    method_a_name: str = "Metadynamics",
    method_b_name: str = "WHAM",
) -> dict[str, float]:
    """Compare two free energy surfaces on potentially different grids.

    Both profiles are interpolated onto a common grid (the finer of the two),
    aligned by setting G(xi_max) = 0, and compared via max absolute deviation,
    RMSD, and Pearson correlation coefficient.

    Parameters
    ----------
    fes_a_grid, fes_b_grid : np.ndarray
        1D arrays of grid points for each FES.
    fes_a_values, fes_b_values : np.ndarray
        1D arrays of free energy values (kJ/mol) for each FES.
    method_a_name, method_b_name : str
        Labels for logging.

    Returns
    -------
    dict with keys: max_absolute_deviation_kj_mol, rmsd_kj_mol, pearson_r
    """
    grid_a = np.asarray(fes_a_grid, dtype=float)
    grid_b = np.asarray(fes_b_grid, dtype=float)
    vals_a = np.asarray(fes_a_values, dtype=float)
    vals_b = np.asarray(fes_b_values, dtype=float)

    # Determine common grid range
    lo = max(grid_a[0], grid_b[0])
    hi = min(grid_a[-1], grid_b[-1])
    n_points = max(len(grid_a), len(grid_b))
    common_grid = np.linspace(lo, hi, n_points)

    # Interpolate both onto common grid
    interp_a = np.interp(common_grid, grid_a, vals_a)
    interp_b = np.interp(common_grid, grid_b, vals_b)

    # Align: set G(xi_max) = 0
    interp_a = interp_a - interp_a[-1]
    interp_b = interp_b - interp_b[-1]

    # Compute metrics
    diff = interp_a - interp_b
    max_abs_dev = float(np.max(np.abs(diff)))
    rmsd = float(np.sqrt(np.mean(diff ** 2)))

    # Pearson correlation
    std_a = float(np.std(interp_a))
    std_b = float(np.std(interp_b))
    if std_a > 0 and std_b > 0:
        pearson_r = float(np.corrcoef(interp_a, interp_b)[0, 1])
    else:
        pearson_r = 1.0 if np.allclose(interp_a, interp_b) else 0.0

    logger.info(
        "FES comparison (%s vs %s): max_dev=%.3f kJ/mol, RMSD=%.3f kJ/mol, r=%.4f",
        method_a_name, method_b_name, max_abs_dev, rmsd, pearson_r,
    )

    return {
        "max_absolute_deviation_kj_mol": max_abs_dev,
        "rmsd_kj_mol": rmsd,
        "pearson_r": pearson_r,
    }