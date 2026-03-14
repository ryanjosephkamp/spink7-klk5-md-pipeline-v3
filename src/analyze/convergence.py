"""Statistical convergence diagnostics for simulation timeseries."""

from __future__ import annotations

import numpy as np


def _validate_timeseries(timeseries: np.ndarray) -> np.ndarray:
    """Validate a one-dimensional finite timeseries at the public API boundary."""

    values = np.asarray(timeseries, dtype=float)
    assert values.ndim == 1, "timeseries must have shape [N]"
    if values.ndim != 1 or values.size == 0:
        raise ValueError("timeseries must be a non-empty one-dimensional array")
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