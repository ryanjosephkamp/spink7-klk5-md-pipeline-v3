"""Tests for statistical convergence diagnostics."""

from __future__ import annotations

import numpy as np
import pytest

from src.analyze.convergence import autocorrelation_time, block_average, effective_sample_size


def _make_ar1_series(phi: float, n_samples: int, seed: int) -> np.ndarray:
    """Generate an AR(1) process with unit-variance driving noise."""

    rng = np.random.default_rng(seed)
    noise = rng.normal(loc=0.0, scale=1.0, size=n_samples)
    series = np.empty(n_samples, dtype=float)
    series[0] = noise[0]
    for index in range(1, n_samples):
        series[index] = phi * series[index - 1] + noise[index]
    return series


def test_block_average_matches_manual_block_statistics() -> None:
    """Chunk 20 gate: block averaging should reproduce direct block-mean statistics."""

    timeseries = np.arange(1.0, 13.0, dtype=float)

    result = block_average(timeseries, n_blocks=4)

    block_means = np.asarray([2.0, 5.0, 8.0, 11.0], dtype=float)
    assert result["mean"] == pytest.approx(np.mean(block_means), abs=1e-12)
    assert result["sem"] == pytest.approx(np.std(block_means, ddof=1) / np.sqrt(4), abs=1e-12)
    assert result["block_size"] == pytest.approx(3.0, abs=1e-12)


def test_autocorrelation_time_and_effective_sample_size_distinguish_correlated_series() -> None:
    """A strongly correlated AR(1) series should have larger τ and smaller N_eff than white noise."""

    white_noise = _make_ar1_series(phi=0.0, n_samples=20_000, seed=7)
    correlated = _make_ar1_series(phi=0.9, n_samples=20_000, seed=11)

    tau_white = autocorrelation_time(white_noise)
    tau_correlated = autocorrelation_time(correlated)
    n_eff_white = effective_sample_size(white_noise)
    n_eff_correlated = effective_sample_size(correlated)

    assert tau_white < 1.5
    assert tau_correlated > 5.0
    assert tau_correlated > tau_white
    assert n_eff_white > 0.5 * white_noise.size
    assert n_eff_correlated < 0.2 * correlated.size
    assert n_eff_correlated < n_eff_white


def test_convergence_public_api_rejects_invalid_inputs() -> None:
    """Public API should reject empty timeseries and invalid block counts."""

    with pytest.raises(ValueError, match="non-empty"):
        autocorrelation_time(np.asarray([], dtype=float))

    with pytest.raises(ValueError, match="greater than 1"):
        block_average(np.asarray([1.0, 2.0, 3.0], dtype=float), n_blocks=1)

    with pytest.raises(ValueError, match="must not exceed"):
        block_average(np.asarray([1.0, 2.0], dtype=float), n_blocks=3)