"""Tests for automated equilibration detection."""

from __future__ import annotations

import numpy as np
import pytest

from src.analyze.equilibration import detect_equilibration, equilibrated_subseries


def test_detect_equilibration_pure_stationary():
    """A purely stationary timeseries should have t0 near zero."""
    rng = np.random.default_rng(99)
    stationary = rng.normal(5.0, 0.3, 500)

    result = detect_equilibration(stationary)

    assert result["t0"] < 50, f"Expected t0 near 0, got {result['t0']}"
    assert result["equilibrated"] is True


def test_detect_equilibration_drift_then_stationary():
    """Detection should find the transition point in a drift+stationary signal."""
    rng = np.random.default_rng(42)
    drift = np.linspace(10.0, 0.0, 100)
    stationary = rng.normal(0.0, 0.5, 900)
    timeseries = np.concatenate([drift, stationary])

    result = detect_equilibration(timeseries)

    assert abs(result["t0"] - 100) < 100, (
        f"Expected t0 near 100, got {result['t0']}"
    )
    assert result["equilibrated"] is True
    assert result["n_eff"] > 100


def test_detect_equilibration_warns_when_mostly_transient():
    """A warning should be emitted when equilibration consumes >80% of samples."""
    rng = np.random.default_rng(7)
    drift = np.linspace(100.0, 0.0, 900)
    stationary = rng.normal(0.0, 0.2, 100)
    timeseries = np.concatenate([drift, stationary])

    with pytest.warns(UserWarning, match="threshold"):
        result = detect_equilibration(timeseries)

    assert result["equilibrated"] is False


def test_equilibrated_subseries_returns_correct_slice():
    """equilibrated_subseries should return values[t0:]."""
    values = np.arange(100, dtype=float)
    sub = equilibrated_subseries(values, t0=30)
    assert sub.size == 70
    assert sub[0] == 30.0


def test_equilibrated_subseries_rejects_invalid_t0():
    """t0 outside [0, N) should raise ValueError."""
    values = np.arange(10, dtype=float)
    with pytest.raises(ValueError, match="t0 must be in"):
        equilibrated_subseries(values, t0=10)


def test_detect_equilibration_short_series():
    """Very short timeseries should return t0=0 without error."""
    short = np.array([1.0, 2.0, 3.0])
    result = detect_equilibration(short)
    assert result["t0"] == 0
    assert result["equilibrated"] is True
