"""Automated equilibration detection using the Chodera method."""

from __future__ import annotations

import logging
import warnings

import numpy as np

from src.analyze.convergence import autocorrelation_time

logger = logging.getLogger(__name__)


def _validate_timeseries(timeseries: np.ndarray) -> np.ndarray:
    """Validate a one-dimensional finite timeseries."""
    values = np.asarray(timeseries, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("timeseries must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("timeseries must contain only finite values")
    return values


def detect_equilibration(
    timeseries: np.ndarray,
    *,
    warn_fraction: float = 0.8,
) -> dict[str, float | int | bool]:
    """Detect the equilibration time using maximum effective samples.

    Implements the Chodera method [1]: for each candidate discard index
    t0, compute the statistical inefficiency g(t0) of the remaining
    sub-series and select the t0 that maximizes N_eff = (N - t0) / g(t0).

    Args:
        timeseries: Observable values sampled at uniform intervals.
            Shape: [N_samples].
        warn_fraction: Emit a warning if t_eq exceeds this fraction
            of the total timeseries length.

    Returns:
        dict with keys:
            - t0: optimal discard index (int)
            - n_eff: number of effective samples after discarding (float)
            - g: statistical inefficiency of the equilibrated segment (float)
            - equilibrated: True if t0 < warn_fraction * N (bool)

    References:
        [1] J. D. Chodera, J. Chem. Theory Comput. 12, 1799-1805 (2016).
    """
    values = _validate_timeseries(timeseries)
    n_total = values.size

    if n_total < 4:
        return {"t0": 0, "n_eff": float(n_total), "g": 1.0, "equilibrated": True}

    best_t0 = 0
    best_n_eff = 0.0
    best_g = 1.0

    # Minimum sub-series length for meaningful autocorrelation estimation
    min_subseries = max(4, n_total // 20)

    for t0 in range(n_total - min_subseries):
        subseries = values[t0:]
        g = 2.0 * autocorrelation_time(subseries)
        g = max(g, 1.0)
        n_eff = (n_total - t0) / g
        if n_eff > best_n_eff:
            best_n_eff = n_eff
            best_t0 = t0
            best_g = g

    equilibrated = best_t0 < warn_fraction * n_total

    if not equilibrated:
        warnings.warn(
            f"Equilibration consumes {best_t0}/{n_total} samples "
            f"({100.0 * best_t0 / n_total:.1f}%), exceeding the "
            f"{100.0 * warn_fraction:.0f}% threshold. Consider extending "
            f"the simulation.",
            stacklevel=2,
        )

    logger.info(
        "Equilibration detection: t0=%d, n_eff=%.1f, g=%.2f, equilibrated=%s",
        best_t0, best_n_eff, best_g, equilibrated,
    )

    return {
        "t0": int(best_t0),
        "n_eff": float(best_n_eff),
        "g": float(best_g),
        "equilibrated": bool(equilibrated),
    }


def equilibrated_subseries(
    timeseries: np.ndarray,
    t0: int,
) -> np.ndarray:
    """Return the post-equilibration portion of a timeseries.

    Args:
        timeseries: Full observable timeseries. Shape: [N_samples].
        t0: Discard index from detect_equilibration().

    Returns:
        np.ndarray: Equilibrated sub-series. Shape: [N_samples - t0].
    """
    values = _validate_timeseries(timeseries)
    discard = int(t0)
    if discard < 0 or discard >= values.size:
        raise ValueError("t0 must be in [0, N_samples)")
    return values[discard:]
