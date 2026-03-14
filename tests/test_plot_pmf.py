"""Tests for PMF plotting utilities."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from src.visualization.plot_pmf import plot_pmf


def _synthetic_pmf_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate compact synthetic PMF inputs for plotting tests."""

    xi_bins = np.linspace(1.5, 4.0, 25, dtype=float)
    pmf = 2.5 * (xi_bins - 2.1) ** 2 - 1.8
    pmf_std = np.full_like(xi_bins, 0.25)
    return xi_bins, pmf, pmf_std


def test_plot_pmf_returns_figure_without_uncertainty() -> None:
    """Chunk 22 gate: PMF plotting should return a matplotlib Figure without error."""

    xi_bins, pmf, _ = _synthetic_pmf_inputs()

    figure = plot_pmf(xi_bins, pmf)

    assert isinstance(figure, Figure)
    axis = figure.axes[0]
    assert axis.get_xlabel() == "COM Distance (nm)"
    assert axis.get_ylabel() == "PMF (kcal/mol)"
    plt.close(figure)


def test_plot_pmf_returns_figure_with_uncertainty_and_saves_output(tmp_path: Path) -> None:
    """PMF plotting should support an uncertainty band and save to disk when requested."""

    xi_bins, pmf, pmf_std = _synthetic_pmf_inputs()
    output_path = tmp_path / "figures" / "pmf_profile.png"

    figure = plot_pmf(xi_bins, pmf, pmf_std_kcal_mol=pmf_std, output_path=output_path)

    assert isinstance(figure, Figure)
    assert output_path.exists()
    assert output_path.stat().st_size > 0
    plt.close(figure)


def test_plot_pmf_rejects_length_mismatch() -> None:
    """Public API should reject PMF arrays whose lengths do not match the xi grid."""

    xi_bins, pmf, _ = _synthetic_pmf_inputs()

    try:
        plot_pmf(xi_bins, pmf[:-1])
    except ValueError as exc:
        assert "same length" in str(exc)
    else:
        raise AssertionError("plot_pmf should reject xi/PMF length mismatches")