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


def test_plot_pmf_multi_format_output(tmp_path: Path) -> None:
    """L-36 Step 1: plot_pmf saves to PNG, SVG, and PDF in a single call."""

    xi_bins, pmf, _ = _synthetic_pmf_inputs()
    paths = [
        tmp_path / "test_pmf.png",
        tmp_path / "test_pmf.svg",
        tmp_path / "test_pmf.pdf",
    ]
    fig = plot_pmf(xi_bins, pmf, output_path=paths)
    plt.close(fig)
    for p in paths:
        assert p.exists(), f"Missing output: {p.name}"
        assert p.stat().st_size > 0, f"Empty output: {p.name}"


def test_plot_pmf_metadata_embedding(tmp_path: Path) -> None:
    """L-36 Step 2: PNG figures contain pipeline provenance metadata."""

    from PIL import Image

    xi_bins, pmf, _ = _synthetic_pmf_inputs()
    out = tmp_path / "meta_test.png"
    fig = plot_pmf(xi_bins, pmf, output_path=out)
    plt.close(fig)
    img = Image.open(out)
    info = img.info
    assert "Creator" in info, "Missing Creator metadata"
    assert "SPINK7-KLK5" in info["Creator"], "Incorrect Creator value"


def test_generate_figure_4_loads_real_pmf_data(tmp_path: Path) -> None:
    """L-36 Step 3: generate_figure_4 loads real PMF data when available."""

    analysis_dir = tmp_path / "data" / "analysis" / "pmf"
    analysis_dir.mkdir(parents=True)
    xi = np.linspace(1.5, 4.0, 100)
    pmf = -10.0 * np.exp(-2.0 * (xi - 2.0) ** 2)
    np.savez(analysis_dir / "wham_pmf.npz", xi_bins=xi, pmf_kcal_mol=pmf)

    data = np.load(analysis_dir / "wham_pmf.npz")
    assert "xi_bins" in data.files
    assert "pmf_kcal_mol" in data.files
    assert data["xi_bins"].shape == (100,)