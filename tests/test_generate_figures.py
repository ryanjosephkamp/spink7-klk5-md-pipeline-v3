"""Tests for the automated figure generation pipeline."""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest


EXPECTED_FIGURES = [
    "fig3_protein_complex.png",
    "fig4_pmf_profile.png",
    "fig5_simulation_timeseries.png",
]


@pytest.mark.slow
def test_generate_figures_produces_all_outputs(tmp_path: Path) -> None:
    """Running generate_figures.py produces all expected PNG files."""

    import scripts.generate_figures as gf

    original_dir = gf.FIGURES_DIR
    gf.FIGURES_DIR = tmp_path
    try:
        gf.generate_figure_4(formats=("png",))
        gf.generate_figure_5(formats=("png",))
    finally:
        gf.FIGURES_DIR = original_dir
    plt.close("all")
    for name in ["fig4_pmf_profile.png", "fig5_simulation_timeseries.png"]:
        fig_path = tmp_path / name
        assert fig_path.exists(), f"Missing figure: {name}"
        assert fig_path.stat().st_size > 0, f"Empty figure: {name}"


def test_png_magic_bytes(tmp_path: Path) -> None:
    """Verify that generated PNG files have correct PNG magic bytes."""

    PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    out = tmp_path / "test.png"
    fig.savefig(out)
    plt.close(fig)
    with open(out, "rb") as f:
        header = f.read(8)
    assert header == PNG_MAGIC, "Invalid PNG header"


def test_generate_figure_4_synthetic_fallback(tmp_path: Path) -> None:
    """generate_figure_4 falls back to synthetic data when WHAM output is missing."""

    import scripts.generate_figures as gf

    original_dir = gf.FIGURES_DIR
    gf.FIGURES_DIR = tmp_path
    try:
        saved = gf.generate_figure_4(formats=("png",))
    finally:
        gf.FIGURES_DIR = original_dir
    plt.close("all")
    assert len(saved) == 1
    assert saved[0].exists()
    assert saved[0].stat().st_size > 0


def test_generate_figure_5_synthetic_fallback(tmp_path: Path) -> None:
    """generate_figure_5 falls back to synthetic data when production data is missing."""

    import scripts.generate_figures as gf

    original_dir = gf.FIGURES_DIR
    gf.FIGURES_DIR = tmp_path
    try:
        saved = gf.generate_figure_5(formats=("png",))
    finally:
        gf.FIGURES_DIR = original_dir
    plt.close("all")
    assert len(saved) == 1
    assert saved[0].exists()
    assert saved[0].stat().st_size > 0
