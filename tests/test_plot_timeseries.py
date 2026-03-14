"""Tests for timeseries plotting utilities."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from src.visualization.plot_timeseries import (
    plot_energy_timeseries,
    plot_rmsd_timeseries,
    plot_temperature_timeseries,
)


def _synthetic_timeseries_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Create compact synthetic timeseries arrays for plotting tests."""

    time_ps = np.linspace(0.0, 100.0, 51, dtype=float)
    potential_kj_mol = -850.0 + 5.0 * np.sin(time_ps / 20.0)
    kinetic_kj_mol = 320.0 + 2.5 * np.cos(time_ps / 18.0)
    temperature_k = 310.0 + 3.0 * np.sin(time_ps / 25.0)
    time_ns = time_ps / 1000.0
    rmsd_nm = 0.08 + 0.02 * np.sin(time_ps / 30.0)
    return time_ps, potential_kj_mol, kinetic_kj_mol, temperature_k, time_ns, rmsd_nm


def test_plot_energy_timeseries_returns_figure() -> None:
    """Chunk 23 gate: energy plotting should return a matplotlib Figure without crashing."""

    time_ps, potential_kj_mol, kinetic_kj_mol, _, _, _ = _synthetic_timeseries_inputs()

    figure = plot_energy_timeseries(time_ps, potential_kj_mol, kinetic_kj_mol)

    assert isinstance(figure, Figure)
    axis = figure.axes[0]
    assert axis.get_xlabel() == "Time (ps)"
    assert axis.get_ylabel() == "Energy (kJ/mol)"
    plt.close(figure)


def test_plot_temperature_timeseries_saves_output(tmp_path: Path) -> None:
    """Temperature plotting should save a figure when requested."""

    time_ps, _, _, temperature_k, _, _ = _synthetic_timeseries_inputs()
    output_path = tmp_path / "figures" / "temperature_timeseries.png"

    figure = plot_temperature_timeseries(time_ps, temperature_k, output_path=output_path)

    assert isinstance(figure, Figure)
    assert output_path.exists()
    assert output_path.stat().st_size > 0
    axis = figure.axes[0]
    assert axis.get_xlabel() == "Time (ps)"
    assert axis.get_ylabel() == "Temperature (K)"
    plt.close(figure)


def test_plot_rmsd_timeseries_returns_figure() -> None:
    """RMSD plotting should return a matplotlib Figure without crashing."""

    _, _, _, _, time_ns, rmsd_nm = _synthetic_timeseries_inputs()

    figure = plot_rmsd_timeseries(time_ns, rmsd_nm)

    assert isinstance(figure, Figure)
    axis = figure.axes[0]
    assert axis.get_xlabel() == "Time (ns)"
    assert axis.get_ylabel() == "RMSD (nm)"
    plt.close(figure)


def test_plot_energy_timeseries_rejects_length_mismatch() -> None:
    """Public API should reject mismatched time and observable array lengths."""

    time_ps, potential_kj_mol, kinetic_kj_mol, _, _, _ = _synthetic_timeseries_inputs()

    try:
        plot_energy_timeseries(time_ps, potential_kj_mol[:-1], kinetic_kj_mol)
    except ValueError as exc:
        assert "same length" in str(exc)
    else:
        raise AssertionError("plot_energy_timeseries should reject length mismatches")