"""Timeseries plotting utilities for simulation observables."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from src.visualization._metadata import metadata_for_format


def _validate_series_pair(time_values: np.ndarray, observable: np.ndarray, *, time_name: str, observable_name: str) -> tuple[np.ndarray, np.ndarray]:
    """Validate a time axis and matching observable series at the public API boundary."""

    time_array = np.asarray(time_values, dtype=float)
    observable_array = np.asarray(observable, dtype=float)
    assert time_array.ndim == 1, f"{time_name} must have shape [N]"
    assert observable_array.ndim == 1, f"{observable_name} must have shape [N]"
    if time_array.ndim != 1 or observable_array.ndim != 1 or time_array.size == 0 or observable_array.size == 0:
        raise ValueError(f"{time_name} and {observable_name} must be non-empty one-dimensional arrays")
    if time_array.size != observable_array.size:
        raise ValueError(f"{time_name} and {observable_name} must have the same length")
    if not np.all(np.isfinite(time_array)) or not np.all(np.isfinite(observable_array)):
        raise ValueError(f"{time_name} and {observable_name} must contain only finite values")
    return time_array, observable_array


def _validate_output_path(
    output_path: Path | list[Path] | None,
) -> list[Path] | None:
    """Validate optional figure output paths.

    Accepts a single ``Path``, a list of ``Path`` objects, or ``None``.
    Returns a list of validated paths or ``None``.
    """

    if output_path is None:
        return None
    _SUPPORTED_EXT = {".png", ".svg", ".pdf", ".eps"}
    paths = output_path if isinstance(output_path, list) else [Path(output_path)]
    for p in paths:
        p = Path(p)
        if p.suffix.lower() not in _SUPPORTED_EXT:
            raise ValueError(
                f"output_path must use a supported image extension: {', '.join(sorted(_SUPPORTED_EXT))}"
            )
    return [Path(p) for p in paths]


def _finalize_figure(figure: Figure, output_path: list[Path] | None) -> Figure:
    """Persist a figure when requested and return it."""

    if output_path is not None:
        for _out in output_path:
            _out.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(_out, dpi=300, metadata=metadata_for_format(_out))
    return figure


def plot_energy_timeseries(
    time_ps: np.ndarray,
    potential_kj_mol: np.ndarray,
    kinetic_kj_mol: np.ndarray,
    output_path: Path | list[Path] | None = None,
) -> Figure:
    """Plot potential and kinetic energy as a function of time.

    Invariants: None.

    Args:
        time_ps: Simulation time in ps. Shape: [N_frames].
        potential_kj_mol: Potential energy in kJ/mol. Shape: [N_frames].
        kinetic_kj_mol: Kinetic energy in kJ/mol. Shape: [N_frames].
        output_path: Optional figure-save path.

    Returns:
        Figure: A matplotlib figure containing the energy timeseries.
    """

    time_array, potential_array = _validate_series_pair(
        time_ps,
        potential_kj_mol,
        time_name="time_ps",
        observable_name="potential_kj_mol",
    )
    _, kinetic_array = _validate_series_pair(
        time_array,
        kinetic_kj_mol,
        time_name="time_ps",
        observable_name="kinetic_kj_mol",
    )
    validated_output_path = _validate_output_path(output_path)

    figure, axis = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    axis.plot(time_array, potential_array, color="#0b3c5d", linewidth=2.0, label="Potential")
    axis.plot(time_array, kinetic_array, color="#d94f04", linewidth=2.0, label="Kinetic")
    axis.set_xlabel("Time (ps)")
    axis.set_ylabel("Energy (kJ/mol)")
    axis.set_title("Energy Timeseries")
    axis.grid(True, alpha=0.2)
    axis.legend(frameon=False)
    return _finalize_figure(figure, validated_output_path)


def plot_temperature_timeseries(
    time_ps: np.ndarray,
    temperature_k: np.ndarray,
    output_path: Path | list[Path] | None = None,
) -> Figure:
    """Plot temperature as a function of time.

    Invariants: None.

    Args:
        time_ps: Simulation time in ps. Shape: [N_frames].
        temperature_k: Temperature in Kelvin. Shape: [N_frames].
        output_path: Optional figure-save path.

    Returns:
        Figure: A matplotlib figure containing the temperature timeseries.
    """

    time_array, temperature_array = _validate_series_pair(
        time_ps,
        temperature_k,
        time_name="time_ps",
        observable_name="temperature_k",
    )
    validated_output_path = _validate_output_path(output_path)

    figure, axis = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    axis.plot(time_array, temperature_array, color="#328cc1", linewidth=2.0)
    axis.set_xlabel("Time (ps)")
    axis.set_ylabel("Temperature (K)")
    axis.set_title("Temperature Timeseries")
    axis.grid(True, alpha=0.2)
    return _finalize_figure(figure, validated_output_path)


def plot_rmsd_timeseries(
    time_ns: np.ndarray,
    rmsd_nm: np.ndarray,
    output_path: Path | list[Path] | None = None,
) -> Figure:
    """Plot RMSD as a function of time.

    Invariants: None.

    Args:
        time_ns: Simulation time in ns. Shape: [N_frames].
        rmsd_nm: RMSD values in nm. Shape: [N_frames].
        output_path: Optional figure-save path.

    Returns:
        Figure: A matplotlib figure containing the RMSD timeseries.
    """

    time_array, rmsd_array = _validate_series_pair(
        time_ns,
        rmsd_nm,
        time_name="time_ns",
        observable_name="rmsd_nm",
    )
    validated_output_path = _validate_output_path(output_path)

    figure, axis = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    axis.plot(time_array, rmsd_array, color="#6c4f77", linewidth=2.0)
    axis.set_xlabel("Time (ns)")
    axis.set_ylabel("RMSD (nm)")
    axis.set_title("RMSD Timeseries")
    axis.grid(True, alpha=0.2)
    return _finalize_figure(figure, validated_output_path)