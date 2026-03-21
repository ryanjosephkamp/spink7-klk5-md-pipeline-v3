"""Publication-style PMF plotting utilities."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from src.visualization._metadata import metadata_for_format


def _validate_bins_and_pmf(xi_bins_nm: np.ndarray, pmf_kcal_mol: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Validate the core PMF plotting arrays at the public API boundary."""

    xi_bins = np.asarray(xi_bins_nm, dtype=float)
    pmf = np.asarray(pmf_kcal_mol, dtype=float)
    assert xi_bins.ndim == 1, "xi_bins_nm must have shape [N_bins]"
    assert pmf.ndim == 1, "pmf_kcal_mol must have shape [N_bins]"
    if xi_bins.ndim != 1 or pmf.ndim != 1 or xi_bins.size == 0 or pmf.size == 0:
        raise ValueError("xi_bins_nm and pmf_kcal_mol must be non-empty one-dimensional arrays")
    if xi_bins.size != pmf.size:
        raise ValueError("xi_bins_nm and pmf_kcal_mol must have the same length")
    if not np.all(np.isfinite(xi_bins)) or not np.all(np.isfinite(pmf)):
        raise ValueError("xi_bins_nm and pmf_kcal_mol must contain only finite values")
    return xi_bins, pmf


def _validate_uncertainty(pmf_std_kcal_mol: np.ndarray | None, n_bins: int) -> np.ndarray | None:
    """Validate optional PMF uncertainty arrays."""

    if pmf_std_kcal_mol is None:
        return None

    pmf_std = np.asarray(pmf_std_kcal_mol, dtype=float)
    assert pmf_std.ndim == 1, "pmf_std_kcal_mol must have shape [N_bins]"
    if pmf_std.ndim != 1 or pmf_std.size != n_bins:
        raise ValueError("pmf_std_kcal_mol must be a one-dimensional array with one value per PMF bin")
    if not np.all(np.isfinite(pmf_std)):
        raise ValueError("pmf_std_kcal_mol must contain only finite values")
    if np.any(pmf_std < 0.0):
        raise ValueError("pmf_std_kcal_mol must be non-negative")
    return pmf_std


def _validate_output_path(
    output_path: Path | list[Path] | None,
) -> list[Path] | None:
    """Validate the optional figure output path(s).

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


def plot_pmf(
    xi_bins_nm: np.ndarray,
    pmf_kcal_mol: np.ndarray,
    pmf_std_kcal_mol: np.ndarray | None = None,
    output_path: Path | list[Path] | None = None,
) -> Figure:
    """Plot a PMF profile with optional uncertainty and bound-state annotation.

    Invariants: None.

    Args:
        xi_bins_nm: Reaction-coordinate bin centers in nm. Shape: [N_bins].
        pmf_kcal_mol: PMF values in kcal/mol. Shape: [N_bins].
        pmf_std_kcal_mol: Optional PMF standard deviations in kcal/mol. Shape: [N_bins].
        output_path: Optional image-save path.

    Returns:
        Figure: A matplotlib figure for the PMF profile.
    """

    xi_bins, pmf = _validate_bins_and_pmf(xi_bins_nm, pmf_kcal_mol)
    pmf_std = _validate_uncertainty(pmf_std_kcal_mol, xi_bins.size)
    validated_output_path = _validate_output_path(output_path)

    figure, axis = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    axis.plot(xi_bins, pmf, color="#0b3c5d", linewidth=2.2)

    if pmf_std is not None:
        axis.fill_between(
            xi_bins,
            pmf - pmf_std,
            pmf + pmf_std,
            color="#328cc1",
            alpha=0.25,
            linewidth=0.0,
        )

    minimum_index = int(np.argmin(pmf))
    minimum_x = float(xi_bins[minimum_index])
    minimum_y = float(pmf[minimum_index])
    axis.scatter([minimum_x], [minimum_y], color="#d94f04", s=35, zorder=3)
    axis.annotate(
        f"ΔG_bind ≈ {minimum_y:.2f} kcal/mol",
        xy=(minimum_x, minimum_y),
        xytext=(12, 12),
        textcoords="offset points",
        fontsize=10,
        color="#222222",
        arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.0},
    )

    axis.set_xlabel("COM Distance (nm)")
    axis.set_ylabel("PMF (kcal/mol)")
    axis.set_title("Potential of Mean Force")
    axis.grid(True, alpha=0.2)

    if validated_output_path is not None:
        for _out in validated_output_path:
            _out.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(_out, dpi=300, metadata=metadata_for_format(_out))

    return figure