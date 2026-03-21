"""Tests for the cross-validation script (scripts/cross_validate.py)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from scripts.cross_validate import build_parser, cross_validate, main
from src.config import BOLTZMANN_KJ


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_jarzynski_npz(path: Path, delta_g_kj_mol: float) -> Path:
    """Create a minimal Jarzynski summary NPZ."""
    np.savez(
        path,
        delta_g_kj_mol=delta_g_kj_mol,
        delta_g_kcal_mol=delta_g_kj_mol / 4.184,
        mean_work_kj_mol=delta_g_kj_mol + 5.0,
    )
    return path


def _make_wham_npz(path: Path, dg_kj_mol: float) -> Path:
    """Create a minimal WHAM PMF NPZ with a known ΔG.

    The PMF has a minimum at ``dg_kj_mol`` below the dissociated plateau (0).
    """
    xi_bins = np.linspace(1.5, 4.0, 50)
    pmf_kj = np.zeros(50, dtype=float)
    # Set the minimum at the first bin to simulate a binding well.
    pmf_kj[0] = dg_kj_mol  # negative value = binding
    np.savez(path, xi_bins=xi_bins, pmf_kj_mol=pmf_kj, pmf_kcal_mol=pmf_kj / 4.184)
    return path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_cross_validate_passing(tmp_path: Path) -> None:
    """Cross-validation passes when ΔG estimates agree within 2 k_BT."""
    dg = -42.0  # kJ/mol
    jarz_npz = _make_jarzynski_npz(tmp_path / "jarz.npz", dg)
    wham_npz = _make_wham_npz(tmp_path / "wham.npz", dg)

    result = cross_validate(jarz_npz, wham_npz, temperature_k=310.0)

    assert result["passed"] == 1.0
    assert result["discrepancy_kj_mol"] < 2.0 * BOLTZMANN_KJ * 310.0


def test_cross_validate_failing(tmp_path: Path) -> None:
    """Cross-validation fails when ΔG estimates differ by more than 2 k_BT."""
    kbt = BOLTZMANN_KJ * 310.0
    jarz_npz = _make_jarzynski_npz(tmp_path / "jarz.npz", -42.0)
    wham_npz = _make_wham_npz(tmp_path / "wham.npz", -42.0 - 3.0 * kbt)

    result = cross_validate(jarz_npz, wham_npz, temperature_k=310.0)

    assert result["passed"] == 0.0
    assert result["discrepancy_kj_mol"] > 2.0 * kbt


def test_cross_validate_cli_parser() -> None:
    """CLI parser accepts expected arguments."""
    parser = build_parser()
    args = parser.parse_args([
        "--jarzynski-npz", "/tmp/jarz.npz",
        "--wham-npz", "/tmp/wham.npz",
        "--temperature-k", "310.0",
        "--threshold-kbt", "2.0",
    ])
    assert args.temperature_k == 310.0
    assert args.threshold_kbt == 2.0


def test_cross_validate_main_exit_code(tmp_path: Path) -> None:
    """CLI returns exit code 0 on pass and 1 on fail."""
    dg = -42.0
    jarz_npz = _make_jarzynski_npz(tmp_path / "jarz.npz", dg)
    wham_npz = _make_wham_npz(tmp_path / "wham.npz", dg)

    rc = main([
        "--jarzynski-npz", str(jarz_npz),
        "--wham-npz", str(wham_npz),
        "--temperature-k", "310.0",
        "--output", str(tmp_path / "result.npz"),
    ])
    assert rc == 0
    assert (tmp_path / "result.npz").exists()
