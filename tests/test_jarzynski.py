"""Tests for Jarzynski free-energy estimators."""

from __future__ import annotations

import numpy as np
import pytest

from src.analyze.jarzynski import evaluate_convergence, jarzynski_free_energy
from src.config import BOLTZMANN_KJ


def _make_gaussian_work_values(
    *,
    n_samples: int = 6000,
    temperature_k: float = 310.0,
    delta_g_kj_mol: float = 12.0,
    sigma_kj_mol: float = 4.0,
) -> tuple[np.ndarray, float, float, float]:
    """Generate Gaussian work data with a known analytical Jarzynski free energy."""

    rng = np.random.default_rng(20260312)
    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)
    mean_work_kj_mol = delta_g_kj_mol + 0.5 * beta * sigma_kj_mol**2
    work_values = rng.normal(loc=mean_work_kj_mol, scale=sigma_kj_mol, size=n_samples)
    return work_values, temperature_k, delta_g_kj_mol, sigma_kj_mol


def test_jarzynski_free_energy_recovers_analytical_gaussian_result() -> None:
    """Chunk 19 gate: both Jarzynski estimators should recover the analytical Gaussian free energy."""

    work_values, temperature_k, analytical_delta_g_kj_mol, sigma_kj_mol = _make_gaussian_work_values()
    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)

    result = jarzynski_free_energy(work_values, temperature_k)

    assert result["mean_work_kj_mol"] == pytest.approx(np.mean(work_values), abs=1e-10)
    assert result["work_variance_kj_mol2"] == pytest.approx(np.var(work_values), rel=1e-12)
    assert result["delta_g_kj_mol"] == pytest.approx(analytical_delta_g_kj_mol, abs=0.2)
    assert result["delta_g_cumulant2_kj_mol"] == pytest.approx(analytical_delta_g_kj_mol, abs=0.2)
    expected_cumulant = np.mean(work_values) - 0.5 * beta * np.var(work_values)
    assert result["delta_g_cumulant2_kj_mol"] == pytest.approx(expected_cumulant, abs=1e-10)
    assert abs(result["delta_g_kj_mol"] - result["delta_g_cumulant2_kj_mol"]) < 0.2


def test_evaluate_convergence_returns_expected_schema_and_final_estimate() -> None:
    """Convergence diagnostics should return bounded arrays and approach the analytical free energy."""

    work_values, temperature_k, analytical_delta_g_kj_mol, _ = _make_gaussian_work_values(n_samples=5000)

    result = evaluate_convergence(work_values, temperature_k, n_subsets=12)

    assert result["subset_sizes"].ndim == 1
    assert result["delta_g_vs_n"].shape == result["subset_sizes"].shape
    assert result["std_vs_n"].shape == result["subset_sizes"].shape
    assert np.all(result["subset_sizes"][:-1] < result["subset_sizes"][1:])
    assert np.all(result["std_vs_n"] >= 0.0)
    assert result["subset_sizes"][-1] == work_values.size
    assert result["delta_g_vs_n"][-1] == pytest.approx(analytical_delta_g_kj_mol, abs=0.25)


def test_jarzynski_free_energy_rejects_invalid_inputs() -> None:
    """Public API should reject empty work arrays and non-physical temperatures."""

    with pytest.raises(ValueError, match="non-empty"):
        jarzynski_free_energy(np.asarray([], dtype=float), 310.0)

    with pytest.raises(ValueError, match="positive"):
        evaluate_convergence(np.asarray([1.0, 2.0, 3.0], dtype=float), 0.0, n_subsets=5)