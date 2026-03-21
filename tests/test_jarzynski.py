"""Tests for Jarzynski free-energy estimators."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from src.analyze.jarzynski import bar_free_energy, diagnose_dissipation, evaluate_convergence, jarzynski_free_energy
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


def test_diagnose_dissipation_low_regime() -> None:
    """Dissipation diagnostic should not warn for near-equilibrium work data."""

    rng = np.random.default_rng(42)
    beta = 1.0 / (BOLTZMANN_KJ * 310.0)
    delta_g = 10.0  # kJ/mol
    w_diss_target = 1.0  # kJ/mol, well below 3 kBT ~ 7.73 kJ/mol
    mean_work = delta_g + w_diss_target
    sigma = 2.0
    work = rng.normal(loc=mean_work, scale=sigma, size=50_000)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = diagnose_dissipation(work, 310.0)
        assert not result["warning_triggered"]
        assert len(w) == 0

    assert result["w_diss_over_kbt"] == pytest.approx(beta * w_diss_target, abs=0.5)


def test_diagnose_dissipation_high_regime_warns() -> None:
    """Dissipation diagnostic must warn when W_diss > 3 kBT."""

    rng = np.random.default_rng(43)
    # For Gaussian work, W_diss / kBT = beta^2 * sigma^2 / 2.
    # With beta ~ 0.388 mol/kJ at 310 K and sigma = 10 kJ/mol,
    # W_diss / kBT ~ 7.5.
    sigma = 10.0
    mean_work = 30.0
    work = rng.normal(loc=mean_work, scale=sigma, size=50_000)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = diagnose_dissipation(work, 310.0)
        assert result["warning_triggered"]
        assert len(w) == 1
        assert "dissipation" in str(w[0].message).lower()

    assert result["w_diss_over_kbt"] > 3.0
    assert 0.0 < result["efficiency_ratio"] < 1.0


def test_normality_test_passes_for_gaussian_work() -> None:
    """Normality test should pass for genuinely Gaussian work data."""

    work, temp, _, _ = _make_gaussian_work_values(n_samples=6000)
    result = jarzynski_free_energy(work, temp)
    assert result["cumulant2_normality_passed"] is True


def test_normality_test_rejects_bimodal_work() -> None:
    """Normality test must reject bimodal (non-Gaussian) work distributions."""

    rng = np.random.default_rng(44)
    mode_a = rng.normal(loc=30.0, scale=2.0, size=3000)
    mode_b = rng.normal(loc=60.0, scale=2.0, size=3000)
    bimodal_work = np.concatenate([mode_a, mode_b])

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = jarzynski_free_energy(bimodal_work, 310.0)
        assert result["cumulant2_normality_passed"] is False
        normality_warnings = [x for x in w if "cumulant" in str(x.message).lower()]
        assert len(normality_warnings) >= 1


def test_bar_recovers_exact_free_energy_gaussian() -> None:
    """BAR must recover the exact free energy from Gaussian forward+reverse work."""

    rng = np.random.default_rng(45)
    temperature_k = 310.0
    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)
    delta_g_true = 20.0  # kJ/mol
    sigma = 4.0

    # Crooks-consistent Gaussian distributions: mu_F = dG + beta*sigma^2/2,
    # mu_R = -dG + beta*sigma^2/2.
    w_diss = 0.5 * beta * sigma**2
    forward_work = rng.normal(loc=delta_g_true + w_diss, scale=sigma, size=5000)
    reverse_work = rng.normal(loc=-delta_g_true + w_diss, scale=sigma, size=5000)

    result = bar_free_energy(forward_work, reverse_work, temperature_k)
    assert result["converged"] is True
    assert result["delta_g_bar_kj_mol"] == pytest.approx(delta_g_true, abs=0.5)


def test_bar_outperforms_jarzynski_at_high_dissipation() -> None:
    """At high dissipation, BAR should be closer to the true dG than Jarzynski."""

    rng = np.random.default_rng(46)
    temperature_k = 310.0
    beta = 1.0 / (BOLTZMANN_KJ * temperature_k)
    delta_g_true = 20.0
    # Large sigma gives high dissipation: W_diss/kBT = beta^2*sigma^2/2 ~ 7.5
    sigma = 10.0
    n = 500  # Deliberately small N

    w_diss = 0.5 * beta * sigma**2
    forward_work = rng.normal(loc=delta_g_true + w_diss, scale=sigma, size=n)
    reverse_work = rng.normal(loc=-delta_g_true + w_diss, scale=sigma, size=n)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        jar = jarzynski_free_energy(forward_work, temperature_k)
    bar = bar_free_energy(forward_work, reverse_work, temperature_k)

    jar_error = abs(jar["delta_g_kj_mol"] - delta_g_true)
    bar_error = abs(bar["delta_g_bar_kj_mol"] - delta_g_true)
    assert bar_error < jar_error, (
        f"BAR error ({bar_error:.2f}) should be less than Jarzynski error ({jar_error:.2f})"
    )


def test_jarzynski_free_energy_returns_dissipation_keys() -> None:
    """The public API must include dissipation diagnostic fields."""

    work, temp, _, _ = _make_gaussian_work_values()
    result = jarzynski_free_energy(work, temp)
    assert "w_diss_kj_mol" in result
    assert "w_diss_over_kbt" in result
    assert "efficiency_ratio" in result
    assert "cumulant2_normality_passed" in result
    assert isinstance(result["w_diss_over_kbt"], float)
    assert isinstance(result["cumulant2_normality_passed"], bool)


def test_bar_symmetric_distributions_yield_zero_delta_g() -> None:
    """Symmetric forward/reverse work around dG=0 must yield BAR estimate near zero."""

    rng = np.random.default_rng(47)
    sigma = 3.0
    beta = 1.0 / (BOLTZMANN_KJ * 310.0)
    # Crooks-consistent with dG=0: mu_F = mu_R = beta*sigma^2/2
    mu = 0.5 * beta * sigma**2
    forward = rng.normal(loc=mu, scale=sigma, size=5000)
    reverse = rng.normal(loc=mu, scale=sigma, size=5000)
    result = bar_free_energy(forward, reverse, 310.0)
    assert result["delta_g_bar_kj_mol"] == pytest.approx(0.0, abs=0.5)