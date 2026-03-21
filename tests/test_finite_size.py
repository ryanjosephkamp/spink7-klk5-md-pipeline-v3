"""Tests for finite-size electrostatic correction utilities."""

from __future__ import annotations

import numpy as np
import pytest

from src.config import FiniteSizeCorrectionConfig
from src.physics.finite_size import (
    extrapolate_to_infinite_box,
    finite_size_correction,
    hunenberger_mccammon_correction,
    run_box_size_convergence_study,
)


def test_correction_scales_as_inverse_L() -> None:
    """Finite-size correction must scale as 1/L for constant charge."""

    charge = 3.0
    box_lengths = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    corrections = [hunenberger_mccammon_correction(charge, L) for L in box_lengths]

    products = [c * L for c, L in zip(corrections, box_lengths)]
    for product in products:
        assert product == pytest.approx(products[0], rel=1e-10), (
            "Correction * L must be constant for 1/L scaling"
        )


def test_correction_scales_as_Q_squared() -> None:
    """Finite-size correction must scale as Q^2 for constant box size."""

    box_length = 5.0
    charges = [1.0, 2.0, 3.0, 4.0]
    corrections = [hunenberger_mccammon_correction(q, box_length) for q in charges]

    for q, c in zip(charges, corrections):
        expected = corrections[0] * q**2
        assert c == pytest.approx(expected, rel=1e-10), (
            f"Correction for Q={q} should be Q^2 * correction(Q=1)"
        )


def test_zero_charge_gives_zero_correction() -> None:
    """Zero net charge should produce exactly zero correction."""

    assert hunenberger_mccammon_correction(0.0, 5.0) == pytest.approx(0.0, abs=1e-15)
    assert hunenberger_mccammon_correction(0.0, 10.0) == pytest.approx(0.0, abs=1e-15)


def test_non_cubic_box_uses_geometric_mean() -> None:
    """Non-cubic boxes should compute effective length as geometric mean."""

    box_vectors = np.array([
        [4.0, 0.0, 0.0],
        [0.0, 5.0, 0.0],
        [0.0, 0.0, 6.0],
    ])
    result = finite_size_correction(1.0, box_vectors)

    expected_length = float(np.cbrt(4.0 * 5.0 * 6.0))
    assert result["effective_box_length_nm"] == pytest.approx(expected_length, rel=1e-10)

    expected_correction = hunenberger_mccammon_correction(1.0, expected_length)
    assert result["correction_kj_mol"] == pytest.approx(expected_correction, rel=1e-10)
    assert result["correction_kcal_mol"] == pytest.approx(
        expected_correction / 4.184, rel=1e-10
    )


def test_extrapolation_with_noisy_data() -> None:
    """Extrapolation should recover the correct value within tolerance for noisy 1/L data."""

    rng = np.random.default_rng(42)
    box_lengths = np.array([3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0])
    true_e_inf = -200.0
    true_slope = 8.0
    energies = true_e_inf + true_slope / box_lengths + rng.normal(0, 0.02, len(box_lengths))

    result = extrapolate_to_infinite_box(box_lengths, energies)
    assert result["energy_infinite_kj_mol"] == pytest.approx(true_e_inf, abs=0.5), (
        f"Extrapolated energy should be near {true_e_inf}, got {result['energy_infinite_kj_mol']}"
    )
    assert result["r_squared"] > 0.95, (
        f"R^2 should be high even with noise, got {result['r_squared']}"
    )


def test_correction_validation_rejects_non_positive_box_length() -> None:
    """Correction functions should reject non-positive box lengths."""

    with pytest.raises(ValueError, match="box_length_nm must be positive"):
        hunenberger_mccammon_correction(1.0, 0.0)

    with pytest.raises(ValueError, match="box_length_nm must be positive"):
        hunenberger_mccammon_correction(1.0, -5.0)


def test_hunenberger_mccammon_known_values() -> None:
    """HM correction should match analytically computed values for known inputs."""

    correction = hunenberger_mccammon_correction(net_charge_e=1.0, box_length_nm=5.0)
    expected = -(1.0**2 * (-2.837297) * 138.935456) / (2.0 * 80.0 * 5.0)
    assert correction == pytest.approx(expected, rel=1e-6)


def test_extrapolation_rejects_fewer_than_three_points() -> None:
    """Extrapolation requires at least 3 data points."""

    with pytest.raises(ValueError, match="At least 3 data points"):
        extrapolate_to_infinite_box(np.array([3.0, 4.0]), np.array([-100.0, -99.0]))


def test_convergence_study_synthetic_1_over_L() -> None:
    """Convergence study should recover infinite-box energy from synthetic 1/L data."""

    energies_dict = {L: -100.0 + 5.0 / L for L in [3.0, 4.0, 5.0, 6.0, 8.0, 10.0]}
    study = run_box_size_convergence_study(energies_dict)
    assert study["energy_infinite_kj_mol"] == pytest.approx(-100.0, abs=0.1)
    assert study["r_squared"] > 0.99
    assert study["scaling_exponent"] == pytest.approx(1.0, abs=0.3)
