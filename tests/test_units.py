"""Tests for physics unit conversion utilities."""

import pytest

from src.physics.units import kbt, kj_to_kcal, nm_to_angstrom, ps_to_ns


def test_kj_to_kcal_converts_to_six_decimal_places() -> None:
    """Energy conversion should match the configured kcal-to-kJ factor."""

    assert kj_to_kcal(4.184) == pytest.approx(1.0, abs=1e-6)
    assert kj_to_kcal(41.84) == pytest.approx(10.0, abs=1e-6)


def test_nm_to_angstrom_converts_to_six_decimal_places() -> None:
    """Length conversion should preserve the exact nanometer-to-angstrom ratio."""

    assert nm_to_angstrom(1.0) == pytest.approx(10.0, abs=1e-6)
    assert nm_to_angstrom(0.14) == pytest.approx(1.4, abs=1e-6)


def test_ps_to_ns_converts_to_six_decimal_places() -> None:
    """Time conversion should preserve the exact picosecond-to-nanosecond ratio."""

    assert ps_to_ns(1000.0) == pytest.approx(1.0, abs=1e-6)
    assert ps_to_ns(250.0) == pytest.approx(0.25, abs=1e-6)


def test_kbt_returns_thermal_energy_to_six_decimal_places() -> None:
    """Thermal energy should equal the Boltzmann constant times temperature."""

    assert kbt(310.0) == pytest.approx(2.57748341158, abs=1e-6)
    assert kbt(273.15) == pytest.approx(2.2710954641066996, abs=1e-6)


def test_kbt_rejects_non_physical_temperature() -> None:
    """Temperature input must be strictly positive at the public API boundary."""

    with pytest.raises(ValueError, match="temperature_k must be positive"):
        kbt(0.0)

    with pytest.raises(ValueError, match="temperature_k must be positive"):
        kbt(-1.0)
