"""Unit conversion helpers for the SPINK7-KLK5 MD pipeline."""

from src.config import BOLTZMANN_KJ, KCAL_TO_KJ


def kj_to_kcal(value: float) -> float:
    """Convert an energy value from kJ/mol to kcal/mol.

    Invariants: None.

    Args:
        value: Energy value in kJ/mol.

    Returns:
        Energy value in kcal/mol.
    """

    return value / KCAL_TO_KJ


def nm_to_angstrom(value: float) -> float:
    """Convert a length value from nanometers to angstroms.

    Invariants: None.

    Args:
        value: Length value in nm.

    Returns:
        Length value in angstroms.
    """

    return value * 10.0


def ps_to_ns(value: float) -> float:
    """Convert a time value from picoseconds to nanoseconds.

    Invariants: None.

    Args:
        value: Time value in ps.

    Returns:
        Time value in ns.
    """

    return value / 1000.0


def kbt(temperature_k: float) -> float:
    """Compute the thermal energy $kT$ in kJ/mol at the given temperature.

    Invariants: None.

    Args:
        temperature_k: Absolute temperature in K.

    Returns:
        Thermal energy in kJ/mol.

    Raises:
        ValueError: If the temperature is not strictly positive.
    """

    if temperature_k <= 0.0:
        raise ValueError("temperature_k must be positive")

    return BOLTZMANN_KJ * temperature_k
