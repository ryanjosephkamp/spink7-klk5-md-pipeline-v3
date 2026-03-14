"""Tests for collective variable utilities."""

import numpy as np
import pytest

from src.physics.collective_variables import com_distance, com_vector


@pytest.fixture
def simple_positions() -> np.ndarray:
    """Provide a small coordinate array with analytically tractable COM values."""

    return np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ],
        dtype=float,
    )


@pytest.fixture
def simple_masses() -> np.ndarray:
    """Provide masses paired with the simple test coordinates."""

    return np.array([1.0, 3.0, 2.0, 4.0], dtype=float)


def test_com_vector_matches_manual_calculation(simple_positions: np.ndarray, simple_masses: np.ndarray) -> None:
    """COM vector should agree with a manual mass-weighted average."""

    result = com_vector(simple_positions, simple_masses, np.array([0, 1]))
    expected = np.array([1.5, 0.0, 0.0], dtype=float)

    assert np.allclose(result, expected, atol=1e-12)


def test_com_distance_matches_manual_calculation(simple_positions: np.ndarray, simple_masses: np.ndarray) -> None:
    """COM distance should match the Euclidean distance between manual COM vectors."""

    result = com_distance(simple_positions, simple_masses, np.array([0, 1]), np.array([2, 3]))
    expected_group_a = np.array([1.5, 0.0, 0.0], dtype=float)
    expected_group_b = np.array([0.0, 1.0, 8.0 / 3.0], dtype=float)
    expected_distance = float(np.linalg.norm(expected_group_a - expected_group_b))

    assert result == pytest.approx(expected_distance, abs=1e-12)


def test_com_distance_is_symmetric(simple_positions: np.ndarray, simple_masses: np.ndarray) -> None:
    """COM distance should be invariant to swapping the input groups."""

    forward = com_distance(simple_positions, simple_masses, np.array([0, 1]), np.array([2, 3]))
    reverse = com_distance(simple_positions, simple_masses, np.array([2, 3]), np.array([0, 1]))

    assert forward == pytest.approx(reverse, abs=1e-12)


@pytest.mark.parametrize(
    "group_indices,error_message",
    [
        (np.array([], dtype=int), "group_indices must be non-empty"),
        (np.array([4], dtype=int), r"group_indices must be within \[0, N_atoms\)"),
        (np.array([-1], dtype=int), r"group_indices must be within \[0, N_atoms\)"),
    ],
)
def test_com_vector_rejects_invalid_group_indices(
    simple_positions: np.ndarray,
    simple_masses: np.ndarray,
    group_indices: np.ndarray,
    error_message: str,
) -> None:
    """Public API should reject empty or out-of-bounds atom selections."""

    with pytest.raises(ValueError, match=error_message):
        com_vector(simple_positions, simple_masses, group_indices)


def test_com_vector_rejects_non_positive_masses(simple_positions: np.ndarray) -> None:
    """Public API should reject non-physical atomic masses."""

    invalid_masses = np.array([1.0, 0.0, 2.0, 4.0], dtype=float)

    with pytest.raises(ValueError, match="masses must be strictly positive"):
        com_vector(simple_positions, invalid_masses, np.array([0, 1]))
