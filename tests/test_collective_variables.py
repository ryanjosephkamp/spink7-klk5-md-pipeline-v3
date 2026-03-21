"""Tests for collective variable utilities."""

import numpy as np
import pytest

from src.physics.collective_variables import (
    com_distance,
    com_vector,
    com_angle,
    contact_fraction,
    CollectiveVariableSpec,
    compute_cv,
)


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


def test_com_vector_pbc_single_axis_wrap() -> None:
    """COM should be correct when a group wraps around one periodic boundary."""
    box_lengths = np.array([5.0, 5.0, 5.0])
    # Two equal-mass atoms at x=0.1 and x=4.9  ->  true COM at x=0.0 (i.e., 5.0 equiv 0.0)
    positions = np.array([[0.1, 1.0, 1.0],
                          [4.9, 1.0, 1.0]])
    masses = np.array([1.0, 1.0])
    group = np.array([0, 1])

    com = com_vector(positions, masses, group, box_lengths=box_lengths)

    # The true COM x-coordinate should be at 0.0 (equivalently 5.0), not 2.5
    # Due to modular arithmetic, result should be near 0.0 or 5.0
    actual_x = com[0] % box_lengths[0]
    assert abs(actual_x) < 1e-6 or abs(actual_x - box_lengths[0]) < 1e-6, \
        f"PBC COM x={com[0]:.6f}, expected near 0.0 or 5.0"
    assert abs(com[1] - 1.0) < 1e-6, "y-coordinate should be unchanged"
    assert abs(com[2] - 1.0) < 1e-6, "z-coordinate should be unchanged"


def test_com_distance_pbc_cross_boundary() -> None:
    """COM distance must be correct when groups are on opposite sides of a box edge."""
    box_lengths = np.array([5.0, 5.0, 5.0])
    # Group A: single atom at (4.8, 2.5, 2.5)
    # Group B: single atom at (0.2, 2.5, 2.5)
    # Naive distance in x: |4.8 - 0.2| = 4.6
    # PBC-corrected distance in x: 5.0 - 4.6 = 0.4
    positions = np.array([[4.8, 2.5, 2.5],
                          [0.2, 2.5, 2.5]])
    masses = np.array([1.0, 1.0])
    group_a = np.array([0])
    group_b = np.array([1])

    dist_pbc = com_distance(positions, masses, group_a, group_b, box_lengths=box_lengths)
    dist_naive = com_distance(positions, masses, group_a, group_b, box_lengths=None)

    assert dist_pbc == pytest.approx(0.4, abs=1e-6), \
        f"PBC distance should be 0.4 nm, got {dist_pbc}"
    assert dist_naive == pytest.approx(4.6, abs=1e-6), \
        f"Naive distance should be 4.6 nm, got {dist_naive}"
    assert dist_naive > box_lengths[0] / 4, \
        "Naive distance should be demonstrably incorrect (> L/4)"


def test_com_vector_pbc_multi_axis_wrap() -> None:
    """COM must be correct when a group wraps across two or more periodic boundaries."""
    box_lengths = np.array([4.0, 4.0, 4.0])
    # Atom at (0.1, 0.1, 2.0) and atom at (3.9, 3.9, 2.0)
    # True COM should be near (0.0, 0.0, 2.0) -- wrapping in x and y
    positions = np.array([[0.1, 0.1, 2.0],
                          [3.9, 3.9, 2.0]])
    masses = np.array([1.0, 1.0])
    group = np.array([0, 1])

    com = com_vector(positions, masses, group, box_lengths=box_lengths)
    com_wrapped = com % box_lengths

    assert abs(com_wrapped[0]) < 0.15 or abs(com_wrapped[0] - 4.0) < 0.15, \
        f"PBC COM x={com_wrapped[0]:.6f}, expected near 0.0 or 4.0"
    assert abs(com_wrapped[1]) < 0.15 or abs(com_wrapped[1] - 4.0) < 0.15, \
        f"PBC COM y={com_wrapped[1]:.6f}, expected near 0.0 or 4.0"
    assert abs(com_wrapped[2] - 2.0) < 1e-6, \
        f"z-coordinate should be 2.0, got {com_wrapped[2]:.6f}"


def test_com_vector_pbc_no_wrap_matches_naive() -> None:
    """When no atoms wrap, PBC-aware COM should match naive COM."""
    # Box must be much larger than position spread for angular mean to converge
    # to the arithmetic mean within tight tolerance (error is O(theta^2)).
    box_lengths = np.array([1000.0, 1000.0, 1000.0])
    positions = np.array([[1.0, 2.0, 3.0],
                          [2.0, 3.0, 4.0],
                          [1.5, 2.5, 3.5]])
    masses = np.array([12.0, 16.0, 14.0])
    group = np.array([0, 1, 2])

    com_pbc = com_vector(positions, masses, group, box_lengths=box_lengths)
    com_naive = com_vector(positions, masses, group, box_lengths=None)

    np.testing.assert_allclose(com_pbc, com_naive, atol=1e-6)


def test_com_vector_pbc_unequal_masses() -> None:
    """PBC COM must correctly weight by mass during circular averaging."""
    box_lengths = np.array([6.0, 6.0, 6.0])
    # Heavy atom at x=0.2, light atom at x=5.8 -> COM should be closer to x=0.2
    positions = np.array([[0.2, 3.0, 3.0],
                          [5.8, 3.0, 3.0]])
    masses = np.array([10.0, 1.0])
    group = np.array([0, 1])

    com = com_vector(positions, masses, group, box_lengths=box_lengths)
    com_x = com[0] % box_lengths[0]

    # With 10:1 mass ratio, COM should be much closer to x=0.2 than x=5.8
    assert com_x < 0.5 or com_x > 5.5, \
        f"Mass-weighted PBC COM should be near x=0.2, got x={com_x:.4f}"


def test_com_vector_rejects_invalid_box_lengths() -> None:
    """box_lengths must be shape [3] with all positive values."""
    positions = np.array([[1.0, 2.0, 3.0]])
    masses = np.array([1.0])
    group = np.array([0])

    with pytest.raises(ValueError, match="box_lengths must have shape"):
        com_vector(positions, masses, group, box_lengths=np.array([5.0, 5.0]))

    with pytest.raises(ValueError, match="box_lengths must be strictly positive"):
        com_vector(positions, masses, group, box_lengths=np.array([5.0, 0.0, 5.0]))


# ---------- L-08 Step 1: com_angle tests ----------


def test_com_angle_aligned_returns_zero() -> None:
    """When COM axis is aligned with reference, angle should be zero."""
    positions = np.array([
        [2.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=float)
    masses = np.array([1.0, 1.0], dtype=float)
    ref_dir = np.array([1.0, 0.0, 0.0], dtype=float)
    angle = com_angle(positions, masses, np.array([0]), np.array([1]), ref_dir)
    assert angle == pytest.approx(0.0, abs=1e-12)


def test_com_angle_perpendicular_returns_pi_over_2() -> None:
    """When COM axis is perpendicular to reference, angle should be pi/2."""
    positions = np.array([
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 0.0],
    ], dtype=float)
    masses = np.array([1.0, 1.0], dtype=float)
    ref_dir = np.array([1.0, 0.0, 0.0], dtype=float)
    angle = com_angle(positions, masses, np.array([0]), np.array([1]), ref_dir)
    assert angle == pytest.approx(np.pi / 2, abs=1e-12)


def test_com_angle_antiparallel_returns_pi() -> None:
    """When COM axis is antiparallel to reference, angle should be pi."""
    positions = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
    ], dtype=float)
    masses = np.array([1.0, 1.0], dtype=float)
    ref_dir = np.array([1.0, 0.0, 0.0], dtype=float)
    angle = com_angle(positions, masses, np.array([0]), np.array([1]), ref_dir)
    assert angle == pytest.approx(np.pi, abs=1e-12)


def test_com_angle_rejects_zero_reference_direction() -> None:
    """Public API should reject zero-norm reference direction."""
    positions = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=float)
    masses = np.array([1.0, 1.0], dtype=float)
    with pytest.raises(ValueError, match="non-zero"):
        com_angle(positions, masses, np.array([0]), np.array([1]),
                  np.array([0.0, 0.0, 0.0]))


# ---------- L-08 Step 2: contact_fraction tests ----------


def test_contact_fraction_all_contacts_intact() -> None:
    """When all pairs are at their reference distance, Q should be close to 1."""
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.3, 0.0, 0.0],
        [0.0, 0.4, 0.0],
        [0.3, 0.4, 0.0],
    ], dtype=float)
    pairs = np.array([[0, 1], [2, 3]], dtype=int)
    ref_dists = np.array([0.3, 0.3], dtype=float)
    q = contact_fraction(positions, pairs, ref_dists)
    assert q > 0.95, f"Expected Q > 0.95 for intact contacts, got {q}"


def test_contact_fraction_all_contacts_broken() -> None:
    """When all pairs are far apart, Q should be close to 0."""
    positions = np.array([
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [0.0, 10.0, 0.0],
        [10.0, 10.0, 0.0],
    ], dtype=float)
    pairs = np.array([[0, 1], [2, 3]], dtype=int)
    ref_dists = np.array([0.3, 0.3], dtype=float)
    q = contact_fraction(positions, pairs, ref_dists)
    assert q < 0.01, f"Expected Q < 0.01 for broken contacts, got {q}"


def test_contact_fraction_known_analytic_value() -> None:
    """Verify Q against a manually computed switching value."""
    # Single pair at exactly r_0 distance: switching = 1/(1 + 1^12) = 0.5
    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.45, 0.0, 0.0],
    ], dtype=float)
    pairs = np.array([[0, 1]], dtype=int)
    ref_dists = np.array([0.3], dtype=float)  # r_0 = 1.5 * 0.3 = 0.45
    q = contact_fraction(positions, pairs, ref_dists, distance_factor=1.5)
    assert q == pytest.approx(0.5, abs=1e-10)


def test_contact_fraction_rejects_empty_pairs() -> None:
    """Public API should reject empty contact pair list."""
    positions = np.array([[0.0, 0.0, 0.0]], dtype=float)
    with pytest.raises(ValueError):
        contact_fraction(positions, np.empty((0, 2), dtype=int),
                         np.empty(0, dtype=float))


def test_com_distance_backward_compatible_without_box_lengths(
    simple_positions: np.ndarray, simple_masses: np.ndarray
) -> None:
    """Existing calls without box_lengths must produce identical results."""
    dist_explicit_none = com_distance(
        simple_positions, simple_masses,
        np.array([0, 1]), np.array([2, 3]),
        box_lengths=None,
    )
    dist_default = com_distance(
        simple_positions, simple_masses,
        np.array([0, 1]), np.array([2, 3]),
    )
    assert dist_explicit_none == pytest.approx(dist_default, abs=1e-15)


def test_com_distance_pbc_is_symmetric() -> None:
    """PBC-corrected COM distance must be symmetric: d(A,B) == d(B,A)."""
    box_lengths = np.array([5.0, 5.0, 5.0])
    positions = np.array([[4.8, 2.5, 2.5],
                          [0.2, 2.5, 2.5]])
    masses = np.array([1.0, 1.0])

    d_ab = com_distance(positions, masses, np.array([0]), np.array([1]), box_lengths)
    d_ba = com_distance(positions, masses, np.array([1]), np.array([0]), box_lengths)
    assert d_ab == pytest.approx(d_ba, abs=1e-15)


# ---------- L-08 Step 3: CV registry / dispatcher tests ----------


def test_compute_cv_dispatches_com_distance() -> None:
    """compute_cv with 'com_distance' matches direct com_distance call."""
    positions = np.array([
        [0.0, 0.0, 0.0], [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0], [0.0, 0.0, 4.0],
    ], dtype=float)
    masses = np.array([1.0, 3.0, 2.0, 4.0], dtype=float)
    spec = CollectiveVariableSpec(name="com_distance")
    result = compute_cv(spec, positions, masses, np.array([0, 1]), np.array([2, 3]))
    expected = com_distance(positions, masses, np.array([0, 1]), np.array([2, 3]))
    assert result == pytest.approx(expected, abs=1e-12)


def test_compute_cv_dispatches_com_angle() -> None:
    """compute_cv with 'com_angle' matches direct com_angle call."""
    positions = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=float)
    masses = np.array([1.0, 1.0], dtype=float)
    ref = np.array([1.0, 0.0, 0.0], dtype=float)
    spec = CollectiveVariableSpec(name="com_angle", reference_direction=ref)
    result = compute_cv(spec, positions, masses, np.array([0]), np.array([1]))
    assert result == pytest.approx(0.0, abs=1e-12)


def test_compute_cv_dispatches_contact_fraction() -> None:
    """compute_cv with 'contact_fraction' matches direct contact_fraction call."""
    positions = np.array([
        [0.0, 0.0, 0.0], [0.45, 0.0, 0.0],
    ], dtype=float)
    masses = np.array([1.0, 1.0], dtype=float)
    pairs = np.array([[0, 1]], dtype=int)
    ref_dists = np.array([0.3], dtype=float)
    spec = CollectiveVariableSpec(
        name="contact_fraction",
        contact_pairs=pairs,
        contact_reference_distances_nm=ref_dists,
    )
    result = compute_cv(spec, positions, masses, np.array([0]), np.array([1]))
    assert result == pytest.approx(0.5, abs=1e-10)


def test_compute_cv_rejects_unknown_name() -> None:
    """compute_cv should raise ValueError for unrecognized CV names."""
    positions = np.array([[0.0, 0.0, 0.0]], dtype=float)
    masses = np.array([1.0], dtype=float)
    spec = CollectiveVariableSpec(name="nonexistent_cv")
    with pytest.raises(ValueError, match="Unknown collective variable"):
        compute_cv(spec, positions, masses, np.array([0]), np.array([0]))
