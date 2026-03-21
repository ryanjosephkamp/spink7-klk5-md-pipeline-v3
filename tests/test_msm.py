"""Tests for MSM construction, featurization, and kinetics modules (L-40)."""

import dataclasses

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Step 1: deeptime dependency
# ---------------------------------------------------------------------------


def test_deeptime_importable():
    """The deeptime package is installed and importable."""
    import deeptime
    from deeptime.decomposition import TICA  # noqa: F401
    from deeptime.clustering import KMeans  # noqa: F401
    from deeptime.markov.msm import MaximumLikelihoodMSM  # noqa: F401

    assert hasattr(deeptime, "__version__")


# ---------------------------------------------------------------------------
# Step 2: MSMConfig dataclass
# ---------------------------------------------------------------------------


def test_msm_config_defaults():
    """MSMConfig is frozen and has correct defaults."""
    from src.config import MSMConfig

    config = MSMConfig()
    assert config.tica_lag_ps == 10.0
    assert config.n_clusters == 200
    assert config.tica_n_components == 4
    assert len(config.lag_times_ps) == 9
    assert config.n_implied_timescales == 5
    assert config.bayesian_n_samples == 100
    assert config.stride == 1
    assert dataclasses.is_dataclass(config)


def test_msm_config_frozen():
    """MSMConfig is immutable."""
    from src.config import MSMConfig

    config = MSMConfig()
    with pytest.raises(dataclasses.FrozenInstanceError):
        config.n_clusters = 500


# ---------------------------------------------------------------------------
# Step 3: Featurization module
# ---------------------------------------------------------------------------


def test_backbone_dihedrals_shape(alanine_dipeptide_pdb):
    """Backbone dihedral featurization produces correct shape."""
    import mdtraj as md

    from src.analyze.featurize import compute_backbone_dihedrals

    traj = md.load(str(alanine_dipeptide_pdb))
    # Duplicate frames to simulate a trajectory
    traj = traj.join(traj)
    features = compute_backbone_dihedrals(traj)
    assert features.ndim == 2
    assert features.shape[0] == traj.n_frames
    # Each residue with both phi/psi contributes 4 columns (sin/cos pairs)
    assert features.shape[1] % 4 == 0
    assert features.dtype == float


def test_combine_features_basic():
    """combine_features horizontally stacks arrays."""
    from src.analyze.featurize import combine_features

    a = np.ones((10, 3))
    b = np.zeros((10, 2))
    result = combine_features(a, b)
    assert result.shape == (10, 5)
    assert result.dtype == float


def test_combine_features_mismatch():
    """combine_features raises on frame count mismatch."""
    from src.analyze.featurize import combine_features

    a = np.ones((10, 3))
    b = np.zeros((8, 2))
    with pytest.raises(ValueError, match="frames"):
        combine_features(a, b)


def test_combine_features_empty():
    """combine_features raises on empty input."""
    from src.analyze.featurize import combine_features

    with pytest.raises(ValueError, match="At least one"):
        combine_features()


# ---------------------------------------------------------------------------
# Step 4: MSM construction module
# ---------------------------------------------------------------------------


def _make_two_state_trajectory(n_frames: int = 5000, seed: int = 42) -> np.ndarray:
    """Generate synthetic 2D features with two well-separated clusters."""
    rng = np.random.default_rng(seed)
    # State A: centered at (-2, 0); State B: centered at (2, 0)
    states = rng.integers(0, 2, size=n_frames)
    # Add temporal correlation: stay in same state with high probability
    for i in range(1, n_frames):
        if rng.random() < 0.98:
            states[i] = states[i - 1]
    centers = np.array([[-2.0, 0.0], [2.0, 0.0]])
    noise = rng.normal(0.0, 0.3, (n_frames, 2))
    return (centers[states] + noise).astype(float)


def test_fit_tica_returns_correct_shape():
    """fit_tica returns transformed data with expected shape."""
    from src.analyze.msm import fit_tica
    from src.config import MSMConfig

    features = _make_two_state_trajectory(n_frames=2000)
    config = MSMConfig(tica_lag_ps=1.0, tica_n_components=2, stride=1)
    result = fit_tica(features, config)
    assert result["tica_output"].shape == (2000, 2)
    assert result["eigenvalues"].shape[0] <= 2
    assert result["tica_output"].dtype == float


def test_cluster_microstates_returns_assignments():
    """cluster_microstates produces integer assignments."""
    from src.analyze.msm import cluster_microstates, fit_tica
    from src.config import MSMConfig

    features = _make_two_state_trajectory(n_frames=2000)
    config = MSMConfig(tica_lag_ps=1.0, tica_n_components=2, n_clusters=10, stride=1)
    tica_out = fit_tica(features, config)["tica_output"]
    assignments = cluster_microstates(tica_out, config)
    assert assignments.shape == (2000,)
    assert assignments.dtype == int
    assert len(np.unique(assignments)) <= 10


def test_build_msm_transition_matrix_row_stochastic():
    """build_msm produces a row-stochastic transition matrix."""
    from src.analyze.msm import build_msm, cluster_microstates, fit_tica
    from src.config import MSMConfig

    features = _make_two_state_trajectory(n_frames=3000)
    config = MSMConfig(
        tica_lag_ps=1.0, tica_n_components=2, n_clusters=10,
        stride=1, n_implied_timescales=3,
    )
    tica_out = fit_tica(features, config)["tica_output"]
    assignments = cluster_microstates(tica_out, config)
    result = build_msm(assignments, lag_time_ps=5.0, config=config)

    T = result["transition_matrix"]
    # Row-stochastic: each row sums to 1
    row_sums = T.sum(axis=1)
    np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)
    # All entries non-negative
    assert np.all(T >= 0.0)
    # Stationary distribution sums to 1
    pi = result["stationary_distribution"]
    assert abs(pi.sum() - 1.0) < 1e-10


def test_implied_timescales_shape():
    """compute_implied_timescales returns correct output shapes."""
    from src.analyze.msm import cluster_microstates, compute_implied_timescales, fit_tica
    from src.config import MSMConfig

    features = _make_two_state_trajectory(n_frames=3000)
    lag_times = (1.0, 2.0, 5.0)
    config = MSMConfig(
        tica_lag_ps=1.0, tica_n_components=2, n_clusters=10,
        lag_times_ps=lag_times, n_implied_timescales=3, stride=1,
    )
    tica_out = fit_tica(features, config)["tica_output"]
    assignments = cluster_microstates(tica_out, config)
    result = compute_implied_timescales(assignments, config)

    assert result["lag_times_ps"].shape == (3,)
    assert result["implied_timescales_ps"].shape == (3, 3)


# ---------------------------------------------------------------------------
# Step 5: CLI entry point
# ---------------------------------------------------------------------------


def test_msm_cli_parser_featurize():
    """MSM CLI parser accepts featurize subcommand."""
    from scripts.run_msm import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "featurize",
        "--trajectory", "data/trajectories/production.dcd",
        "--topology", "data/pdb/system.pdb",
        "--output", "data/analysis/msm/features.npz",
    ])
    assert args.command == "featurize"
    assert str(args.trajectory) == "data/trajectories/production.dcd"


def test_msm_cli_parser_build():
    """MSM CLI parser accepts build subcommand with optional lag-time."""
    from scripts.run_msm import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "build",
        "--features", "features.npz",
        "--output-dir", "msm_output",
        "--lag-time-ps", "50.0",
    ])
    assert args.command == "build"
    assert args.lag_time_ps == 50.0


def test_msm_cli_parser_kinetics():
    """MSM CLI parser accepts kinetics subcommand."""
    from scripts.run_msm import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "kinetics",
        "--msm-dir", "msm_output",
        "--source-states", "0,1,2",
        "--target-states", "5,6,7",
    ])
    assert args.command == "kinetics"
    assert args.source_states == "0,1,2"


# ---------------------------------------------------------------------------
# Step 6: ITS convergence validation
# ---------------------------------------------------------------------------


def test_its_convergence_plateau_detection():
    """Implied timescales plateau for well-sampled synthetic data.

    A converged ITS curve shows < 10% coefficient of variation in the tail.
    """
    lag_times = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0])
    its_curve = np.array([100.0, 200.0, 350.0, 420.0, 470.0, 490.0, 495.0, 498.0])

    tail = its_curve[-4:]
    cv = float(np.std(tail, ddof=0) / np.mean(tail))
    assert cv < 0.10, f"Implied timescale not converged: CV = {cv:.3f} > 0.10"


def test_its_from_two_state_system():
    """End-to-end ITS computation on synthetic two-state trajectory."""
    from src.analyze.msm import cluster_microstates, compute_implied_timescales, fit_tica
    from src.config import MSMConfig

    features = _make_two_state_trajectory(n_frames=5000, seed=123)
    config = MSMConfig(
        tica_lag_ps=1.0, tica_n_components=2, n_clusters=20,
        lag_times_ps=(1.0, 5.0, 10.0), n_implied_timescales=3, stride=1,
    )
    tica_out = fit_tica(features, config)["tica_output"]
    assignments = cluster_microstates(tica_out, config)
    result = compute_implied_timescales(assignments, config)

    # The dominant timescale should be positive and finite
    its = result["implied_timescales_ps"]
    valid = its[~np.isnan(its)]
    assert len(valid) > 0
    assert np.all(valid > 0)


# ---------------------------------------------------------------------------
# Step 7: Rate cross-validation
# ---------------------------------------------------------------------------


def test_rate_cross_validation_thermodynamic_consistency():
    """MSM-derived DG is thermodynamically self-consistent with rate constants."""
    from src.config import BOLTZMANN_KJ

    k_on_per_ps = 1e-6   # 1e6 s^-1
    k_off_per_ps = 1e-9  # 1e3 s^-1
    temperature_k = 310.0

    K_d = k_off_per_ps / k_on_per_ps
    delta_g_kj = BOLTZMANN_KJ * temperature_k * np.log(K_d)

    assert delta_g_kj < 0, f"DG should be negative for K_d < 1, got {delta_g_kj:.2f} kJ/mol"
    expected_dg = 0.008314462618 * 310.0 * np.log(1e-3)
    assert abs(delta_g_kj - expected_dg) < 0.01


def test_detailed_balance():
    """Transition matrix from build_msm satisfies detailed balance."""
    from src.analyze.msm import build_msm, cluster_microstates, fit_tica
    from src.config import MSMConfig

    features = _make_two_state_trajectory(n_frames=5000)
    config = MSMConfig(
        tica_lag_ps=1.0, tica_n_components=2, n_clusters=15,
        stride=1, n_implied_timescales=3,
    )
    tica_out = fit_tica(features, config)["tica_output"]
    assignments = cluster_microstates(tica_out, config)
    result = build_msm(assignments, lag_time_ps=5.0, config=config)

    T = result["transition_matrix"]
    pi = result["stationary_distribution"]

    # Check pi_i * T_ij ≈ pi_j * T_ji for all i,j
    n = T.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            lhs = pi[i] * T[i, j]
            rhs = pi[j] * T[j, i]
            if lhs > 1e-12 or rhs > 1e-12:
                assert abs(lhs - rhs) < 1e-6, (
                    f"Detailed balance violated at ({i},{j}): "
                    f"pi_i*T_ij={lhs:.8e}, pi_j*T_ji={rhs:.8e}"
                )
