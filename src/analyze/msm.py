"""Markov State Model construction and analysis.

This module implements the MSM pipeline: TICA dimensionality reduction,
k-means clustering, transition matrix estimation, implied timescale
convergence analysis, and kinetic rate extraction.
"""

from __future__ import annotations

import logging

import numpy as np

from src.config import MSMConfig

__all__ = [
    "fit_tica",
    "cluster_microstates",
    "compute_implied_timescales",
    "build_msm",
    "compute_mfpt",
]

logger = logging.getLogger(__name__)


def fit_tica(
    features: np.ndarray,
    config: MSMConfig,
) -> dict[str, np.ndarray]:
    """Fit TICA and transform features to the slow subspace.

    Args:
        features: Feature matrix. Shape: [N_frames, N_features].
        config: MSM configuration parameters.

    Returns:
        dict with keys:
            "tica_output": Transformed data. Shape: [N_frames, tica_n_components].
            "eigenvalues": Singular values. Shape: [tica_n_components].
    """
    from deeptime.decomposition import TICA

    tica_lag_steps = max(1, round(config.tica_lag_ps / max(config.stride, 1)))
    n_components = min(config.tica_n_components, features.shape[1])
    tica = TICA(lagtime=tica_lag_steps, dim=n_components)
    model = tica.fit_fetch(features)
    transformed = model.transform(features)

    n_sv = min(n_components, len(model.singular_values))
    singular_values = model.singular_values[:n_sv]

    logger.info(
        "TICA: %d features -> %d components, explained kinetic variance: %.3f",
        features.shape[1],
        n_components,
        float(np.sum(singular_values**2)),
    )

    return {
        "tica_output": np.asarray(transformed, dtype=float),
        "eigenvalues": np.asarray(singular_values, dtype=float),
    }


def cluster_microstates(
    tica_output: np.ndarray,
    config: MSMConfig,
) -> np.ndarray:
    """Cluster TICA-transformed features into microstates via k-means.

    Args:
        tica_output: TICA-transformed features. Shape: [N_frames, N_tica].
        config: MSM configuration parameters.

    Returns:
        np.ndarray: Microstate assignments. Shape: [N_frames]. dtype: int.
    """
    from deeptime.clustering import KMeans

    n_clusters = min(config.n_clusters, tica_output.shape[0])
    kmeans = KMeans(
        n_clusters=n_clusters,
        max_iter=300,
        init_strategy="kmeans++",
    )
    model = kmeans.fit_fetch(tica_output)
    assignments = model.transform(tica_output)

    logger.info(
        "k-Means: %d frames -> %d clusters, inertia: %.4f",
        tica_output.shape[0],
        n_clusters,
        float(model.inertia),
    )

    return np.asarray(assignments, dtype=int)


def compute_implied_timescales(
    assignments: np.ndarray,
    config: MSMConfig,
) -> dict[str, np.ndarray]:
    """Compute implied timescales at multiple lag times for convergence analysis.

    For a well-constructed MSM, implied timescales should plateau as lag time
    increases, indicating that the Markov assumption is satisfied.

    Args:
        assignments: Microstate assignments. Shape: [N_frames].
        config: MSM configuration parameters.

    Returns:
        dict with keys:
            "lag_times_ps": Lag times used. Shape: [N_lags].
            "implied_timescales_ps": ITS matrix. Shape: [N_lags, n_implied_timescales].
    """
    from deeptime.markov.msm import MaximumLikelihoodMSM

    lag_times_ps = np.array(config.lag_times_ps, dtype=float)
    n_its = config.n_implied_timescales
    its_matrix = np.full((len(lag_times_ps), n_its), np.nan, dtype=float)

    for i, lag_ps in enumerate(lag_times_ps):
        lag_steps = max(1, round(lag_ps / max(config.stride, 1)))
        try:
            estimator = MaximumLikelihoodMSM(lagtime=lag_steps)
            msm = estimator.fit_fetch(assignments)
            timescales = msm.timescales(k=n_its)
            n_available = min(len(timescales), n_its)
            its_matrix[i, :n_available] = timescales[:n_available] * max(config.stride, 1)
        except Exception as exc:
            logger.warning("MSM estimation failed at lag=%.1f ps: %s", lag_ps, exc)

    return {
        "lag_times_ps": lag_times_ps,
        "implied_timescales_ps": its_matrix,
    }


def build_msm(
    assignments: np.ndarray,
    lag_time_ps: float,
    config: MSMConfig,
) -> dict[str, object]:
    """Build a maximum-likelihood MSM at the selected lag time.

    Args:
        assignments: Microstate assignments. Shape: [N_frames].
        lag_time_ps: Selected lag time in ps (from implied timescale convergence).
        config: MSM configuration parameters.

    Returns:
        dict with keys:
            "transition_matrix": Row-stochastic matrix. Shape: [N_active, N_active].
            "stationary_distribution": Equilibrium probabilities. Shape: [N_active].
            "timescales_ps": Relaxation timescales. Shape: [n_its].
            "active_set": Active state indices. Shape: [N_active].
            "msm_model": The fitted deeptime MSM model.
    """
    from deeptime.markov.msm import MaximumLikelihoodMSM

    lag_steps = max(1, round(lag_time_ps / max(config.stride, 1)))
    estimator = MaximumLikelihoodMSM(lagtime=lag_steps)
    msm = estimator.fit_fetch(assignments)

    logger.info(
        "MSM built: lag=%.1f ps, %d active states, %d timescales",
        lag_time_ps,
        msm.n_states,
        len(msm.timescales()),
    )

    return {
        "transition_matrix": np.asarray(msm.transition_matrix, dtype=float),
        "stationary_distribution": np.asarray(msm.stationary_distribution, dtype=float),
        "timescales_ps": np.asarray(msm.timescales() * max(config.stride, 1), dtype=float),
        "active_set": np.asarray(msm.count_model.state_symbols, dtype=int),
        "msm_model": msm,
    }


def compute_mfpt(
    msm_model: object,
    source_states: np.ndarray,
    target_states: np.ndarray,
    stride: float = 1.0,
) -> float:
    """Compute the mean first passage time from source to target macrostates.

    Args:
        msm_model: Fitted deeptime MSM model.
        source_states: Microstate indices belonging to the source macrostate.
        target_states: Microstate indices belonging to the target macrostate.
        stride: Timestep between frames in ps, for converting to physical time.

    Returns:
        float: MFPT in ps.
    """
    mfpt_steps = msm_model.mfpt(source_states, target_states)
    return float(mfpt_steps) * stride
