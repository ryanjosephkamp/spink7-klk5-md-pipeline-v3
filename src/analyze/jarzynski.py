"""Jarzynski free-energy estimators for nonequilibrium SMD work data."""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import brentq
from scipy.stats import anderson

from src.config import BOLTZMANN_KJ, KCAL_TO_KJ


def _validate_work_values(work_values: np.ndarray) -> np.ndarray:
    """Validate a one-dimensional work array at the public API boundary."""

    values = np.asarray(work_values, dtype=float)
    assert values.ndim == 1, "work_values must have shape [N_replicates]"
    if values.ndim != 1 or values.size == 0:
        raise ValueError("work_values must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("work_values must contain only finite values")
    return values


def _validate_temperature(temperature_k: float) -> float:
    """Validate the thermodynamic temperature input."""

    temperature = float(temperature_k)
    if temperature <= 0.0:
        raise ValueError("temperature_k must be positive")
    return temperature


def _beta(temperature_k: float) -> float:
    """Return the inverse thermal energy in mol/kJ."""

    return 1.0 / (BOLTZMANN_KJ * temperature_k)


def _test_work_normality(
    work_values: np.ndarray,
    significance_level: float = 0.05,
) -> tuple[bool, float]:
    """Test whether the work distribution is approximately Gaussian.

    Uses the Anderson-Darling test, which is valid for large sample sizes
    (N > 5000) unlike Shapiro-Wilk.

    Returns:
        (is_normal, statistic): Whether normality is accepted and the AD statistic.
    """

    result = anderson(work_values, dist="norm", method="interpolate")
    is_normal = bool(result.pvalue >= significance_level)
    return is_normal, float(result.statistic)


def _jarzynski_exact_kj_mol(work_values: np.ndarray, temperature_k: float) -> float:
    """Compute the numerically stable exact Jarzynski estimator in kJ/mol."""

    # Numerically stable Jarzynski: shift by max(-beta * W) before
    # exponentiation to avoid overflow in the exponential average.
    beta = _beta(temperature_k)
    scaled_work = -beta * work_values
    maximum = np.max(scaled_work)
    log_average = maximum + np.log(np.mean(np.exp(scaled_work - maximum), dtype=float))
    return -(1.0 / beta) * log_average


def jarzynski_free_energy(
    work_values: np.ndarray,
    temperature_k: float,
) -> dict[str, float]:
    """Estimate free-energy differences from nonequilibrium work values.

    Args:
        work_values: Total work values in kJ/mol. Shape: [N_replicates].
        temperature_k: Temperature in Kelvin.

    Returns:
        dict[str, float]: Exact Jarzynski and second-order cumulant estimates.
    """

    values = _validate_work_values(work_values)
    temperature = _validate_temperature(temperature_k)
    beta = _beta(temperature)

    mean_work_kj_mol = float(np.mean(values, dtype=float))
    work_variance_kj_mol2 = float(np.var(values, ddof=0, dtype=float))
    delta_g_kj_mol = float(_jarzynski_exact_kj_mol(values, temperature))
    delta_g_cumulant2_kj_mol = float(mean_work_kj_mol - 0.5 * beta * work_variance_kj_mol2)

    # Dissipation diagnostics (backward-compatible additions).
    w_diss_kj_mol = float(mean_work_kj_mol - delta_g_kj_mol)
    w_diss_over_kbt = float(beta * w_diss_kj_mol)
    efficiency_ratio = float(np.exp(-beta * w_diss_kj_mol))

    # Normality test for cumulant expansion validity.
    normality_passed, ad_statistic = _test_work_normality(values)
    if not normality_passed:
        warnings.warn(
            "Work distribution failed the Anderson-Darling normality test. "
            "The second-order cumulant expansion is unreliable; "
            "use the BAR estimator (bar_free_energy) for robust results.",
            UserWarning,
            stacklevel=2,
        )

    return {
        "delta_g_kj_mol": delta_g_kj_mol,
        "delta_g_kcal_mol": delta_g_kj_mol / KCAL_TO_KJ,
        "mean_work_kj_mol": mean_work_kj_mol,
        "work_variance_kj_mol2": work_variance_kj_mol2,
        "delta_g_cumulant2_kj_mol": delta_g_cumulant2_kj_mol,
        "delta_g_cumulant2_kcal_mol": delta_g_cumulant2_kj_mol / KCAL_TO_KJ,
        "w_diss_kj_mol": w_diss_kj_mol,
        "w_diss_over_kbt": w_diss_over_kbt,
        "efficiency_ratio": efficiency_ratio,
        "cumulant2_normality_passed": normality_passed,
        "anderson_darling_statistic": ad_statistic,
    }


def diagnose_dissipation(
    work_values: np.ndarray,
    temperature_k: float,
    threshold_kbt: float = 3.0,
) -> dict[str, float]:
    """Compute dissipation diagnostics and warn when estimator is unreliable.

    Args:
        work_values: Total work values in kJ/mol. Shape: [N_replicates].
        temperature_k: Temperature in Kelvin.
        threshold_kbt: Warn when W_diss exceeds this many kBT units.

    Returns:
        dict with w_diss_kj_mol, w_diss_over_kbt, efficiency_ratio,
        and warning_triggered.
    """

    values = _validate_work_values(work_values)
    temperature = _validate_temperature(temperature_k)
    beta_val = _beta(temperature)

    delta_g = _jarzynski_exact_kj_mol(values, temperature)
    mean_work = float(np.mean(values, dtype=float))
    w_diss = mean_work - delta_g
    w_diss_over_kbt = beta_val * w_diss
    efficiency_ratio = float(np.exp(-beta_val * w_diss))

    warning_triggered = w_diss_over_kbt > threshold_kbt
    if warning_triggered:
        warnings.warn(
            f"High dissipation detected: W_diss / kBT = {w_diss_over_kbt:.1f} "
            f"(threshold = {threshold_kbt:.1f}). The Jarzynski exponential "
            f"average requires exponentially more trajectories for convergence. "
            f"Consider increasing N_traj or decreasing pulling velocity.",
            UserWarning,
            stacklevel=2,
        )

    return {
        "w_diss_kj_mol": float(w_diss),
        "w_diss_over_kbt": float(w_diss_over_kbt),
        "efficiency_ratio": efficiency_ratio,
        "warning_triggered": warning_triggered,
    }


def bar_free_energy(
    forward_work: np.ndarray,
    reverse_work: np.ndarray,
    temperature_k: float,
    max_iterations: int = 1000,
    tolerance_kj_mol: float = 1e-6,
) -> dict[str, float]:
    """Compute the Bennett Acceptance Ratio (BAR) free energy estimate.

    Uses forward and reverse nonequilibrium work distributions to solve the
    BAR equation via Brent's method, yielding the statistically optimal
    free energy estimate from the Crooks Fluctuation Theorem.

    Args:
        forward_work: Forward (pulling) work values in kJ/mol.
        reverse_work: Reverse (reinsertion) work values in kJ/mol.
        temperature_k: Temperature in Kelvin.
        max_iterations: Maximum Brent solver iterations.
        tolerance_kj_mol: Convergence tolerance in kJ/mol.

    Returns:
        dict with BAR free energy estimate and convergence info.
    """

    fwd = _validate_work_values(forward_work)
    rev = _validate_work_values(reverse_work)
    temperature = _validate_temperature(temperature_k)
    beta_val = _beta(temperature)
    n_f = fwd.size
    n_r = rev.size
    ln_ratio_fr = np.log(n_f / n_r)

    def _bar_objective(delta_g: float) -> float:
        # BAR equation: sum of Fermi functions must balance.
        # Fermi: sigma(-x) = 1/(1+exp(x)), computed with clipping for stability.
        # Forward: x_i = ln(N_F/N_R) + beta*(W_i^F - dG)
        x_fwd = ln_ratio_fr + beta_val * (fwd - delta_g)
        lhs = float(np.sum(1.0 / (1.0 + np.exp(np.clip(x_fwd, -500.0, 500.0)))))

        # Reverse: x_j = ln(N_R/N_F) + beta*(W_j^R + dG)
        x_rev = -ln_ratio_fr + beta_val * (rev + delta_g)
        rhs = float(np.sum(1.0 / (1.0 + np.exp(np.clip(x_rev, -500.0, 500.0)))))

        return lhs - rhs

    # Bracket the root: use a wide bracket around the cumulant expansion estimate.
    bracket_half = 50.0 / beta_val
    dg_init = float(np.mean(fwd)) - 0.5 * beta_val * float(np.var(fwd, ddof=0))
    lo = dg_init - bracket_half
    hi = dg_init + bracket_half

    try:
        delta_g_bar = brentq(
            _bar_objective, lo, hi,
            xtol=tolerance_kj_mol, maxiter=max_iterations,
        )
        converged = True
    except ValueError:
        # Bracket does not contain a root; fall back to cumulant estimate.
        delta_g_bar = dg_init
        converged = False

    return {
        "delta_g_bar_kj_mol": float(delta_g_bar),
        "delta_g_bar_kcal_mol": float(delta_g_bar) / KCAL_TO_KJ,
        "n_forward": n_f,
        "n_reverse": n_r,
        "converged": converged,
    }


def evaluate_convergence(
    work_values: np.ndarray,
    temperature_k: float,
    n_subsets: int = 10,
) -> dict[str, np.ndarray]:
    """Evaluate Jarzynski estimator convergence as a function of trajectory count.

    Args:
        work_values: Total work values in kJ/mol. Shape: [N_replicates].
        temperature_k: Temperature in Kelvin.
        n_subsets: Number of prefix subset sizes to evaluate.

    Returns:
        dict[str, np.ndarray]: Prefix sizes, exact Jarzynski estimates, and bootstrap SE.
    """

    values = _validate_work_values(work_values)
    temperature = _validate_temperature(temperature_k)

    subset_count = int(n_subsets)
    if subset_count <= 0:
        raise ValueError("n_subsets must be positive")

    # Evaluate convergence by computing the Jarzynski estimate for
    # growing prefix subsets and estimating bootstrap standard errors.
    subset_sizes = np.unique(np.linspace(1, values.size, subset_count, dtype=int))
    delta_g_vs_n = np.empty(subset_sizes.size, dtype=float)
    std_vs_n = np.empty(subset_sizes.size, dtype=float)
    rng = np.random.default_rng(0)
    n_bootstrap = 200

    for subset_index, subset_size in enumerate(subset_sizes):
        subset_values = values[:subset_size]
        delta_g_vs_n[subset_index] = _jarzynski_exact_kj_mol(subset_values, temperature)

        bootstrap_estimates = np.empty(n_bootstrap, dtype=float)
        for bootstrap_index in range(n_bootstrap):
            sampled_indices = rng.integers(0, subset_size, size=subset_size)
            bootstrap_estimates[bootstrap_index] = _jarzynski_exact_kj_mol(subset_values[sampled_indices], temperature)
        std_vs_n[subset_index] = np.std(bootstrap_estimates, ddof=1, dtype=float)

    return {
        "subset_sizes": subset_sizes,
        "delta_g_vs_n": delta_g_vs_n,
        "std_vs_n": std_vs_n,
    }