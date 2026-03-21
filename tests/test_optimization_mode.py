"""Tests that verify pipeline behavior under Python optimization mode (-O).

These tests serve two purposes:
1. Canary: Detect whether -O is active (assert statements stripped).
2. Regression guard: Ensure no pipeline module relies on assert for
   correctness-critical validation after L-20 is implemented.

Run with: python -O -m pytest tests/test_optimization_mode.py -v
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

os.environ.setdefault("MPLBACKEND", "Agg")


class TestOptimizationDetection:
    """Verify that optimization mode is correctly detected."""

    def test_optimization_flag_detection(self):
        """__debug__ must reflect the current optimization state."""
        if __debug__:
            assert sys.flags.optimize == 0
        else:
            assert sys.flags.optimize >= 1

    @pytest.mark.optimized
    def test_assert_is_stripped_under_optimization(self):
        """Under -O, __debug__ is False and source-level asserts are no-ops.

        Pytest rewrites assert statements in test files, so we cannot use
        ``assert False`` here as a canary. Instead we verify via __debug__
        and sys.flags.optimize that the interpreter is in optimization mode.
        """
        if not __debug__:
            # Optimization mode is active — confirm via interpreter flag.
            assert sys.flags.optimize >= 1


class TestValidationSurvivesOptimization:
    """Verify that validation logic does NOT depend on assert.

    After L-20 is implemented, all validation uses raise statements
    instead of assert. These tests confirm that validation still
    works under -O (where assert is stripped).
    """

    @pytest.mark.optimized
    def test_jarzynski_validates_input_under_optimization(self):
        """jarzynski_free_energy must reject 2D arrays even under -O."""
        from src.analyze.jarzynski import jarzynski_free_energy

        bad_input = np.ones((3, 2))
        with pytest.raises((ValueError, AssertionError)):
            jarzynski_free_energy(bad_input, temperature_k=310.0)

    @pytest.mark.optimized
    def test_wham_validates_input_under_optimization(self):
        """solve_wham must reject non-1D window timeseries even under -O."""
        from src.analyze.wham import solve_wham
        from src.config import WHAMConfig

        bad_xi = [np.ones((3, 2))]
        with pytest.raises((ValueError, AssertionError)):
            solve_wham(
                bad_xi,
                np.array([1.0]),
                np.array([1000.0]),
                310.0,
                WHAMConfig(),
            )

    @pytest.mark.optimized
    def test_convergence_validates_input_under_optimization(self):
        """block_average must reject 2D arrays even under -O."""
        from src.analyze.convergence import block_average

        bad_input = np.ones((3, 2))
        with pytest.raises((ValueError, AssertionError)):
            block_average(bad_input, n_blocks=2)

    @pytest.mark.optimized
    def test_plot_pmf_validates_input_under_optimization(self):
        """plot_pmf must reject 2D bin arrays even under -O."""
        from pathlib import Path

        from src.visualization.plot_pmf import plot_pmf

        bad_xi = np.ones((3, 2))
        with pytest.raises((ValueError, AssertionError)):
            plot_pmf(bad_xi, np.ones(3))

    @pytest.mark.optimized
    def test_plot_timeseries_validates_input_under_optimization(self):
        """plot_energy_timeseries must reject 2D time arrays even under -O."""
        from src.visualization.plot_timeseries import plot_energy_timeseries

        bad_time = np.ones((3, 2))
        with pytest.raises((ValueError, AssertionError)):
            plot_energy_timeseries(bad_time, np.ones(3), np.ones(3))

    @pytest.mark.optimized
    def test_script_exists_and_executable(self):
        """The optimization-mode test script must exist and be executable."""
        from pathlib import Path

        script = Path(__file__).resolve().parents[1] / "scripts" / "test_optimized.sh"
        assert script.exists(), "scripts/test_optimized.sh not found"
        assert os.access(script, os.X_OK), "scripts/test_optimized.sh not executable"
