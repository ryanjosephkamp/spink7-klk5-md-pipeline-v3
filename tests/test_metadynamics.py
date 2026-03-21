"""Tests for metadynamics implementation (L-13)."""

from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Step 1: MetadynamicsConfig
# ---------------------------------------------------------------------------


def test_step_1_metadynamics_config_instantiates_and_is_frozen():
    """MetadynamicsConfig should instantiate with defaults and be immutable."""
    from src.config import MetadynamicsConfig

    config = MetadynamicsConfig()
    assert config.gaussian_height_kj_mol == 1.0
    assert config.bias_factor == 15.0
    assert config.temperature_k == 310.0
    assert config.grid_num_bins == 200
    assert config.gaussian_width_nm == 0.05
    assert config.deposition_interval_ps == 1.0
    assert config.simulation_duration_ns == 50.0
    assert config.save_interval_ps == 1.0
    assert config.convergence_tolerance_kj_mol == 1.0
    assert config.convergence_check_interval_ns == 1.0
    assert config.grid_min_nm == 1.0
    assert config.grid_max_nm == 5.0

    with pytest.raises(FrozenInstanceError):
        config.bias_factor = 10.0


# ---------------------------------------------------------------------------
# Step 2: Core module validation
# ---------------------------------------------------------------------------


def test_step_2_metadynamics_module_imports():
    """Metadynamics module should import cleanly."""
    from src.simulate.metadynamics import (
        _validate_config,
        _create_com_distance_bias_variable,
        _extract_fes,
        _check_fes_convergence,
        run_metadynamics,
    )
    assert callable(_validate_config)
    assert callable(_create_com_distance_bias_variable)
    assert callable(_extract_fes)
    assert callable(_check_fes_convergence)
    assert callable(run_metadynamics)


def test_step_2_validate_config_accepts_defaults():
    """Validation should accept default MetadynamicsConfig."""
    from src.config import MetadynamicsConfig
    from src.simulate.metadynamics import _validate_config

    _validate_config(MetadynamicsConfig())


def test_step_2_validate_config_rejects_invalid():
    """Validation should reject invalid parameter values."""
    from src.config import MetadynamicsConfig
    from src.simulate.metadynamics import _validate_config

    with pytest.raises(ValueError, match="gaussian_height_kj_mol must be positive"):
        _validate_config(MetadynamicsConfig(gaussian_height_kj_mol=-1.0))

    with pytest.raises(ValueError, match="bias_factor must be greater than 1.0"):
        _validate_config(MetadynamicsConfig(bias_factor=0.5))

    with pytest.raises(ValueError, match="grid_min_nm must be less than.*grid_max_nm"):
        _validate_config(MetadynamicsConfig(grid_min_nm=5.0, grid_max_nm=1.0))

    with pytest.raises(ValueError, match="grid_num_bins must be at least 2"):
        _validate_config(MetadynamicsConfig(grid_num_bins=1))


# ---------------------------------------------------------------------------
# Step 3: Convergence diagnostics
# ---------------------------------------------------------------------------


def test_step_3_fes_convergence_check_detects_convergence():
    """Convergence check should return True when FES difference is below tolerance."""
    from src.simulate.metadynamics import _check_fes_convergence

    fes_a = np.array([0.0, -5.0, -10.0, -8.0, -3.0])
    fes_b = np.array([0.0, -5.1, -10.2, -8.1, -3.0])

    # Max deviation is 0.2 kJ/mol — should converge at tolerance 0.5
    converged, max_delta = _check_fes_convergence(fes_a, fes_b, tolerance_kj_mol=0.5)
    assert converged is True
    assert max_delta == pytest.approx(0.2, abs=1e-10)

    # Should NOT converge at tolerance 0.1
    converged, max_delta = _check_fes_convergence(fes_a, fes_b, tolerance_kj_mol=0.1)
    assert converged is False


def test_step_3_fes_convergence_identical():
    """Identical FES should yield zero deviation."""
    from src.simulate.metadynamics import _check_fes_convergence

    fes = np.array([0.0, -5.0, -10.0, -5.0, 0.0])
    converged, max_delta = _check_fes_convergence(fes, fes, tolerance_kj_mol=0.001)
    assert converged is True
    assert max_delta == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# Step 4: Integration test with alanine dipeptide simulation
# ---------------------------------------------------------------------------


def test_step_4_metadynamics_runs_on_alanine_dipeptide(
    tmp_path, alanine_dipeptide_simulation,
):
    """Metadynamics should run without error on a minimal solvated system."""
    from src.config import MetadynamicsConfig
    from src.simulate.metadynamics import run_metadynamics

    simulation = alanine_dipeptide_simulation

    # Infer two groups: solute heavy atoms vs. a small water subset
    # For a single-chain system we split solute into two halves as
    # a minimal exercise of the COM-distance CV.
    solute_indices = []
    for chain in simulation.topology.chains():
        for atom in chain.atoms():
            solute_indices.append(atom.index)
        break  # only first chain

    mid = max(1, len(solute_indices) // 2)
    group_a = solute_indices[:mid]
    group_b = solute_indices[mid:]

    config = MetadynamicsConfig(
        simulation_duration_ns=0.001,   # 1 ps — extremely short for CI
        deposition_interval_ps=0.1,
        save_interval_ps=0.1,
        convergence_check_interval_ns=0.001,
        convergence_tolerance_kj_mol=50.0,  # relaxed for short run
        grid_min_nm=0.0,
        grid_max_nm=2.0,
        grid_num_bins=50,
        gaussian_width_nm=0.05,
        gaussian_height_kj_mol=0.5,
        bias_factor=5.0,
    )

    result = run_metadynamics(
        simulation, config, tmp_path / "metad_output",
        pull_group_1=group_a,
        pull_group_2=group_b,
    )

    assert result["fes_grid_nm"] is not None
    assert result["fes_kj_mol"] is not None
    assert len(result["xi_timeseries"]) > 0
    assert result["output_dir"].exists()
    assert (result["output_dir"] / "metadynamics_fes.npy").exists()
    assert (result["output_dir"] / "metadynamics_xi_timeseries.npy").exists()


# ---------------------------------------------------------------------------
# Step 5: Cross-validation utility (compare_fes_profiles)
# ---------------------------------------------------------------------------


def test_step_5_compare_fes_profiles_identical_curves():
    """Identical FES profiles should yield zero deviation and perfect correlation."""
    from src.analyze.convergence import compare_fes_profiles

    grid = np.linspace(1.0, 4.0, 100)
    fes = -10.0 * np.exp(-((grid - 2.5) ** 2) / 0.5)

    result = compare_fes_profiles(grid, fes, grid, fes)
    assert result["max_absolute_deviation_kj_mol"] == pytest.approx(0.0, abs=1e-10)
    assert result["rmsd_kj_mol"] == pytest.approx(0.0, abs=1e-10)
    assert result["pearson_r"] == pytest.approx(1.0, abs=1e-10)


def test_step_5_compare_fes_profiles_shifted_curves():
    """Profiles with constant offset should agree after alignment."""
    from src.analyze.convergence import compare_fes_profiles

    grid = np.linspace(1.0, 4.0, 100)
    fes_a = -10.0 * np.exp(-((grid - 2.5) ** 2) / 0.5)
    fes_b = fes_a + 5.0  # constant offset

    result = compare_fes_profiles(grid, fes_a, grid, fes_b)
    # After alignment (set G(xi_max)=0), curves should be identical
    assert result["max_absolute_deviation_kj_mol"] == pytest.approx(0.0, abs=1e-8)
    assert result["pearson_r"] == pytest.approx(1.0, abs=1e-10)


def test_step_5_compare_fes_profiles_different_grids():
    """Profiles on different grids should be interpolated and compared."""
    from src.analyze.convergence import compare_fes_profiles

    grid_a = np.linspace(1.0, 4.0, 50)
    grid_b = np.linspace(1.0, 4.0, 100)
    fes_a = -10.0 * np.exp(-((grid_a - 2.5) ** 2) / 0.5)
    fes_b = -10.0 * np.exp(-((grid_b - 2.5) ** 2) / 0.5)

    result = compare_fes_profiles(grid_a, fes_a, grid_b, fes_b)
    assert result["max_absolute_deviation_kj_mol"] < 0.5  # interpolation error
    assert result["pearson_r"] > 0.99
