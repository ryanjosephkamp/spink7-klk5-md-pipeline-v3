"""Tests for the central configuration module."""

from dataclasses import FrozenInstanceError, is_dataclass
from pathlib import Path

import pytest

from src.config import (
    AVOGADRO,
    BOLTZMANN_KJ,
    KCAL_TO_KJ,
    DATA_DIR,
    PROJECT_ROOT,
    TRAJ_DIR,
    EquilibrationConfig,
    MinimizationConfig,
    ProductionConfig,
    SMDConfig,
    SystemConfig,
    UmbrellaConfig,
    WHAMConfig,
)


@pytest.mark.parametrize(
    "config_cls",
    [
        SystemConfig,
        MinimizationConfig,
        EquilibrationConfig,
        ProductionConfig,
        SMDConfig,
        UmbrellaConfig,
        WHAMConfig,
    ],
)
def test_config_dataclasses_instantiate(config_cls: type) -> None:
    """All configuration dataclasses should instantiate successfully."""

    instance = config_cls()

    assert is_dataclass(instance)


@pytest.mark.parametrize(
    ("config_instance", "field_name", "new_value"),
    [
        (SystemConfig(), "ph", 6.5),
        (MinimizationConfig(), "max_iterations", 5000),
        (EquilibrationConfig(), "temperature_k", 300.0),
        (ProductionConfig(), "duration_ns", 50.0),
        (SMDConfig(), "n_replicates", 10),
        (UmbrellaConfig(), "xi_min_nm", 1.0),
        (WHAMConfig(), "histogram_bins", 100),
    ],
)
def test_config_dataclasses_are_frozen(config_instance: object, field_name: str, new_value: object) -> None:
    """Chunk 1 gate: config dataclasses must be immutable."""

    with pytest.raises(FrozenInstanceError):
        setattr(config_instance, field_name, new_value)


def test_config_constants_match_blueprint_defaults() -> None:
    """Physical constants and default parameters should match the blueprint."""

    assert BOLTZMANN_KJ == pytest.approx(0.008314462618)
    assert AVOGADRO == pytest.approx(6.02214076e23)
    assert KCAL_TO_KJ == pytest.approx(4.184)

    assert SystemConfig().force_field == "amber14-all.xml"
    assert SystemConfig().water_model == "amber14/tip3p.xml"
    assert EquilibrationConfig().temperature_k == pytest.approx(310.0)
    assert ProductionConfig().checkpoint_interval_ps == pytest.approx(100.0)
    assert SMDConfig().pulling_velocity_nm_per_ps == pytest.approx(0.001)
    assert UmbrellaConfig().window_spacing_nm == pytest.approx(0.05)
    assert WHAMConfig().tolerance == pytest.approx(1e-6)


def test_project_paths_are_defined_relative_to_module_root() -> None:
    """Root paths should resolve relative to the project directory."""

    expected_root = Path(__file__).resolve().parent.parent

    assert PROJECT_ROOT == expected_root
    assert DATA_DIR == PROJECT_ROOT / "data"
    assert TRAJ_DIR == DATA_DIR / "trajectories"
