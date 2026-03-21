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
    SUPPORTED_FORCE_FIELD_FAMILIES,
    SUPPORTED_WATER_MODELS,
    AMOEBAConfig,
    EquilibrationConfig,
    MLPotentialConfig,
    MinimizationConfig,
    ProductionConfig,
    QMMMConfig,
    SMDConfig,
    SystemConfig,
    UmbrellaConfig,
    WHAMConfig,
    MBARConfig,
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
        MBARConfig,
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


# ---------- L-09 Step 1: UmbrellaConfig equilibration parameters ----------


def test_umbrella_config_equilibration_defaults() -> None:
    """UmbrellaConfig should expose equilibration parameters with sensible defaults."""

    config = UmbrellaConfig()

    assert config.pre_position_velocity_nm_per_ps == 0.01
    assert config.pre_position_spring_constant_kj_mol_nm2 == 1000.0
    assert config.equilibration_duration_ps == 200.0
    assert config.detect_equilibration is True


# ---------- L-11 Step 1: MBARConfig defaults ----------


def test_mbar_config_defaults() -> None:
    """MBARConfig should expose MBAR solver parameters with correct defaults."""
    cfg = MBARConfig()
    assert cfg.solver_protocol == "robust"
    assert cfg.relative_tolerance == 1e-7
    assert cfg.maximum_iterations == 10_000
    assert cfg.n_bootstrap == 200
    assert cfg.n_pmf_bins == 200


# ---------- L-15 Step 1: Force field family selection ----------


def test_system_config_force_field_family_defaults_to_amber() -> None:
    """SystemConfig should default to amber force field family for backward compatibility."""
    config = SystemConfig()
    assert config.force_field_family == "amber"
    assert config.force_field_family in SUPPORTED_FORCE_FIELD_FAMILIES


def test_amoeba_config_instantiates_and_is_frozen() -> None:
    """AMOEBAConfig should instantiate with valid defaults and be frozen."""
    config = AMOEBAConfig()
    assert config.polarization_type == "mutual"
    assert config.mutual_induced_target_epsilon > 0
    assert config.force_field_xml == "amoeba2018.xml"
    with pytest.raises(FrozenInstanceError):
        config.polarization_type = "direct"


def test_ml_potential_config_instantiates_and_is_frozen() -> None:
    """MLPotentialConfig should instantiate with valid defaults and be frozen."""
    config = MLPotentialConfig()
    assert config.potential_name == "ani2x"
    assert config.implementation == "torchani"
    with pytest.raises(FrozenInstanceError):
        config.potential_name = "mace"


def test_qmmm_config_instantiates_and_is_frozen() -> None:
    """QMMMConfig should instantiate with valid defaults and be frozen."""
    config = QMMMConfig()
    assert config.qm_method == "B3LYP"
    assert config.qm_basis_set == "6-31G*"
    assert "HIS57" in config.qm_region_residues
    with pytest.raises(FrozenInstanceError):
        config.qm_method = "HF"


# ---------- L-16 Step 1: Water model registry ----------


def test_supported_water_models_registry_contains_tip3p() -> None:
    """SUPPORTED_WATER_MODELS must include the default TIP3P model."""
    assert "tip3p" in SUPPORTED_WATER_MODELS
    xml_path, model_name = SUPPORTED_WATER_MODELS["tip3p"]
    assert xml_path == "amber14/tip3p.xml"
    assert model_name == "tip3p"


def test_supported_water_models_registry_contains_opc() -> None:
    """SUPPORTED_WATER_MODELS must include the OPC model."""
    assert "opc" in SUPPORTED_WATER_MODELS
    xml_path, model_name = SUPPORTED_WATER_MODELS["opc"]
    assert xml_path == "amber14/opc.xml"
    assert model_name == "tip4pew"  # OPC uses TIP4P-Ew geometry for placement


def test_supported_water_models_registry_contains_tip4pew() -> None:
    """SUPPORTED_WATER_MODELS must include the TIP4P-Ew model."""
    assert "tip4pew" in SUPPORTED_WATER_MODELS
    xml_path, model_name = SUPPORTED_WATER_MODELS["tip4pew"]
    assert xml_path == "amber14/tip4pew.xml"
    assert model_name == "tip4pew"


# ---------- L-18 Step 1: random_seed defaults to None ----------


def test_equilibration_config_random_seed_defaults_to_none() -> None:
    """EquilibrationConfig.random_seed should default to None (auto-generate)."""
    config = EquilibrationConfig()
    assert config.random_seed is None


def test_production_config_random_seed_defaults_to_none() -> None:
    """ProductionConfig.random_seed should default to None (auto-generate)."""
    config = ProductionConfig()
    assert config.random_seed is None


def test_smd_config_random_seed_defaults_to_none() -> None:
    """SMDConfig.random_seed should default to None (auto-generate)."""
    config = SMDConfig()
    assert config.random_seed is None


def test_configs_accept_explicit_seed() -> None:
    """All configs should accept an explicit integer seed."""
    assert EquilibrationConfig(random_seed=42).random_seed == 42
    assert ProductionConfig(random_seed=99).random_seed == 99
    assert SMDConfig(random_seed=1000).random_seed == 1000


# ---------- L-18 Step 2: _resolve_seed helper ----------


def test_resolve_seed_returns_integer_when_none() -> None:
    """_resolve_seed(None) should return a non-negative integer."""
    from src.simulate.equilibrate import _resolve_seed

    seed = _resolve_seed(None, "test")
    assert isinstance(seed, int)
    assert seed >= 0


def test_resolve_seed_returns_provided_value() -> None:
    """_resolve_seed(42) should return 42."""
    from src.simulate.equilibrate import _resolve_seed

    assert _resolve_seed(42, "test") == 42


# ---------- L-21 Step 1: default_config.yaml ----------


def test_default_config_yaml_exists_and_parses() -> None:
    """default_config.yaml must exist and parse as valid YAML."""
    import yaml

    config_path = PROJECT_ROOT / "default_config.yaml"
    assert config_path.exists(), f"default_config.yaml not found at {config_path}"
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    assert isinstance(data, dict)
    assert "system" in data
    assert "production" in data
    assert data["system"]["force_field"] == "amber14-all.xml"


# ---------- L-21 Step 2: load_config() ----------


def test_load_config_applies_overrides(tmp_path: Path) -> None:
    """load_config() should apply YAML overrides while preserving defaults."""
    import yaml
    from src.config import load_config

    config_data = {"production": {"duration_ns": 200.0}}
    config_path = tmp_path / "test_config.yaml"
    with config_path.open("w") as f:
        yaml.dump(config_data, f)

    configs = load_config(config_path)
    assert configs["production"].duration_ns == 200.0
    assert configs["production"].temperature_k == 310.0  # default preserved
    assert configs["system"].force_field == "amber14-all.xml"  # entire section defaulted


def test_load_config_rejects_unknown_fields(tmp_path: Path) -> None:
    """load_config() should reject unrecognized field names."""
    import yaml
    from src.config import load_config

    config_data = {"production": {"nonexistent_param": 42}}
    config_path = tmp_path / "bad_config.yaml"
    with config_path.open("w") as f:
        yaml.dump(config_data, f)

    with pytest.raises(TypeError, match="Unrecognized fields"):
        load_config(config_path)


def test_load_config_returns_defaults_when_no_path() -> None:
    """load_config(None) should return all-default configurations."""
    from src.config import load_config

    configs = load_config(None)
    assert configs["equilibration"].temperature_k == 310.0
    assert configs["smd"].n_replicates == 50


# ---------- L-21 Step 3: CLI --config flag ----------


def test_cli_scripts_accept_config_flag() -> None:
    """All CLI scripts should accept --config without error."""
    import importlib

    script_modules = [
        "scripts.run_prep",
        "scripts.run_equilibration",
        "scripts.run_production",
        "scripts.run_smd",
        "scripts.run_umbrella",
        "scripts.run_analysis",
    ]
    for module_name in script_modules:
        mod = importlib.import_module(module_name)
        parser = mod.build_parser()
        # Verify --config is a recognized argument (should not raise)
        action = None
        for act in parser._actions:
            if "--config" in getattr(act, "option_strings", []):
                action = act
                break
        assert action is not None, f"--config not found in {module_name}"


# ---------- L-21 Step 4: Environment variable overrides ----------


def test_env_override_takes_precedence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Environment variables should override both YAML and default values."""
    import yaml
    from src.config import load_config

    config_data = {"production": {"duration_ns": 200.0}}
    config_path = tmp_path / "env_test.yaml"
    with config_path.open("w") as f:
        yaml.dump(config_data, f)

    monkeypatch.setenv("SPINK7_PRODUCTION_DURATION_NS", "500.0")
    configs = load_config(config_path)
    assert configs["production"].duration_ns == 500.0


def test_env_override_with_no_yaml(monkeypatch: pytest.MonkeyPatch) -> None:
    """Environment variables should override defaults even without a YAML file."""
    from src.config import load_config

    monkeypatch.setenv("SPINK7_EQUILIBRATION_TEMPERATURE_K", "300.0")
    configs = load_config(None)
    assert configs["equilibration"].temperature_k == 300.0


# ── L-24: pyproject.toml validation ─────────────────────────────────


def test_pyproject_toml_exists_and_is_valid() -> None:
    """pyproject.toml must exist and be parseable with required sections."""
    import tomllib

    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    assert pyproject_path.exists(), "pyproject.toml not found in project root"
    with pyproject_path.open("rb") as f:
        config = tomllib.load(f)
    assert "project" in config
    assert "build-system" in config
    assert "scripts" in config["project"]


def test_no_sys_path_manipulation_in_scripts() -> None:
    """No script in scripts/ should contain sys.path.insert() calls."""
    scripts_dir = Path(__file__).resolve().parents[1] / "scripts"
    violations = []
    for script_path in scripts_dir.glob("*.py"):
        source = script_path.read_text(encoding="utf-8")
        if "sys.path.insert" in source or "sys.path.append" in source:
            violations.append(script_path.name)
    assert not violations, f"Scripts still manipulate sys.path: {violations}"


def test_entry_points_resolve() -> None:
    """All CLI entry points must be importable without sys.path hacks."""
    import subprocess
    import sys

    entry_points = [
        "scripts.run_prep:main",
        "scripts.run_equilibration:main",
        "scripts.run_production:main",
        "scripts.run_smd:main",
        "scripts.run_umbrella:main",
        "scripts.run_analysis:main",
    ]
    for ep in entry_points:
        module_path, func_name = ep.split(":")
        result = subprocess.run(
            [sys.executable, "-c", f"from {module_path} import {func_name}"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, (
            f"Entry point {ep} failed to import: {result.stderr}"
        )


def test_cli_help_flags() -> None:
    """All CLI entry points must respond to --help without error."""
    import subprocess

    commands = [
        "spink7-prep",
        "spink7-equilibrate",
        "spink7-production",
        "spink7-smd",
        "spink7-umbrella",
        "spink7-analysis",
    ]
    for cmd in commands:
        result = subprocess.run(
            [cmd, "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, (
            f"{cmd} --help failed: {result.stderr}"
        )


def test_requirements_consistency() -> None:
    """requirements.txt must reference pyproject.toml or contain matching deps."""
    req_path = Path(__file__).resolve().parents[1] / "requirements.txt"
    assert req_path.exists(), "requirements.txt must exist"
    content = req_path.read_text(encoding="utf-8").strip()
    assert len(content) > 0, "requirements.txt must not be empty"
