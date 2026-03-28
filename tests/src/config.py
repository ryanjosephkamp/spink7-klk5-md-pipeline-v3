"""Central configuration for the SPINK7-KLK5 MD pipeline.

All physical constants and simulation parameters are defined here.
Units follow the OpenMM internal convention (nm, ps, kJ/mol, K, e).
"""

import os
from dataclasses import dataclass, fields as dataclass_fields
from pathlib import Path

import yaml


# Physical constants
BOLTZMANN_KJ: float = 0.008314462618
AVOGADRO: float = 6.02214076e23
KCAL_TO_KJ: float = 4.184


# Root paths
PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
DATA_DIR: Path = PROJECT_ROOT / "data"
TRAJ_DIR: Path = DATA_DIR / "trajectories"


SUPPORTED_FORCE_FIELD_FAMILIES: tuple[str, ...] = ("amber", "amoeba", "ml_potential", "qmmm")

SUPPORTED_WATER_MODELS: dict[str, tuple[str, str]] = {
    "tip3p": ("amber14/tip3p.xml", "tip3p"),
    "opc": ("amber14/opc.xml", "tip4pew"),
    "tip4pew": ("amber14/tip4pew.xml", "tip4pew"),
}

SUPPORTED_BOX_SHAPES: tuple[str, ...] = ("cubic", "dodecahedron")


@dataclass(frozen=True)
class SystemConfig:
    """Immutable system preparation parameters."""

    force_field: str = "amber14-all.xml"
    water_model: str = "amber14/tip3p.xml"
    force_field_family: str = "amber"
    ph: float = 7.4
    box_padding_nm: float = 1.2
    ionic_strength_molar: float = 0.15
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"
    box_shape: str = "cubic"


@dataclass(frozen=True)
class MinimizationConfig:
    """Energy minimization parameters."""

    max_iterations: int = 10_000
    tolerance_kj_mol_nm: float = 10.0


@dataclass(frozen=True)
class EquilibrationConfig:
    """NVT/NPT equilibration parameters."""

    nvt_duration_ps: float = 500.0
    npt_duration_ps: float = 1000.0
    temperature_k: float = 310.0
    friction_per_ps: float = 1.0
    timestep_ps: float = 0.002
    pressure_atm: float = 1.0
    barostat_interval: int = 25
    restraint_k_kj_mol_nm2: float = 1000.0
    save_interval_ps: float = 10.0
    random_seed: int | None = None


@dataclass(frozen=True)
class ProductionConfig:
    """Production MD parameters."""

    duration_ns: float = 100.0
    temperature_k: float = 310.0
    friction_per_ps: float = 1.0
    timestep_ps: float = 0.002
    pressure_atm: float = 1.0
    save_interval_ps: float = 10.0
    checkpoint_interval_ps: float = 100.0
    random_seed: int | None = None


@dataclass(frozen=True)
class SMDConfig:
    """Steered Molecular Dynamics parameters."""

    spring_constant_kj_mol_nm2: float = 1000.0
    pulling_velocity_nm_per_ps: float = 0.001
    pull_distance_nm: float = 3.0
    n_replicates: int = 50
    save_interval_ps: float = 1.0
    random_seed: int | None = None


@dataclass(frozen=True)
class UmbrellaConfig:
    """Umbrella Sampling parameters."""

    xi_min_nm: float = 1.5
    xi_max_nm: float = 4.0
    window_spacing_nm: float = 0.05
    spring_constant_kj_mol_nm2: float = 1000.0
    per_window_duration_ns: float = 10.0
    save_interval_ps: float = 1.0
    pre_position_velocity_nm_per_ps: float = 0.01
    pre_position_spring_constant_kj_mol_nm2: float = 1000.0
    equilibration_duration_ps: float = 200.0
    detect_equilibration: bool = True


@dataclass(frozen=True)
class WHAMConfig:
    """WHAM solver parameters."""

    tolerance: float = 1e-6
    max_iterations: int = 100_000
    n_bootstrap: int = 200
    histogram_bins: int = 200


@dataclass(frozen=True)
class MetadynamicsConfig:
    """Well-tempered metadynamics parameters."""

    gaussian_height_kj_mol: float = 1.0
    gaussian_width_nm: float = 0.05
    deposition_interval_ps: float = 1.0
    bias_factor: float = 15.0
    temperature_k: float = 310.0
    simulation_duration_ns: float = 50.0
    save_interval_ps: float = 1.0
    convergence_tolerance_kj_mol: float = 1.0
    convergence_check_interval_ns: float = 1.0
    grid_min_nm: float = 1.0
    grid_max_nm: float = 5.0
    grid_num_bins: int = 200


@dataclass(frozen=True)
class MBARConfig:
    """MBAR solver parameters."""

    solver_protocol: str = "robust"
    relative_tolerance: float = 1e-7
    maximum_iterations: int = 10_000
    n_bootstrap: int = 200
    n_pmf_bins: int = 200


@dataclass(frozen=True)
class FEPConfig:
    """Alchemical free energy perturbation parameters."""

    n_lambda_windows: int = 20
    per_window_duration_ns: float = 2.0
    temperature_k: float = 310.0
    soft_core_alpha: float = 0.5
    soft_core_power: int = 1
    save_interval_ps: float = 1.0
    n_equilibration_steps: int = 5000
    annihilate_electrostatics: bool = True
    annihilate_sterics: bool = False


@dataclass(frozen=True)
class AMOEBAConfig:
    """AMOEBA polarizable force field parameters."""

    force_field_xml: str = "amoeba2018.xml"
    water_model_xml: str = "amoeba2018_gk.xml"
    mutual_induced_target_epsilon: float = 1e-5
    polarization_type: str = "mutual"
    pme_grid_spacing_nm: float = 0.06
    ewald_error_tolerance: float = 1e-5


@dataclass(frozen=True)
class MLPotentialConfig:
    """Machine-learned force field parameters."""

    potential_name: str = "ani2x"
    implementation: str = "torchani"
    mixed_precision: bool = True


@dataclass(frozen=True)
class QMMMConfig:
    """QM/MM hybrid potential parameters (stub)."""

    qm_method: str = "B3LYP"
    qm_basis_set: str = "6-31G*"
    qm_region_residues: tuple[str, ...] = ("HIS57", "ASP102", "SER195")


@dataclass(frozen=True)
class FiniteSizeCorrectionConfig:
    """Finite-size electrostatic correction parameters."""

    ewald_self_energy_constant: float = -2.837297
    solvent_dielectric: float = 80.0
    coulomb_constant_kj_nm_per_mol_e2: float = 138.935456
    apply_correction: bool = True


@dataclass(frozen=True)
class MSMConfig:
    """Markov State Model construction parameters."""

    tica_lag_ps: float = 10.0
    tica_n_components: int = 4
    n_clusters: int = 200
    lag_times_ps: tuple[float, ...] = (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0)
    n_implied_timescales: int = 5
    bayesian_n_samples: int = 100
    stride: int = 1


# ---------------------------------------------------------------------------
# Configuration loading (L-21)
# ---------------------------------------------------------------------------

_SECTION_TO_DATACLASS: dict[str, type] = {
    "system": SystemConfig,
    "minimization": MinimizationConfig,
    "equilibration": EquilibrationConfig,
    "production": ProductionConfig,
    "smd": SMDConfig,
    "umbrella": UmbrellaConfig,
    "wham": WHAMConfig,
    "mbar": MBARConfig,
    "amoeba": AMOEBAConfig,
    "ml_potential": MLPotentialConfig,
    "qmmm": QMMMConfig,
    "finite_size_correction": FiniteSizeCorrectionConfig,
    "fep": FEPConfig,
    "metadynamics": MetadynamicsConfig,
    "msm": MSMConfig,
}


def _build_dataclass_from_dict(dataclass_type: type, overrides: dict) -> object:
    """Construct a frozen dataclass, applying overrides on top of defaults."""
    valid_field_names = {f.name for f in dataclass_fields(dataclass_type)}
    unknown_keys = set(overrides.keys()) - valid_field_names
    if unknown_keys:
        raise TypeError(
            f"Unrecognized fields for {dataclass_type.__name__}: {sorted(unknown_keys)}"
        )

    kwargs = {}
    for field in dataclass_fields(dataclass_type):
        if field.name not in overrides:
            continue
        value = overrides[field.name]
        if value is None:
            kwargs[field.name] = None
            continue
        default_val = field.default
        if default_val is not dataclass_fields.__class__:
            if isinstance(default_val, float) and isinstance(value, (int, float)):
                value = float(value)
            elif isinstance(default_val, tuple) and isinstance(value, list):
                value = tuple(value)
        kwargs[field.name] = value

    return dataclass_type(**kwargs)


_ENV_PREFIX = "SPINK7"


def _apply_env_overrides(section_name: str, config_instance: object) -> object:
    """Override dataclass fields from environment variables."""
    from dataclasses import replace

    overrides = {}
    for field in dataclass_fields(config_instance):
        env_key = f"{_ENV_PREFIX}_{section_name}_{field.name}".upper()
        env_value = os.environ.get(env_key)
        if env_value is not None:
            default_val = field.default
            if isinstance(default_val, float):
                overrides[field.name] = float(env_value)
            elif isinstance(default_val, int):
                overrides[field.name] = int(env_value)
            elif isinstance(default_val, bool):
                overrides[field.name] = env_value.lower() in ("true", "1", "yes")
            elif default_val is None:
                try:
                    overrides[field.name] = int(env_value)
                except ValueError:
                    try:
                        overrides[field.name] = float(env_value)
                    except ValueError:
                        overrides[field.name] = env_value
            else:
                overrides[field.name] = env_value

    if overrides:
        return replace(config_instance, **overrides)
    return config_instance


def load_config(config_path: str | Path | None = None) -> dict[str, object]:
    """Load simulation configuration from a YAML file."""
    if config_path is None:
        configs = {name: cls() for name, cls in _SECTION_TO_DATACLASS.items()}
        for section_name in configs:
            configs[section_name] = _apply_env_overrides(section_name, configs[section_name])
        return configs

    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    if raw is None:
        raw = {}

    configs: dict[str, object] = {}
    for section_name, dataclass_type in _SECTION_TO_DATACLASS.items():
        section_data = raw.get(section_name, {})
        if section_data is None:
            section_data = {}
        configs[section_name] = _build_dataclass_from_dict(dataclass_type, section_data)

    for section_name in configs:
        configs[section_name] = _apply_env_overrides(section_name, configs[section_name])

    return configs
