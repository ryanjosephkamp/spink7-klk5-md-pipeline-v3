"""Central configuration for the SPINK7-KLK5 MD pipeline.

All physical constants and simulation parameters are defined here.
Units follow the OpenMM internal convention (nm, ps, kJ/mol, K, e).
"""

from dataclasses import dataclass
from pathlib import Path


# Physical constants
BOLTZMANN_KJ: float = 0.008314462618
AVOGADRO: float = 6.02214076e23
KCAL_TO_KJ: float = 4.184


# Root paths
PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
DATA_DIR: Path = PROJECT_ROOT / "data"
TRAJ_DIR: Path = DATA_DIR / "trajectories"


@dataclass(frozen=True)
class SystemConfig:
    """Immutable system preparation parameters."""

    force_field: str = "amber14-all.xml"
    water_model: str = "amber14/tip3p.xml"
    ph: float = 7.4
    box_padding_nm: float = 1.2
    ionic_strength_molar: float = 0.15
    positive_ion: str = "Na+"
    negative_ion: str = "Cl-"


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


@dataclass(frozen=True)
class SMDConfig:
    """Steered Molecular Dynamics parameters."""

    spring_constant_kj_mol_nm2: float = 1000.0
    pulling_velocity_nm_per_ps: float = 0.001
    pull_distance_nm: float = 3.0
    n_replicates: int = 50
    save_interval_ps: float = 1.0


@dataclass(frozen=True)
class UmbrellaConfig:
    """Umbrella Sampling parameters."""

    xi_min_nm: float = 1.5
    xi_max_nm: float = 4.0
    window_spacing_nm: float = 0.05
    spring_constant_kj_mol_nm2: float = 1000.0
    per_window_duration_ns: float = 10.0
    save_interval_ps: float = 1.0


@dataclass(frozen=True)
class WHAMConfig:
    """WHAM solver parameters."""

    tolerance: float = 1e-6
    max_iterations: int = 100_000
    n_bootstrap: int = 200
    histogram_bins: int = 200
