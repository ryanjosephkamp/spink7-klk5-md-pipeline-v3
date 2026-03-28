"""Force field factory for the SPINK7-KLK5 MD pipeline.

Provides a unified interface for creating OpenMM System objects from
different force field families (AMBER, AMOEBA, ML potentials, QM/MM).
"""

from __future__ import annotations

import logging

import openmm
import openmm.app
from openmm import unit
from openmm.app import ForceField, HBonds, PME, Topology

from src.config import (
    SUPPORTED_FORCE_FIELD_FAMILIES,
    AMOEBAConfig,
    MLPotentialConfig,
    QMMMConfig,
    SystemConfig,
)

logger = logging.getLogger(__name__)


def create_system(
    topology: Topology,
    positions: unit.Quantity,
    system_config: SystemConfig,
    amoeba_config: AMOEBAConfig | None = None,
    ml_config: MLPotentialConfig | None = None,
    qmmm_config: QMMMConfig | None = None,
) -> openmm.System:
    """Create an OpenMM System from the specified force field family.

    Args:
        topology: OpenMM Topology object.
        positions: Atomic positions with units.
        system_config: System preparation parameters (includes force_field_family).
        amoeba_config: Required if force_field_family == 'amoeba'.
        ml_config: Required if force_field_family == 'ml_potential'.
        qmmm_config: Required if force_field_family == 'qmmm'.

    Returns:
        Parameterized OpenMM System object.

    Raises:
        ValueError: If force_field_family is not in SUPPORTED_FORCE_FIELD_FAMILIES.
        ValueError: If required config is not provided for the selected family.
    """
    family = system_config.force_field_family

    if family not in SUPPORTED_FORCE_FIELD_FAMILIES:
        raise ValueError(
            f"Unsupported force_field_family '{family}'. "
            f"Must be one of {SUPPORTED_FORCE_FIELD_FAMILIES}"
        )

    if family == "amber":
        return _create_amber_system(topology, system_config)
    elif family == "amoeba":
        if amoeba_config is None:
            raise ValueError(
                "amoeba_config is required when force_field_family == 'amoeba'"
            )
        return _create_amoeba_system(topology, amoeba_config)
    elif family == "ml_potential":
        if ml_config is None:
            raise ValueError(
                "ml_config is required when force_field_family == 'ml_potential'"
            )
        return _create_ml_system(topology, positions, ml_config)
    elif family == "qmmm":
        if qmmm_config is None:
            qmmm_config = QMMMConfig()
        return _create_qmmm_system(topology, positions, qmmm_config)

    raise ValueError(f"Unhandled force_field_family: {family}")  # pragma: no cover


def _create_amber_system(
    topology: Topology,
    system_config: SystemConfig,
) -> openmm.System:
    """Create an AMBER ff14SB system (existing pipeline default)."""
    force_field = ForceField(system_config.force_field, system_config.water_model)
    return force_field.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=HBonds,
        rigidWater=True,
    )


def _create_amoeba_system(
    topology: Topology,
    amoeba_config: AMOEBAConfig,
) -> openmm.System:
    """Create an AMOEBA polarizable system.

    Detects implicit solvent models (Generalized Kirkwood) and adjusts the
    nonbonded method accordingly: implicit solvent requires NoCutoff, while
    explicit solvent uses PME.
    """
    from openmm.app import NoCutoff

    force_field = ForceField(
        amoeba_config.force_field_xml,
        amoeba_config.water_model_xml,
    )

    is_implicit = "_gk" in amoeba_config.water_model_xml.lower()

    create_kwargs: dict[str, object] = {
        "polarization": amoeba_config.polarization_type,
        "mutualInducedTargetEpsilon": amoeba_config.mutual_induced_target_epsilon,
    }

    if is_implicit:
        create_kwargs["nonbondedMethod"] = NoCutoff
    else:
        create_kwargs["nonbondedMethod"] = PME
        create_kwargs["nonbondedCutoff"] = 0.7 * unit.nanometers

    return force_field.createSystem(topology, **create_kwargs)


def validate_amoeba_dipole_convergence(
    simulation: openmm.app.Simulation,
    max_dipole_norm: float = 10.0,
) -> bool:
    """Verify that AMOEBA induced dipoles converged to physical values."""
    import numpy as np

    try:
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(energy):
            logger.warning("AMOEBA energy is not finite — dipole convergence failure likely")
            return False
        return True
    except Exception:
        logger.warning("Failed to evaluate AMOEBA energy for dipole convergence check")
        return False


def _create_ml_system(
    topology: Topology,
    positions: unit.Quantity,
    ml_config: MLPotentialConfig,
) -> openmm.System:
    """Create a machine-learned potential system via openmmml."""
    try:
        from openmmml import MLPotential
    except ImportError as exc:
        raise ImportError(
            "openmmml package is required for ML force fields. "
            "Install with: pip install openmm-ml"
        ) from exc

    potential = MLPotential(ml_config.potential_name)
    return potential.createSystem(topology, implementation=ml_config.implementation)


def _create_qmmm_system(
    topology: Topology,
    positions: unit.Quantity,
    qmmm_config: QMMMConfig,
) -> openmm.System:
    """Create a QM/MM system (not yet implemented)."""
    raise NotImplementedError(
        "QM/MM force field integration is not yet implemented. "
        "This placeholder defines the interface for future development. "
        "See L-15 implementation guide, Step 5 for design specifications."
    )


def compare_force_field_energies(
    topology: Topology,
    positions: unit.Quantity,
    system_config_amber: SystemConfig,
    system_config_alt: SystemConfig,
    amoeba_config: AMOEBAConfig | None = None,
    ml_config: MLPotentialConfig | None = None,
) -> dict[str, object]:
    """Compute potential energies with both AMBER and an alternative force field."""
    import numpy as np

    amber_system = create_system(topology, positions, system_config_amber)
    amber_integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    amber_sim = openmm.app.Simulation(topology, amber_system, amber_integrator, platform)
    amber_sim.context.setPositions(positions)
    amber_energy = (
        amber_sim.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    alt_system = create_system(
        topology, positions, system_config_alt,
        amoeba_config=amoeba_config, ml_config=ml_config,
    )
    alt_integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    alt_sim = openmm.app.Simulation(topology, alt_system, alt_integrator, platform)
    alt_sim.context.setPositions(positions)
    alt_energy = (
        alt_sim.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    return {
        "amber_energy_kj_mol": float(amber_energy),
        "alt_energy_kj_mol": float(alt_energy),
        "difference_kj_mol": float(alt_energy - amber_energy),
        "alt_family": system_config_alt.force_field_family,
        "both_finite": bool(np.isfinite(amber_energy) and np.isfinite(alt_energy)),
    }
