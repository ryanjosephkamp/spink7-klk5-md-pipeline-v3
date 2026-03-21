"""Energy-minimization utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
import warnings
from typing import Dict

import numpy as np

from openmm import unit
from openmm.app import Simulation

from src import PhysicalValidityError
from src.config import MinimizationConfig


logger = logging.getLogger(__name__)


def minimize_energy(
    simulation: Simulation,
    config: MinimizationConfig,
) -> Dict[str, float]:
    """Perform steepest-descent energy minimization on an initialized simulation.

    Invariants: Enforces IV-1 by requiring the post-minimization potential energy
    to be strictly lower than the initial potential energy.

    Args:
        simulation: Initialized OpenMM simulation object.
        config: Energy-minimization parameters.

    Returns:
        Dict[str, float]: Initial energy, final energy, and requested iteration count.
    """

    if config.max_iterations <= 0:
        raise ValueError("config.max_iterations must be positive")
    if config.tolerance_kj_mol_nm <= 0.0:
        raise ValueError("config.tolerance_kj_mol_nm must be positive")

    logger.info(
        "Starting energy minimization with max_iterations=%d and tolerance=%.3f kJ/mol/nm",
        config.max_iterations,
        config.tolerance_kj_mol_nm,
    )

    initial_state = simulation.context.getState(getEnergy=True)
    initial_energy_kj_mol = initial_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    simulation.minimizeEnergy(
        tolerance=config.tolerance_kj_mol_nm * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=config.max_iterations,
    )

    final_state = simulation.context.getState(getEnergy=True, getForces=True)
    final_energy_kj_mol = final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    # Compute maximum force magnitude across all particles
    forces = final_state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )
    max_force_kj_mol_nm = float(np.max(np.linalg.norm(forces, axis=1)))

    # Convergence: the minimizer reached the force tolerance
    converged = max_force_kj_mol_nm < config.tolerance_kj_mol_nm

    # Energy reduction ratio
    energy_reduction = (
        (initial_energy_kj_mol - final_energy_kj_mol) / abs(initial_energy_kj_mol)
        if abs(initial_energy_kj_mol) > 1e-10
        else 0.0
    )

    if final_energy_kj_mol >= initial_energy_kj_mol:
        logger.error(
            "IV-1 violated during minimization: final_energy=%.6f kJ/mol, initial_energy=%.6f kJ/mol",
            final_energy_kj_mol,
            initial_energy_kj_mol,
        )
        raise PhysicalValidityError(
            "IV-1 violated: post-minimization potential energy must be lower than initial energy"
        )

    if not converged:
        warnings.warn(
            f"Energy minimization did not converge: max force "
            f"{max_force_kj_mol_nm:.4f} kJ/mol/nm exceeds tolerance "
            f"{config.tolerance_kj_mol_nm:.4f} kJ/mol/nm after "
            f"{config.max_iterations} iterations.",
            stacklevel=2,
        )

    logger.info(
        "Energy minimization complete: initial=%.6f kJ/mol, final=%.6f kJ/mol",
        initial_energy_kj_mol,
        final_energy_kj_mol,
    )
    return {
        "initial_energy_kj_mol": float(initial_energy_kj_mol),
        "final_energy_kj_mol": float(final_energy_kj_mol),
        "n_steps": int(config.max_iterations),
        "converged": converged,
        "max_force_kj_mol_nm": max_force_kj_mol_nm,
        "energy_reduction": energy_reduction,
    }