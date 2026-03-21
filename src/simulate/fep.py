"""Alchemical free energy perturbation for computational mutagenesis.

This module implements lambda-windowed alchemical transformations
for computing relative binding free energies (DDG) of point mutations
in the SPINK7-KLK5 system.

The alchemical approach uses OpenMM's openmmtools.alchemy module to
construct lambda-coupled Hamiltonians and pymbar for BAR/MBAR analysis.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import openmm
from openmm import unit

from src.config import FEPConfig
from src.simulate.platform import select_platform

logger = logging.getLogger(__name__)


def generate_lambda_schedule(n_windows: int) -> np.ndarray:
    """Generate a uniform lambda schedule for alchemical FEP.

    Args:
        n_windows: Number of lambda windows (must be >= 2).

    Returns:
        np.ndarray: Lambda values from 0.0 to 1.0. Shape: [n_windows].

    Raises:
        ValueError: If n_windows < 2.
    """
    if n_windows < 2:
        raise ValueError("n_windows must be at least 2")
    return np.linspace(0.0, 1.0, n_windows, dtype=float)


def create_alchemical_system(
    system: openmm.System,
    mutant_atom_indices: list[int],
    config: FEPConfig,
) -> openmm.System:
    """Construct an alchemical system with lambda-coupled interactions.

    Uses openmmtools.alchemy.AbsoluteAlchemicalFactory to create a system
    where the specified atoms can be alchemically transformed via a
    lambda coupling parameter.

    Args:
        system: The original OpenMM System.
        mutant_atom_indices: Indices of atoms in the mutated residue.
        config: FEP configuration parameters.

    Returns:
        openmm.System: The alchemical system with lambda-controlled forces.

    Raises:
        ValueError: If mutant_atom_indices is empty or contains
            out-of-range indices.
    """
    if not mutant_atom_indices:
        raise ValueError("mutant_atom_indices must not be empty")
    n_particles = system.getNumParticles()
    for idx in mutant_atom_indices:
        if idx < 0 or idx >= n_particles:
            raise ValueError(
                f"Atom index {idx} out of range [0, {n_particles})"
            )

    from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion

    alchemical_region = AlchemicalRegion(
        alchemical_atoms=mutant_atom_indices,
        annihilate_electrostatics=config.annihilate_electrostatics,
        annihilate_sterics=config.annihilate_sterics,
        softcore_alpha=config.soft_core_alpha,
        softcore_a=config.soft_core_power,
        softcore_b=config.soft_core_power,
        softcore_c=config.soft_core_power,
    )
    factory = AbsoluteAlchemicalFactory()
    return factory.create_alchemical_system(system, alchemical_region)


def run_fep_window(
    alchemical_system: openmm.System,
    positions: Any,
    lambda_value: float,
    lambda_schedule: np.ndarray,
    config: FEPConfig,
    output_dir: Path,
) -> dict[str, Any]:
    """Run a single lambda window and collect energy samples.

    Equilibrates the system at the given lambda value, then collects
    energy samples evaluated at all lambda states for MBAR analysis.

    Args:
        alchemical_system: The alchemical OpenMM System.
        positions: Initial atomic positions.
        lambda_value: The lambda value for this window.
        lambda_schedule: Full array of lambda values for energy re-evaluation.
        config: FEP configuration parameters.
        output_dir: Directory to save window output files.

    Returns:
        dict with keys:
            - "lambda": the lambda value for this window
            - "energies": np.ndarray of shape [n_lambda, n_samples],
              reduced potential energies at all lambda states
    """
    from openmmtools.alchemy import AlchemicalState

    output_dir.mkdir(parents=True, exist_ok=True)

    platform = select_platform()
    integrator = openmm.LangevinMiddleIntegrator(
        config.temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    simulation = openmm.app.Simulation(
        openmm.app.Topology(),
        alchemical_system,
        integrator,
        platform,
    )
    simulation.context.setPositions(positions)

    alchemical_state = AlchemicalState.from_system(alchemical_system)
    alchemical_state.lambda_electrostatics = lambda_value
    alchemical_state.lambda_sterics = lambda_value
    alchemical_state.apply_to_context(simulation.context)

    # Equilibration
    logger.info("Equilibrating lambda=%.3f (%d steps)", lambda_value, config.n_equilibration_steps)
    simulation.step(config.n_equilibration_steps)

    # Production sampling
    n_steps_total = int(config.per_window_duration_ns * 1e6 / 2.0)  # 2 fs timestep
    save_interval_steps = int(config.save_interval_ps * 1000 / 2.0)
    n_samples = n_steps_total // save_interval_steps

    beta = 1.0 / (0.008314462618 * config.temperature_k)  # 1/(kB*T) in mol/kJ
    energies = np.zeros((len(lambda_schedule), n_samples), dtype=float)

    for sample_idx in range(n_samples):
        simulation.step(save_interval_steps)
        state = simulation.context.getState(getPositions=True, getEnergy=True)

        # Evaluate energy at all lambda values
        for lam_idx, lam in enumerate(lambda_schedule):
            alchemical_state.lambda_electrostatics = lam
            alchemical_state.lambda_sterics = lam
            alchemical_state.apply_to_context(simulation.context)
            energy_state = simulation.context.getState(getEnergy=True)
            energies[lam_idx, sample_idx] = (
                energy_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole) * beta
            )

        # Reset to this window's lambda
        alchemical_state.lambda_electrostatics = lambda_value
        alchemical_state.lambda_sterics = lambda_value
        alchemical_state.apply_to_context(simulation.context)

    np.savez(
        output_dir / f"window_lambda_{lambda_value:.4f}.npz",
        lambda_value=lambda_value,
        energies=energies,
    )

    logger.info("Window lambda=%.3f complete: %d samples collected", lambda_value, n_samples)
    return {"lambda": lambda_value, "energies": energies}


def run_fep_campaign(
    system: openmm.System,
    positions: Any,
    mutant_atom_indices: list[int],
    config: FEPConfig,
    output_dir: Path,
) -> dict[str, Any]:
    """Orchestrate a full alchemical FEP campaign over all lambda windows.

    Constructs the alchemical system, generates the lambda schedule,
    runs each window, and assembles the full energy matrix for MBAR.

    Args:
        system: The original (wild-type) OpenMM System.
        positions: Initial atomic positions.
        mutant_atom_indices: Atom indices of the residue to transform.
        config: FEP configuration.
        output_dir: Directory for campaign outputs.

    Returns:
        dict with keys:
            - "lambda_schedule": np.ndarray of lambda values
            - "energy_matrix": np.ndarray [n_lambda, n_total_samples]
            - "n_samples_per_state": np.ndarray [n_lambda]
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    lambda_schedule = generate_lambda_schedule(config.n_lambda_windows)
    alchemical_system = create_alchemical_system(system, mutant_atom_indices, config)

    all_energies = []
    n_samples_list = []

    for lam in lambda_schedule:
        window_dir = output_dir / f"window_{lam:.4f}"
        result = run_fep_window(
            alchemical_system, positions, lam, lambda_schedule, config, window_dir,
        )
        all_energies.append(result["energies"])
        n_samples_list.append(result["energies"].shape[1])

    # Assemble full energy matrix: [n_lambda, total_samples]
    energy_matrix = np.concatenate(all_energies, axis=1)
    n_samples_per_state = np.array(n_samples_list, dtype=int)

    np.savez(
        output_dir / "fep_campaign_results.npz",
        lambda_schedule=lambda_schedule,
        energy_matrix=energy_matrix,
        n_samples_per_state=n_samples_per_state,
    )

    logger.info(
        "FEP campaign complete: %d windows, %d total samples",
        len(lambda_schedule),
        energy_matrix.shape[1],
    )

    return {
        "lambda_schedule": lambda_schedule,
        "energy_matrix": energy_matrix,
        "n_samples_per_state": n_samples_per_state,
    }
