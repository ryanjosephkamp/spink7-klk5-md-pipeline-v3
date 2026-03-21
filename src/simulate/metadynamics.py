"""Well-tempered metadynamics utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Any

import numpy as np
import openmm
from openmm import unit
from openmm.app import Simulation
from openmm.app.metadynamics import BiasVariable, Metadynamics

from src.config import BOLTZMANN_KJ, MetadynamicsConfig

logger = logging.getLogger(__name__)

__all__ = [
    "run_metadynamics",
    "_validate_config",
    "_create_com_distance_bias_variable",
    "_extract_fes",
    "_check_fes_convergence",
]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _validate_config(config: MetadynamicsConfig) -> None:
    """Validate public-boundary metadynamics parameters."""
    if config.gaussian_height_kj_mol <= 0.0:
        raise ValueError("config.gaussian_height_kj_mol must be positive")
    if config.gaussian_width_nm <= 0.0:
        raise ValueError("config.gaussian_width_nm must be positive")
    if config.deposition_interval_ps <= 0.0:
        raise ValueError("config.deposition_interval_ps must be positive")
    if config.bias_factor <= 1.0:
        raise ValueError("config.bias_factor must be greater than 1.0")
    if config.temperature_k <= 0.0:
        raise ValueError("config.temperature_k must be positive")
    if config.simulation_duration_ns <= 0.0:
        raise ValueError("config.simulation_duration_ns must be positive")
    if config.save_interval_ps <= 0.0:
        raise ValueError("config.save_interval_ps must be positive")
    if config.grid_min_nm >= config.grid_max_nm:
        raise ValueError("config.grid_min_nm must be less than config.grid_max_nm")
    if config.grid_num_bins < 2:
        raise ValueError("config.grid_num_bins must be at least 2")


# ---------------------------------------------------------------------------
# Bias variable construction
# ---------------------------------------------------------------------------


def _create_com_distance_bias_variable(
    group_a_indices: list[int],
    group_b_indices: list[int],
    config: MetadynamicsConfig,
) -> tuple[openmm.CustomCentroidBondForce, BiasVariable]:
    """Create an OpenMM BiasVariable for the COM distance between two groups.

    Returns the underlying CustomCentroidBondForce and the BiasVariable wrapper.
    The force must be added to the system before constructing the Metadynamics object.
    """
    force = openmm.CustomCentroidBondForce(2, "distance(g1, g2)")
    group_a = force.addGroup(group_a_indices)
    group_b = force.addGroup(group_b_indices)
    force.addBond([group_a, group_b], [])

    bias_var = BiasVariable(
        force,
        minValue=config.grid_min_nm,
        maxValue=config.grid_max_nm,
        biasWidth=config.gaussian_width_nm,
        periodic=False,
        gridWidth=config.grid_num_bins,
    )
    return force, bias_var


# ---------------------------------------------------------------------------
# Free energy surface extraction
# ---------------------------------------------------------------------------


def _extract_fes(
    meta: Metadynamics,
    config: MetadynamicsConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract the corrected free energy surface from accumulated bias.

    Applies the well-tempered correction:
        G(xi) = -(gamma / (gamma - 1)) * V_bias(xi)
    and shifts the FES so that G(xi_max) = 0.
    """
    raw_fes = meta.getFreeEnergy()
    # getFreeEnergy returns a 1D array of -V_bias values on the grid in kJ/mol
    fes_values = np.array(raw_fes, dtype=float)
    # The OpenMM Metadynamics.getFreeEnergy() already applies the
    # well-tempered scaling internally, returning the free energy estimate.
    # We just need to extract and shift.

    grid_nm = np.linspace(config.grid_min_nm, config.grid_max_nm, len(fes_values))
    # Shift so that the value at maximum xi is zero
    fes_values = fes_values - fes_values[-1]
    return grid_nm, fes_values


# ---------------------------------------------------------------------------
# Convergence diagnostics
# ---------------------------------------------------------------------------


def _check_fes_convergence(
    current_fes: np.ndarray,
    previous_fes: np.ndarray,
    tolerance_kj_mol: float,
) -> tuple[bool, float]:
    """Compare two FES snapshots and check for convergence.

    Returns (is_converged, max_absolute_deviation_kj_mol).
    """
    deviation = np.abs(current_fes - previous_fes)
    max_delta = float(np.max(deviation))
    return max_delta <= tolerance_kj_mol, max_delta


# ---------------------------------------------------------------------------
# Main simulation function
# ---------------------------------------------------------------------------


def run_metadynamics(
    simulation: Simulation,
    config: MetadynamicsConfig,
    output_dir: Path,
    pull_group_1: list[int],
    pull_group_2: list[int],
) -> dict[str, Any]:
    """Run well-tempered metadynamics on a COM-distance collective variable.

    Parameters
    ----------
    simulation : Simulation
        Initialized OpenMM Simulation with positions set.
    config : MetadynamicsConfig
        Metadynamics parameters.
    output_dir : Path
        Directory for output files.
    pull_group_1, pull_group_2 : list[int]
        Atom indices defining the two COM groups.

    Returns
    -------
    dict with keys:
        fes_grid_nm, fes_kj_mol, xi_timeseries,
        convergence_history, converged, output_dir
    """
    _validate_config(config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build bias variable
    force, bias_var = _create_com_distance_bias_variable(
        pull_group_1, pull_group_2, config,
    )

    # Compute deposition frequency in simulation steps
    timestep_ps = simulation.integrator.getStepSize().value_in_unit(unit.picoseconds)
    deposition_freq = max(1, int(round(config.deposition_interval_ps / timestep_ps)))
    save_freq = max(1, int(round(config.save_interval_ps / timestep_ps)))
    convergence_check_freq = max(
        1, int(round(config.convergence_check_interval_ns * 1000.0 / timestep_ps))
    )
    total_steps = int(round(config.simulation_duration_ns * 1000.0 / timestep_ps))

    # Construct the Metadynamics object
    meta = Metadynamics(
        simulation.system,
        [bias_var],
        config.temperature_k * unit.kelvin,
        config.bias_factor,
        config.gaussian_height_kj_mol * unit.kilojoules_per_mole,
        deposition_freq,
    )
    # Reinitialize context after adding the bias force
    simulation.context.reinitialize(preserveState=True)

    # Simulation loop
    xi_timeseries: list[float] = []
    convergence_history: list[tuple[float, float]] = []
    converged = False
    previous_fes: np.ndarray | None = None

    steps_done = 0
    while steps_done < total_steps:
        # Advance by save_freq steps
        chunk = min(save_freq, total_steps - steps_done)
        meta.step(simulation, chunk)
        steps_done += chunk

        # Record current xi
        cv_value = meta.getCollectiveVariables(simulation)
        xi_timeseries.append(float(cv_value[0]))

        # Convergence check
        if steps_done % convergence_check_freq == 0 or steps_done >= total_steps:
            time_ns = steps_done * timestep_ps / 1000.0
            grid_nm, fes_kj_mol = _extract_fes(meta, config)
            if previous_fes is not None:
                is_conv, max_delta = _check_fes_convergence(
                    fes_kj_mol, previous_fes, config.convergence_tolerance_kj_mol,
                )
                convergence_history.append((time_ns, max_delta))
                logger.info(
                    "Convergence check at %.1f ns: max delta = %.3f kJ/mol",
                    time_ns, max_delta,
                )
                if is_conv and not converged:
                    converged = True
                    logger.info(
                        "Metadynamics converged at %.1f ns (tolerance: %.3f kJ/mol)",
                        time_ns, config.convergence_tolerance_kj_mol,
                    )
            previous_fes = fes_kj_mol.copy()

    # Final FES extraction
    grid_nm, fes_kj_mol = _extract_fes(meta, config)
    xi_array = np.asarray(xi_timeseries, dtype=float)

    # Save outputs
    np.save(output_dir / "metadynamics_fes.npy", np.column_stack([grid_nm, fes_kj_mol]))
    np.save(output_dir / "metadynamics_xi_timeseries.npy", xi_array)

    if convergence_history:
        with (output_dir / "metadynamics_convergence.csv").open(
            "w", encoding="utf-8", newline="",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(["time_ns", "max_delta_kj_mol"])
            writer.writerows(convergence_history)

    return {
        "fes_grid_nm": grid_nm,
        "fes_kj_mol": fes_kj_mol,
        "xi_timeseries": xi_array,
        "convergence_history": convergence_history,
        "converged": converged,
        "output_dir": output_dir,
    }
