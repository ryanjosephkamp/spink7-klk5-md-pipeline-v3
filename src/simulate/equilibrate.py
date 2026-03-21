"""NVT and NPT equilibration utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import numpy as np
import openmm
from openmm import MonteCarloBarostat, XmlSerializer, unit
from openmm.app import DCDReporter, Simulation

from src import PhysicalValidityError
from src.analyze.equilibration import detect_equilibration
from src.config import BOLTZMANN_KJ, EquilibrationConfig
from src.physics.finite_size import compute_solute_net_charge, finite_size_correction
from src.simulate._topology_io import save_topology_pdb


logger = logging.getLogger(__name__)


def _resolve_seed(seed: int | None, stage_name: str) -> int:
    """Resolve a random seed, generating one from OS entropy if None.

    The resolved seed is always logged for post-hoc reproducibility.

    Args:
        seed: User-provided seed, or None for auto-generation.
        stage_name: Human-readable stage name for log messages.

    Returns:
        The resolved integer seed (32-bit unsigned).
    """
    if seed is None:
        seed = int.from_bytes(os.urandom(4), "big") & 0x7FFFFFFF
    logger.info("Using random seed %d for %s", seed, stage_name)
    return seed


def _validate_config(config: EquilibrationConfig) -> None:
    """Validate public boundary equilibration parameters."""

    if config.nvt_duration_ps <= 0.0:
        raise ValueError("config.nvt_duration_ps must be positive")
    if config.npt_duration_ps <= 0.0:
        raise ValueError("config.npt_duration_ps must be positive")
    if config.temperature_k <= 0.0:
        raise ValueError("config.temperature_k must be positive")
    if config.friction_per_ps <= 0.0:
        raise ValueError("config.friction_per_ps must be positive")
    if config.timestep_ps <= 0.0:
        raise ValueError("config.timestep_ps must be positive")
    if config.pressure_atm <= 0.0:
        raise ValueError("config.pressure_atm must be positive")
    if config.barostat_interval <= 0:
        raise ValueError("config.barostat_interval must be positive")
    if config.save_interval_ps <= 0.0:
        raise ValueError("config.save_interval_ps must be positive")


def _output_paths(output_dir: Path, stage_name: str) -> tuple[Path, Path]:
    """Create stage-specific output paths for trajectory and final state."""

    output_dir.mkdir(parents=True, exist_ok=True)
    trajectory_path = output_dir / f"{stage_name}_equilibration.dcd"
    final_state_path = output_dir / f"{stage_name}_final_state.xml"
    return trajectory_path, final_state_path


def _set_integrator_state(simulation: Simulation, temperature_k: float, friction_per_ps: float, random_seed: int) -> None:
    """Configure the simulation integrator and starting velocities deterministically."""

    integrator = simulation.integrator
    if hasattr(integrator, "setTemperature"):
        integrator.setTemperature(temperature_k * unit.kelvin)
    if hasattr(integrator, "setFriction"):
        integrator.setFriction(friction_per_ps / unit.picosecond)
    if hasattr(integrator, "setRandomNumberSeed"):
        integrator.setRandomNumberSeed(random_seed)
    simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin, random_seed)


def _n_degrees_of_freedom(simulation: Simulation) -> int:
    """Estimate the active translational degrees of freedom for temperature calculation."""

    system = simulation.system
    n_particles = system.getNumParticles()
    n_constraints = system.getNumConstraints()
    has_cmm = any(isinstance(force, openmm.CMMotionRemover) for force in system.getForces())
    dof = 3 * n_particles - n_constraints - (3 if has_cmm else 0)
    if dof <= 0:
        raise PhysicalValidityError("Equilibration requires a system with positive degrees of freedom")
    return dof


def _temperature_from_state(kinetic_energy_kj_mol: float, n_dof: int) -> float:
    """Convert kinetic energy to instantaneous temperature in kelvin."""

    return (2.0 * kinetic_energy_kj_mol) / (n_dof * BOLTZMANN_KJ)


def _system_mass_amu(simulation: Simulation) -> float:
    """Compute the total system mass in atomic mass units."""

    return sum(
        simulation.system.getParticleMass(index).value_in_unit(unit.dalton)
        for index in range(simulation.system.getNumParticles())
    )


def _density_from_state(total_mass_amu: float, state: openmm.State) -> float:
    """Compute density in g/cm^3 from an OpenMM state volume."""

    volume_nm3 = state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
    return (total_mass_amu * 1.66053906660e-3) / volume_nm3


def _run_stage(
    simulation: Simulation,
    total_steps: int,
    report_interval_steps: int,
    trajectory_path: Path,
    include_density: bool,
) -> dict[str, np.ndarray]:
    """Advance a simulation stage while streaming trajectory frames and collecting statistics."""

    # Pre-allocate thermodynamic arrays for the expected number of frames.
    n_frames = total_steps // report_interval_steps
    temperature = np.zeros(n_frames, dtype=float)
    potential_energy = np.zeros(n_frames, dtype=float)
    kinetic_energy = np.zeros(n_frames, dtype=float)
    density = np.zeros(n_frames, dtype=float) if include_density else np.zeros(0, dtype=float)

    # Cache system-level constants used for per-frame temperature and density.
    n_dof = _n_degrees_of_freedom(simulation)
    total_mass_amu = _system_mass_amu(simulation) if include_density else 0.0

    reporter = DCDReporter(str(trajectory_path), report_interval_steps)
    existing_reporters = list(simulation.reporters)
    simulation.reporters.append(reporter)

    try:
        # Advance one report-interval chunk per frame, collecting
        # instantaneous thermodynamic quantities after each chunk.
        for frame_index in range(n_frames):
            simulation.step(report_interval_steps)
            state = simulation.context.getState(getEnergy=True)
            potential_energy[frame_index] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            kinetic_energy[frame_index] = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            temperature[frame_index] = _temperature_from_state(kinetic_energy[frame_index], n_dof)
            if include_density:
                density[frame_index] = _density_from_state(total_mass_amu, state)
    finally:
        simulation.reporters[:] = existing_reporters

    return {
        "temperature": temperature,
        "potential_energy": potential_energy,
        "kinetic_energy": kinetic_energy,
        "density": density,
    }


def _write_final_state(simulation: Simulation, final_state_path: Path) -> None:
    """Serialize the final simulation state to XML."""

    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    final_state_path.write_text(XmlSerializer.serialize(state), encoding="utf-8")


def _equilibrated_segment(values: np.ndarray) -> tuple[np.ndarray, int]:
    """Detect equilibration and return the post-transient segment.

    Uses the Chodera maximum-effective-samples method to identify the
    optimal discard index.  Falls back to returning the full series if
    the timeseries is too short for meaningful detection.

    Returns:
        tuple: (equilibrated_values, t0_discard_index)
    """
    if values.size <= 4:
        return values, 0
    result = detect_equilibration(values)
    t0 = result["t0"]
    return values[t0:], t0


def run_nvt(
    simulation: Simulation,
    config: EquilibrationConfig,
    output_dir: Path,
) -> dict[str, Any]:
    """Run NVT equilibration and verify thermostat-driven temperature stability.

    Invariants: Enforces IV-2 by requiring the average temperature to remain within
    5 K of the 310 K target after equilibration.

    Args:
        simulation: Simulation with positional restraints applied.
        config: NVT/NPT equilibration parameters.
        output_dir: Directory for trajectory and state outputs.

    Returns:
        dict[str, Any]: Trajectory path, average and standard-deviation temperature,
        and the final serialized state path.
    """

    _validate_config(config)
    trajectory_path, final_state_path = _output_paths(Path(output_dir), "nvt")
    total_steps = max(1, round(config.nvt_duration_ps / config.timestep_ps))
    report_interval_steps = max(1, min(total_steps, round(config.save_interval_ps / config.timestep_ps)))
    total_steps = report_interval_steps * max(1, total_steps // report_interval_steps)

    logger.info("Starting NVT equilibration for %d steps", total_steps)
    nvt_seed = _resolve_seed(config.random_seed, "NVT equilibration")
    _set_integrator_state(simulation, config.temperature_k, config.friction_per_ps, nvt_seed)
    warmup_steps = max(report_interval_steps, total_steps // 4)
    simulation.step(warmup_steps)
    observables = _run_stage(simulation, total_steps, report_interval_steps, trajectory_path, include_density=False)

    equilibrated_temperature, t0_temperature = _equilibrated_segment(observables["temperature"])
    avg_temperature_k = float(np.mean(equilibrated_temperature))
    temperature_std_k = float(np.std(equilibrated_temperature))
    _write_final_state(simulation, final_state_path)

    if abs(avg_temperature_k - config.temperature_k) >= 5.0:
        logger.error("IV-2 violated: average temperature %.3f K", avg_temperature_k)
        raise PhysicalValidityError("IV-2 violated: NVT average temperature must remain within 5 K of target")

    return {
        "trajectory_path": trajectory_path,
        "avg_temperature_k": avg_temperature_k,
        "temperature_std_k": temperature_std_k,
        "final_state_path": final_state_path,
        "random_seed": nvt_seed,
        "t0_temperature": t0_temperature,
    }


def run_npt(
    simulation: Simulation,
    config: EquilibrationConfig,
    output_dir: Path,
) -> dict[str, Any]:
    """Run NPT equilibration and verify physical density and temperature.

    Invariants: Enforces IV-2 and IV-3 by requiring the average temperature to
    remain near the target and the average aqueous density to remain within
    [0.95, 1.05] g/cm^3.

    Args:
        simulation: Simulation with a periodic solvated system.
        config: NVT/NPT equilibration parameters.
        output_dir: Directory for trajectory and state outputs.

    Returns:
        dict[str, Any]: Trajectory path, average density and temperature, final
        box vectors in nm, and the final serialized state path.
    """

    _validate_config(config)
    trajectory_path, final_state_path = _output_paths(Path(output_dir), "npt")
    total_steps = max(1, round(config.npt_duration_ps / config.timestep_ps))
    report_interval_steps = max(1, min(total_steps, round(config.save_interval_ps / config.timestep_ps)))
    total_steps = report_interval_steps * max(1, total_steps // report_interval_steps)

    if not any(isinstance(force, MonteCarloBarostat) for force in simulation.system.getForces()):
        simulation.system.addForce(
            MonteCarloBarostat(
                config.pressure_atm * unit.atmosphere,
                config.temperature_k * unit.kelvin,
                config.barostat_interval,
            )
        )
        simulation.context.reinitialize(preserveState=True)

    logger.info("Starting NPT equilibration for %d steps", total_steps)
    npt_seed = _resolve_seed(config.random_seed, "NPT equilibration")
    _set_integrator_state(simulation, config.temperature_k, config.friction_per_ps, npt_seed)
    warmup_steps = max(report_interval_steps, total_steps // 4)
    simulation.step(warmup_steps)
    observables = _run_stage(simulation, total_steps, report_interval_steps, trajectory_path, include_density=True)

    equilibrated_temperature, t0_temperature = _equilibrated_segment(observables["temperature"])
    equilibrated_density, t0_density = _equilibrated_segment(observables["density"])
    avg_temperature_k = float(np.mean(equilibrated_temperature))
    avg_density_g_cm3 = float(np.mean(equilibrated_density))
    final_state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    box_vectors = final_state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    _write_final_state(simulation, final_state_path)
    topology_path = save_topology_pdb(simulation, Path(output_dir) / "topology_reference.pdb")

    if abs(avg_temperature_k - config.temperature_k) >= 5.0:
        logger.error("IV-2 violated during NPT: average temperature %.3f K", avg_temperature_k)
        raise PhysicalValidityError("IV-2 violated: NPT average temperature must remain within 5 K of target")
    if not 0.95 < avg_density_g_cm3 < 1.05:
        logger.error("IV-3 violated: average density %.6f g/cm^3", avg_density_g_cm3)
        raise PhysicalValidityError("IV-3 violated: NPT average density must remain within [0.95, 1.05] g/cm^3")

    # Compute finite-size electrostatic correction diagnostic.
    solute_indices = [
        i for i, atom in enumerate(simulation.topology.atoms())
        if atom.residue.name.upper() not in {"HOH", "WAT", "NA", "CL"}
    ]
    net_charge = compute_solute_net_charge(simulation, solute_indices)
    fs_result = finite_size_correction(net_charge, np.asarray(box_vectors, dtype=float))
    fs_correction = fs_result["correction_kj_mol"]
    logger.info("Finite-size correction: %.4f kJ/mol (Q=%.2f e, L_eff=%.3f nm)",
                fs_correction, net_charge, fs_result["effective_box_length_nm"])

    return {
        "trajectory_path": trajectory_path,
        "avg_density_g_cm3": avg_density_g_cm3,
        "avg_temperature_k": avg_temperature_k,
        "box_vectors_nm": np.asarray(box_vectors, dtype=float),
        "final_state_path": final_state_path,
        "topology_path": topology_path,
        "finite_size_correction_kj_mol": fs_correction,
        "random_seed": npt_seed,
        "t0_temperature": t0_temperature,
        "t0_density": t0_density,
    }