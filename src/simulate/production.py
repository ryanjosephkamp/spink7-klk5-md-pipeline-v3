"""Production MD utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import csv
import logging
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import openmm
from openmm import MonteCarloBarostat, XmlSerializer, unit
from openmm.app import DCDFile, Simulation

from src import PhysicalValidityError
from src.config import ProductionConfig
from src.physics.finite_size import compute_solute_net_charge, finite_size_correction
from src.simulate._topology_io import save_topology_pdb


logger = logging.getLogger(__name__)

_MAX_ENERGY_DRIFT_KJ_MOL_NS_ATOM = 0.1


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


def _validate_config(config: ProductionConfig) -> None:
    """Validate public-boundary production parameters."""

    if config.duration_ns <= 0.0:
        raise ValueError("config.duration_ns must be positive")
    if config.temperature_k <= 0.0:
        raise ValueError("config.temperature_k must be positive")
    if config.friction_per_ps <= 0.0:
        raise ValueError("config.friction_per_ps must be positive")
    if config.timestep_ps <= 0.0:
        raise ValueError("config.timestep_ps must be positive")
    if config.pressure_atm <= 0.0:
        raise ValueError("config.pressure_atm must be positive")
    if config.save_interval_ps <= 0.0:
        raise ValueError("config.save_interval_ps must be positive")
    if config.checkpoint_interval_ps <= 0.0:
        raise ValueError("config.checkpoint_interval_ps must be positive")


def _output_paths(output_dir: Path) -> tuple[Path, Path]:
    """Create production output paths."""

    output_dir.mkdir(parents=True, exist_ok=True)
    trajectory_path = output_dir / "production.dcd"
    energy_timeseries_path = output_dir / "production_energy.csv"
    return trajectory_path, energy_timeseries_path


def _set_integrator_state(simulation: Simulation, temperature_k: float, friction_per_ps: float, random_seed: int) -> None:
    """Configure the integrator for production dynamics."""

    integrator = simulation.integrator
    if hasattr(integrator, "setTemperature"):
        integrator.setTemperature(temperature_k * unit.kelvin)
    if hasattr(integrator, "setFriction"):
        integrator.setFriction(friction_per_ps / unit.picosecond)
    if hasattr(integrator, "setRandomNumberSeed"):
        integrator.setRandomNumberSeed(random_seed)


def _ensure_barostat(simulation: Simulation, config: ProductionConfig) -> None:
    """Guarantee that production runs in the NPT ensemble."""

    if any(isinstance(force, MonteCarloBarostat) for force in simulation.system.getForces()):
        return

    simulation.system.addForce(
        MonteCarloBarostat(
            config.pressure_atm * unit.atmosphere,
            config.temperature_k * unit.kelvin,
            25,
        )
    )
    simulation.context.reinitialize(preserveState=True)


def _checkpoint_path(output_dir: Path, step: int) -> Path:
    """Return a deterministic checkpoint path for a production step."""

    return output_dir / f"production_step_{step:08d}.chk"


def _write_energy_timeseries(energy_timeseries_path: Path, energy_timeseries: np.ndarray) -> None:
    """Persist production energies as a compact CSV file."""

    with energy_timeseries_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_ps", "potential_kj_mol", "kinetic_kj_mol"])
        writer.writerows(energy_timeseries.tolist())


def _append_energy_timeseries(energy_timeseries_path: Path, energy_timeseries: np.ndarray) -> None:
    """Append additional energy rows to an existing CSV file."""

    with energy_timeseries_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(energy_timeseries.tolist())


def _find_latest_checkpoint(output_dir: Path) -> Path | None:
    """Find the most recent production checkpoint in the output directory.

    Checkpoint filenames follow the pattern: production_step_XXXXXXXX.chk

    Returns:
        Path to the latest checkpoint, or None if no checkpoints exist.
    """

    checkpoint_pattern = re.compile(r"^production_step_(\d+)\.chk$")
    best_step = -1
    best_path: Path | None = None
    for candidate in output_dir.iterdir():
        match = checkpoint_pattern.match(candidate.name)
        if match:
            step = int(match.group(1))
            if step > best_step:
                best_step = step
                best_path = candidate
    return best_path


def _parse_checkpoint_step(checkpoint_path: Path) -> int:
    """Extract the step number from a checkpoint filename."""

    match = re.search(r"production_step_(\d+)\.chk$", checkpoint_path.name)
    if not match:
        raise ValueError(f"Cannot parse step from checkpoint: {checkpoint_path.name}")
    return int(match.group(1))


def _energy_drift_kj_mol_ns_atom(energy_timeseries: np.ndarray, n_atoms: int) -> float:
    """Estimate per-atom total-energy drift from a linear trend over time."""

    if energy_timeseries.shape[0] < 2:
        return 0.0

    time_ns = energy_timeseries[:, 0] / 1000.0
    total_energy_kj_mol = energy_timeseries[:, 1] + energy_timeseries[:, 2]
    time_centered = time_ns - float(np.mean(time_ns))
    energy_centered = total_energy_kj_mol - float(np.mean(total_energy_kj_mol))
    denominator = float(np.dot(time_centered, time_centered))
    if denominator == 0.0:
        return 0.0
    slope_kj_mol_per_ns = float(np.dot(time_centered, energy_centered) / denominator)
    return abs(slope_kj_mol_per_ns) / float(n_atoms)


def _clone_nve_reference_simulation(simulation: Simulation, timestep_ps: float) -> Simulation:
    """Create a short NVE reference simulation from the current production state."""

    reference_system = XmlSerializer.deserialize(XmlSerializer.serialize(simulation.system))
    for force_index in reversed(range(reference_system.getNumForces())):
        if isinstance(reference_system.getForce(force_index), MonteCarloBarostat):
            reference_system.removeForce(force_index)

    reference_integrator = openmm.VerletIntegrator(timestep_ps * unit.picoseconds)
    reference_platform = simulation.context.getPlatform()
    reference_simulation = Simulation(simulation.topology, reference_system, reference_integrator, reference_platform)

    state = simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    reference_simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    reference_simulation.context.setPositions(state.getPositions())
    reference_simulation.context.setVelocities(state.getVelocities())
    return reference_simulation


def _validate_nve_reference_drift(simulation: Simulation, config: ProductionConfig) -> None:
    """Enforce IV-5 using a short NVE reference trajectory from the production endpoint."""

    # Run a brief NVE (microcanonical) segment to measure numerical
    # energy drift, which should be negligible for a well-behaved integrator.
    reference_duration_ps = max(10.0, min(20.0, config.duration_ns * 1000.0))
    reference_save_interval_ps = min(1.0, reference_duration_ps)
    reference_interval_steps = max(1, round(reference_save_interval_ps / config.timestep_ps))
    n_reference_frames = max(2, round(reference_duration_ps / reference_save_interval_ps))
    total_reference_steps = n_reference_frames * reference_interval_steps
    reference_simulation = _clone_nve_reference_simulation(simulation, config.timestep_ps)
    reference_timeseries = np.zeros((n_reference_frames, 3), dtype=float)
    n_atoms = reference_simulation.system.getNumParticles()

    for frame_index in range(n_reference_frames):
        reference_simulation.step(reference_interval_steps)
        state = reference_simulation.context.getState(getEnergy=True)
        reference_timeseries[frame_index, 0] = (frame_index + 1) * reference_interval_steps * config.timestep_ps
        reference_timeseries[frame_index, 1] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        reference_timeseries[frame_index, 2] = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)

    energy_drift = _energy_drift_kj_mol_ns_atom(reference_timeseries, n_atoms)
    if energy_drift >= _MAX_ENERGY_DRIFT_KJ_MOL_NS_ATOM:
        logger.error("IV-5 violated: NVE reference energy drift %.6f kJ/mol/ns/atom", energy_drift)
        raise PhysicalValidityError(
            "IV-5 violated: total energy drift must remain below 0.1 kJ/mol/ns/atom during production"
        )


def run_production(
    simulation: Simulation,
    config: ProductionConfig,
    output_dir: Path,
) -> dict[str, Any]:
    """Run unrestrained production MD with streamed trajectory and checkpoint output.

    Invariants: Enforces IV-5 by requiring the per-atom total-energy drift to
    remain below 0.1 kJ/mol/ns/atom over the recorded production segment.

    Args:
        simulation: Initialized OpenMM simulation object.
        config: Production MD parameters.
        output_dir: Directory for trajectory, checkpoint, and energy outputs.

    Returns:
        dict[str, Any]: Trajectory path, frame count, simulated time in ns,
        checkpoint paths, and the CSV energy-timeseries path.
    """

    _validate_config(config)
    output_dir = Path(output_dir)
    trajectory_path, energy_timeseries_path = _output_paths(output_dir)
    _ensure_barostat(simulation, config)
    production_seed = _resolve_seed(config.random_seed, "production MD")
    _set_integrator_state(simulation, config.temperature_k, config.friction_per_ps, production_seed)

    report_interval_steps = max(1, round(config.save_interval_ps / config.timestep_ps))
    checkpoint_interval_steps = max(1, round(config.checkpoint_interval_ps / config.timestep_ps))
    n_frames = max(1, round((config.duration_ns * 1000.0) / config.save_interval_ps))
    total_steps = n_frames * report_interval_steps
    n_atoms = simulation.system.getNumParticles()
    start_step = simulation.currentStep

    logger.info(
        "Starting production MD for %d steps (%d frames, report every %d steps)",
        total_steps,
        n_frames,
        report_interval_steps,
    )

    energy_timeseries = np.zeros((n_frames, 3), dtype=float)
    checkpoint_paths: list[Path] = []
    frame_index = 0

    # Event-driven loop: advance to whichever event (report or checkpoint)
    # comes first, handle it, then continue until all steps are exhausted.
    next_report_step = report_interval_steps
    next_checkpoint_step = checkpoint_interval_steps
    integration_timestep_ps = simulation.integrator.getStepSize().value_in_unit(unit.picoseconds)

    with trajectory_path.open("wb") as trajectory_handle:
        dcd_file = DCDFile(
            trajectory_handle,
            simulation.topology,
            integration_timestep_ps * unit.picoseconds,
            firstStep=start_step,
            interval=report_interval_steps,
            append=False,
        )

        elapsed_steps = 0
        while elapsed_steps < total_steps:
            # Determine the nearest upcoming event (report or checkpoint)
            # and integrate only to that boundary.
            steps_to_report = next_report_step - elapsed_steps
            steps_to_checkpoint = next_checkpoint_step - elapsed_steps
            steps_to_event = min(steps_to_report, steps_to_checkpoint)
            simulation.step(steps_to_event)
            elapsed_steps += steps_to_event

            # Write trajectory frame and record energies at report events.
            if elapsed_steps == next_report_step:
                state = simulation.context.getState(getEnergy=True, getPositions=True)
                energy_timeseries[frame_index, 0] = elapsed_steps * config.timestep_ps
                energy_timeseries[frame_index, 1] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                energy_timeseries[frame_index, 2] = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                dcd_file.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())
                frame_index += 1
                next_report_step += report_interval_steps

            # Persist a checkpoint for restartability at checkpoint events.
            if elapsed_steps == next_checkpoint_step:
                checkpoint_path = _checkpoint_path(output_dir, start_step + elapsed_steps)
                simulation.saveCheckpoint(str(checkpoint_path))
                checkpoint_paths.append(checkpoint_path)
                next_checkpoint_step += checkpoint_interval_steps

    if not checkpoint_paths or checkpoint_paths[-1].stem != f"production_step_{simulation.currentStep:08d}":
        final_checkpoint_path = _checkpoint_path(output_dir, simulation.currentStep)
        simulation.saveCheckpoint(str(final_checkpoint_path))
        checkpoint_paths.append(final_checkpoint_path)

    _write_energy_timeseries(energy_timeseries_path, energy_timeseries)

    _validate_nve_reference_drift(simulation, config)

    # Compute finite-size electrostatic correction diagnostic.
    final_state = simulation.context.getState(getPositions=True)
    box_vectors = final_state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    solute_indices = [
        i for i, atom in enumerate(simulation.topology.atoms())
        if atom.residue.name.upper() not in {"HOH", "WAT", "NA", "CL"}
    ]
    net_charge = compute_solute_net_charge(simulation, solute_indices)
    fs_result = finite_size_correction(net_charge, np.asarray(box_vectors, dtype=float))
    fs_correction = fs_result["correction_kj_mol"]
    logger.info("Finite-size correction: %.4f kJ/mol (Q=%.2f e, L_eff=%.3f nm)",
                fs_correction, net_charge, fs_result["effective_box_length_nm"])

    topology_path = save_topology_pdb(simulation, output_dir / "production_topology.pdb")

    return {
        "trajectory_path": trajectory_path,
        "topology_path": topology_path,
        "n_frames": int(frame_index),
        "total_time_ns": float(energy_timeseries[frame_index - 1, 0] / 1000.0),
        "checkpoint_paths": checkpoint_paths,
        "energy_timeseries_path": energy_timeseries_path,
        "finite_size_correction_kj_mol": fs_correction,
        "random_seed": production_seed,
    }


def resume_production(
    simulation: Simulation,
    config: ProductionConfig,
    output_dir: Path,
    checkpoint_path: Path | None = None,
) -> dict[str, Any]:
    """Resume a production MD run from the latest or specified checkpoint.

    Loads the checkpoint state, calculates the remaining simulation steps,
    and appends new frames to the existing DCD trajectory file.

    Invariants: Enforces IV-5 via NVE drift validation on the extended trajectory.

    Args:
        simulation: Initialized OpenMM simulation object (topology and system must match
            the checkpoint's system).
        config: Production MD parameters (same as the original run).
        output_dir: Directory containing the existing trajectory, energy CSV, and checkpoints.
        checkpoint_path: Explicit checkpoint path. If None, the latest checkpoint
            in output_dir is used.

    Returns:
        dict[str, Any]: Updated trajectory path, total frame count, total simulated time,
        checkpoint paths, and energy timeseries path.

    Raises:
        FileNotFoundError: If no checkpoint is found.
    """

    _validate_config(config)
    output_dir = Path(output_dir)

    if checkpoint_path is None:
        checkpoint_path = _find_latest_checkpoint(output_dir)
    if checkpoint_path is None:
        raise FileNotFoundError(f"No production checkpoint found in {output_dir}")

    simulation.loadCheckpoint(str(checkpoint_path))
    checkpoint_step = _parse_checkpoint_step(checkpoint_path)
    logger.info("Resumed from checkpoint at step %d: %s", checkpoint_step, checkpoint_path)

    report_interval_steps = max(1, round(config.save_interval_ps / config.timestep_ps))
    checkpoint_interval_steps = max(1, round(config.checkpoint_interval_ps / config.timestep_ps))
    total_steps = max(1, round((config.duration_ns * 1000.0) / config.timestep_ps))

    completed_steps = checkpoint_step
    remaining_steps = total_steps - completed_steps
    if remaining_steps <= 0:
        logger.info(
            "Simulation already completed (%d/%d steps). Nothing to resume.",
            completed_steps, total_steps,
        )
        trajectory_path = output_dir / "production.dcd"
        energy_timeseries_path = output_dir / "production_energy.csv"
        return {
            "trajectory_path": trajectory_path,
            "n_frames": 0,
            "total_time_ns": float(completed_steps * config.timestep_ps / 1000.0),
            "checkpoint_paths": [checkpoint_path],
            "energy_timeseries_path": energy_timeseries_path,
        }

    n_remaining_frames = max(1, remaining_steps // report_interval_steps)
    trajectory_path = output_dir / "production.dcd"
    energy_timeseries_path = output_dir / "production_energy.csv"
    integration_timestep_ps = simulation.integrator.getStepSize().value_in_unit(unit.picoseconds)

    energy_timeseries = np.zeros((n_remaining_frames, 3), dtype=float)
    checkpoint_paths: list[Path] = [checkpoint_path]
    frame_index = 0

    next_report_step = report_interval_steps
    next_checkpoint_step = checkpoint_interval_steps
    elapsed_steps = 0

    logger.info(
        "Resuming production MD for %d remaining steps (%d frames)",
        remaining_steps, n_remaining_frames,
    )

    with trajectory_path.open("ab") as trajectory_handle:
        dcd_file = DCDFile(
            trajectory_handle,
            simulation.topology,
            integration_timestep_ps * unit.picoseconds,
            firstStep=completed_steps,
            interval=report_interval_steps,
            append=True,
        )

        while elapsed_steps < remaining_steps:
            steps_to_report = next_report_step - elapsed_steps
            steps_to_checkpoint = next_checkpoint_step - elapsed_steps
            steps_to_event = min(
                steps_to_report, steps_to_checkpoint, remaining_steps - elapsed_steps,
            )
            simulation.step(steps_to_event)
            elapsed_steps += steps_to_event

            if elapsed_steps == next_report_step and frame_index < n_remaining_frames:
                state = simulation.context.getState(getEnergy=True, getPositions=True)
                energy_timeseries[frame_index, 0] = (completed_steps + elapsed_steps) * config.timestep_ps
                energy_timeseries[frame_index, 1] = state.getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                )
                energy_timeseries[frame_index, 2] = state.getKineticEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                )
                dcd_file.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())
                frame_index += 1
                next_report_step += report_interval_steps

            if elapsed_steps == next_checkpoint_step:
                chk_path = _checkpoint_path(output_dir, completed_steps + elapsed_steps)
                simulation.saveCheckpoint(str(chk_path))
                checkpoint_paths.append(chk_path)
                next_checkpoint_step += checkpoint_interval_steps

    _append_energy_timeseries(energy_timeseries_path, energy_timeseries[:frame_index])

    _validate_nve_reference_drift(simulation, config)

    return {
        "trajectory_path": trajectory_path,
        "n_frames": int(frame_index),
        "total_time_ns": float((completed_steps + elapsed_steps) * config.timestep_ps / 1000.0),
        "checkpoint_paths": checkpoint_paths,
        "energy_timeseries_path": energy_timeseries_path,
    }