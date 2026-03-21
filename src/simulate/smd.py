"""Steered molecular dynamics utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import csv
import logging
import os
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import openmm
from openmm import XmlSerializer, unit
from openmm.app import DCDFile, Element, Simulation, Topology

from src import PhysicalValidityError
from src.config import ProductionConfig, SMDConfig
from src.physics.collective_variables import com_distance, com_vector, CollectiveVariableSpec, compute_cv
from src.simulate.platform import select_platform


logger = logging.getLogger(__name__)


def _resolve_seed(seed: int | None, stage_name: str) -> int:
    """Resolve a random seed, generating one from OS entropy if None."""
    if seed is None:
        seed = int.from_bytes(os.urandom(4), "big") & 0x7FFFFFFF
    logger.info("Using base random seed %d for %s", seed, stage_name)
    return seed


def _validate_config(config: SMDConfig) -> None:
    """Validate public-boundary SMD configuration values."""

    if config.spring_constant_kj_mol_nm2 <= 0.0:
        raise ValueError("config.spring_constant_kj_mol_nm2 must be positive")
    if config.pulling_velocity_nm_per_ps <= 0.0:
        raise ValueError("config.pulling_velocity_nm_per_ps must be positive")
    if config.pull_distance_nm <= 0.0:
        raise ValueError("config.pull_distance_nm must be positive")
    if config.n_replicates <= 0:
        raise ValueError("config.n_replicates must be positive")
    if config.save_interval_ps <= 0.0:
        raise ValueError("config.save_interval_ps must be positive")


def _validate_group_indices(group_indices: list[int], n_particles: int, parameter_name: str) -> list[int]:
    """Validate COM-group selections against the system size."""

    validated = [int(index) for index in group_indices]
    if not validated:
        raise ValueError(f"{parameter_name} must be non-empty")
    if any(index < 0 or index >= n_particles for index in validated):
        raise ValueError(f"{parameter_name} must be within [0, N_atoms)")
    return validated


def _validate_pull_direction(pull_direction: np.ndarray) -> np.ndarray:
    """Validate and normalize the pulling direction."""

    direction = np.asarray(pull_direction, dtype=float)
    if direction.shape != (3,):
        raise ValueError("pull_direction must have shape [3]")
    if not np.all(np.isfinite(direction)):
        raise ValueError("pull_direction must contain only finite values")
    norm = float(np.linalg.norm(direction))
    if norm <= 0.0:
        raise ValueError("pull_direction must have non-zero norm")
    return direction / norm


def _particle_masses_amu(system: openmm.System) -> np.ndarray:
    """Extract particle masses in amu."""

    return np.asarray(
        [system.getParticleMass(index).value_in_unit(unit.dalton) for index in range(system.getNumParticles())],
        dtype=float,
    )


def _positions_nm(simulation: Simulation) -> np.ndarray:
    """Extract current Cartesian coordinates in nm."""

    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    return np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)


def _initial_reaction_coordinate(
    simulation: Simulation,
    group_a_indices: list[int],
    group_b_indices: list[int],
    pull_direction: np.ndarray,
) -> float:
    """Compute the initial COM separation and validate the pull direction orientation."""

    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    positions = np.asarray(
        state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float
    )
    box_vectors = state.getPeriodicBoxVectors()
    box_lengths = np.array([
        box_vectors[0][0].value_in_unit(unit.nanometer),
        box_vectors[1][1].value_in_unit(unit.nanometer),
        box_vectors[2][2].value_in_unit(unit.nanometer),
    ])
    masses = _particle_masses_amu(simulation.system)
    com_a = com_vector(positions, masses, np.asarray(group_a_indices, dtype=int), box_lengths=box_lengths)
    com_b = com_vector(positions, masses, np.asarray(group_b_indices, dtype=int), box_lengths=box_lengths)
    com_axis = com_a - com_b
    com_axis = com_axis - box_lengths * np.round(com_axis / box_lengths)
    if float(np.dot(com_axis, pull_direction)) <= 0.0:
        raise ValueError("pull_direction must point from pull_group_2 toward pull_group_1 at t=0")
    return com_distance(
        positions,
        masses,
        np.asarray(group_a_indices, dtype=int),
        np.asarray(group_b_indices, dtype=int),
        box_lengths=box_lengths,
    )


def _create_smd_force(
    simulation: Simulation,
    group_a_indices: list[int],
    group_b_indices: list[int],
    spring_constant_kj_mol_nm2: float,
    initial_distance_nm: float,
) -> openmm.CustomCentroidBondForce:
    """Attach a time-dependent centroid-distance SMD force to the system."""

    smd_force = openmm.CustomCentroidBondForce(2, "0.5*smd_k*(distance(g1, g2)-smd_r_ref)^2")
    smd_force.addGlobalParameter("smd_k", float(spring_constant_kj_mol_nm2))
    smd_force.addGlobalParameter("smd_r_ref", float(initial_distance_nm))
    group_a = smd_force.addGroup(group_a_indices)
    group_b = smd_force.addGroup(group_b_indices)
    smd_force.addBond([group_a, group_b], [])
    simulation.system.addForce(smd_force)
    simulation.context.reinitialize(preserveState=True)
    return simulation.system.getForce(simulation.system.getNumForces() - 1)


def _write_timeseries_csv(output_path: Path, header: list[str], values: np.ndarray) -> None:
    """Persist a two-column SMD timeseries to CSV."""

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(values.tolist())


def _replicate_output_paths(output_dir: Path, replicate_id: int) -> dict[str, Path]:
    """Create replicate-local output paths."""

    replicate_dir = output_dir / f"replicate_{replicate_id:03d}"
    replicate_dir.mkdir(parents=True, exist_ok=True)
    return {
        "replicate_dir": replicate_dir,
        "trajectory_path": replicate_dir / f"smd_replicate_{replicate_id:03d}.dcd",
        "work_path": replicate_dir / f"smd_replicate_{replicate_id:03d}_work.csv",
        "force_path": replicate_dir / f"smd_replicate_{replicate_id:03d}_force.csv",
        "xi_path": replicate_dir / f"smd_replicate_{replicate_id:03d}_xi.csv",
    }


def _generic_topology(n_particles: int) -> Topology:
    """Construct a generic topology when only serialized state and system are available.

    .. deprecated::
        This function produces incorrect chain assignments and element metadata.
        Use _load_topology() with a topology_reference.pdb file from NPT equilibration.
    """

    warnings.warn(
        "_generic_topology() produces incorrect chain assignments. "
        "Provide topology_pdb_path to campaign functions instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("SMD", chain)
    carbon = Element.getByAtomicNumber(6)
    for atom_index in range(n_particles):
        topology.addAtom(f"A{atom_index}", carbon, residue)
    return topology


def _load_topology(topology_pdb_path: Path | None, n_particles: int) -> Topology:
    """Load the authentic topology from a PDB file, falling back to generic if not provided.

    Args:
        topology_pdb_path: Path to the topology reference PDB. If None, falls back
            to _generic_topology() with a deprecation warning.
        n_particles: Expected number of particles for validation.

    Returns:
        OpenMM Topology object.
    """
    if topology_pdb_path is None:
        logger.warning(
            "No topology_pdb_path provided; falling back to _generic_topology(). "
            "Chain assignments and element metadata will be incorrect. "
            "Re-run NPT equilibration to generate a topology_reference.pdb file."
        )
        return _generic_topology(n_particles)

    from openmm.app import PDBFile

    pdb = PDBFile(str(topology_pdb_path))
    loaded_n_atoms = pdb.topology.getNumAtoms()
    if loaded_n_atoms != n_particles:
        raise ValueError(
            f"Topology atom count ({loaded_n_atoms}) does not match "
            f"system particle count ({n_particles})"
        )
    return pdb.topology


def _configure_replicate_velocities(simulation: Simulation, replicate_id: int, base_seed: int) -> None:
    """Initialize replicate-specific velocities from a resolved base seed."""

    temperature_k = ProductionConfig().temperature_k
    random_seed = base_seed + int(replicate_id)
    if hasattr(simulation.integrator, "setRandomNumberSeed"):
        simulation.integrator.setRandomNumberSeed(random_seed)
    simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin, random_seed)


def _validate_work_unimodality(results: list[dict[str, Any]]) -> None:
    """Enforce IV-10 on campaign-level work values when enough replicates are available."""

    if len(results) < 10:
        logger.warning("Skipping IV-10 work-unimodality diagnostic: fewer than 10 SMD replicates")
        return

    total_work_values = np.asarray([result["total_work_kj_mol"] for result in results], dtype=float)
    grid = np.linspace(float(np.min(total_work_values)), float(np.max(total_work_values)), 256)
    if np.allclose(grid[0], grid[-1]):
        return

    bandwidth = float(np.std(total_work_values))
    if bandwidth == 0.0:
        return
    kernel = np.exp(-0.5 * ((grid[:, None] - total_work_values[None, :]) / bandwidth) ** 2)
    density = np.mean(kernel, axis=1)
    peaks = 0
    for index in range(1, density.size - 1):
        if density[index] > density[index - 1] and density[index] > density[index + 1] and density[index] > 0.05 * np.max(density):
            peaks += 1
    if peaks > 1:
        raise PhysicalValidityError("IV-10 violated: SMD work distribution is multimodal")


def _validate_trajectory_chain_integrity(
    topology: Topology,
    expected_min_chains: int = 2,
) -> None:
    """Verify that the trajectory topology contains the expected chain structure.

    Args:
        topology: The OpenMM Topology used to write the trajectory.
        expected_min_chains: Minimum expected number of protein chains.

    Raises:
        PhysicalValidityError: If chain count is below the expected minimum.
    """
    chain_count = sum(1 for _ in topology.chains())
    if chain_count < expected_min_chains:
        raise PhysicalValidityError(
            f"Trajectory topology contains only {chain_count} chain(s), "
            f"expected at least {expected_min_chains}. "
            "This may indicate topology corruption from _generic_topology()."
        )


def run_smd_replicate(
    simulation: Simulation,
    config: SMDConfig,
    replicate_id: int,
    output_dir: Path,
    pull_group_1: list[int],
    pull_group_2: list[int],
    pull_direction: np.ndarray,
    auxiliary_cvs: list[CollectiveVariableSpec] | None = None,
) -> dict[str, Any]:
    """Run a single constant-velocity SMD trajectory along the COM reaction coordinate."""

    _validate_config(config)
    direction = _validate_pull_direction(pull_direction)
    n_particles = simulation.system.getNumParticles()
    group_a = _validate_group_indices(pull_group_1, n_particles, "pull_group_1")
    group_b = _validate_group_indices(pull_group_2, n_particles, "pull_group_2")
    if set(group_a).intersection(group_b):
        raise ValueError("pull_group_1 and pull_group_2 must be disjoint")

    output_paths = _replicate_output_paths(Path(output_dir), replicate_id)
    initial_distance_nm = _initial_reaction_coordinate(simulation, group_a, group_b, direction)
    target_final_distance_nm = initial_distance_nm + config.pull_distance_nm
    total_time_ps = config.pull_distance_nm / config.pulling_velocity_nm_per_ps
    timestep_ps = simulation.integrator.getStepSize().value_in_unit(unit.picoseconds)
    report_interval_steps = max(1, round(config.save_interval_ps / timestep_ps))
    n_samples = max(1, round(total_time_ps / config.save_interval_ps))
    total_time_ps = n_samples * report_interval_steps * timestep_ps
    smd_force = _create_smd_force(simulation, group_a, group_b, config.spring_constant_kj_mol_nm2, initial_distance_nm)

    masses = _particle_masses_amu(simulation.system)
    cumulative_work_kj_mol = 0.0  # Integrated W = integral F * v dt
    final_xi_nm = 0.0

    active_aux_cvs = auxiliary_cvs if auxiliary_cvs else []
    aux_timeseries: dict[str, np.ndarray] = {
        spec.name: np.zeros((n_samples, 2), dtype=float) for spec in active_aux_cvs
    }

    logger.info(
        "Starting SMD replicate %d for %.3f ps with target displacement %.3f nm",
        replicate_id,
        total_time_ps,
        config.pull_distance_nm,
    )

    with output_paths["trajectory_path"].open("wb") as trajectory_handle, \
         output_paths["work_path"].open("w", encoding="utf-8", newline="") as work_handle, \
         output_paths["force_path"].open("w", encoding="utf-8", newline="") as force_handle, \
         output_paths["xi_path"].open("w", encoding="utf-8", newline="") as xi_handle:

        work_writer = csv.writer(work_handle)
        work_writer.writerow(["time_ps", "work_kj_mol"])
        force_writer = csv.writer(force_handle)
        force_writer.writerow(["time_ps", "force_kj_mol_nm"])
        xi_writer = csv.writer(xi_handle)
        xi_writer.writerow(["time_ps", "xi_nm"])

        dcd_file = DCDFile(
            trajectory_handle,
            simulation.topology,
            timestep_ps * unit.picoseconds,
            firstStep=simulation.currentStep,
            interval=report_interval_steps,
            append=False,
        )

        for sample_index in range(n_samples):
            # Advance the reference distance at constant pulling velocity.
            elapsed_time_ps = (sample_index + 1) * report_interval_steps * timestep_ps
            reference_distance_nm = min(
                target_final_distance_nm,
                initial_distance_nm + config.pulling_velocity_nm_per_ps * elapsed_time_ps,
            )
            simulation.context.setParameter("smd_r_ref", reference_distance_nm)
            simulation.step(report_interval_steps)

            # Measure the instantaneous COM separation (reaction coordinate).
            state = simulation.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True)
            positions_nm = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)
            box_vectors = state.getPeriodicBoxVectors()
            box_lengths = np.array([
                box_vectors[0][0].value_in_unit(unit.nanometer),
                box_vectors[1][1].value_in_unit(unit.nanometer),
                box_vectors[2][2].value_in_unit(unit.nanometer),
            ])
            xi_nm = com_distance(
                positions_nm,
                masses,
                np.asarray(group_a, dtype=int),
                np.asarray(group_b, dtype=int),
                box_lengths=box_lengths,
            )

            # Harmonic restoring force: F = -k * (xi - xi_ref).
            force_kj_mol_nm = -config.spring_constant_kj_mol_nm2 * (xi_nm - reference_distance_nm)
            # Accumulate work via the trapezoidal rule: dW = F * v * dt.
            cumulative_work_kj_mol += force_kj_mol_nm * config.pulling_velocity_nm_per_ps * report_interval_steps * timestep_ps

            # Stream time series directly to CSV — O(1) memory per replicate.
            work_writer.writerow([elapsed_time_ps, cumulative_work_kj_mol])
            force_writer.writerow([elapsed_time_ps, force_kj_mol_nm])
            xi_writer.writerow([elapsed_time_ps, xi_nm])
            final_xi_nm = xi_nm

            for aux_spec in active_aux_cvs:
                aux_val = compute_cv(
                    aux_spec, positions_nm, masses,
                    np.asarray(group_a, dtype=int),
                    np.asarray(group_b, dtype=int),
                    box_lengths=box_lengths,
                )
                aux_timeseries[aux_spec.name][sample_index] = [elapsed_time_ps, aux_val]

            dcd_file.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())

    for aux_spec in active_aux_cvs:
        aux_csv_path = output_paths["replicate_dir"] / f"smd_replicate_{replicate_id:03d}_aux_{aux_spec.name}.csv"
        _write_timeseries_csv(aux_csv_path, ["time_ps", aux_spec.name], aux_timeseries[aux_spec.name])

    result: dict[str, Any] = {
        "trajectory_path": output_paths["trajectory_path"],
        "total_work_kj_mol": float(cumulative_work_kj_mol),
        "final_xi_nm": float(final_xi_nm),
        "n_samples": n_samples,
        "work_timeseries_path": output_paths["work_path"],
        "force_timeseries_path": output_paths["force_path"],
        "xi_timeseries_path": output_paths["xi_path"],
    }
    for aux_spec in active_aux_cvs:
        result[f"aux_cv_{aux_spec.name}_timeseries"] = aux_timeseries[aux_spec.name]
    return result


# L-30 Step 1: Process-safe SMD worker for multi-process parallelism.
def _run_smd_worker(
    system_xml: str,
    state_xml: str,
    config: SMDConfig,
    replicate_id: int,
    output_dir: Path,
    pull_group_1: list[int],
    pull_group_2: list[int],
    pull_direction: np.ndarray,
    base_seed: int,
    temperature_k: float,
    friction_per_ps: float,
    timestep_ps: float,
    platform_name: str | None,
    topology_pdb_path: Path | None = None,
    auxiliary_cvs: list[CollectiveVariableSpec] | None = None,
) -> dict[str, Any]:
    """Process-safe SMD replicate worker.

    Reconstructs the OpenMM Simulation from serialized XML strings so that
    this function can be dispatched to a ``ProcessPoolExecutor`` without
    passing unpicklable OpenMM objects across process boundaries.
    """
    system = XmlSerializer.deserialize(system_xml)
    state = XmlSerializer.deserialize(state_xml)
    topology = _load_topology(topology_pdb_path, system.getNumParticles())
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        friction_per_ps / unit.picosecond,
        timestep_ps * unit.picoseconds,
    )
    simulation = Simulation(
        topology, system, integrator,
        select_platform(platform_name),
    )
    simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    simulation.context.setPositions(state.getPositions())
    _configure_replicate_velocities(simulation, replicate_id, base_seed)
    return run_smd_replicate(
        simulation, config, replicate_id,
        Path(output_dir), pull_group_1, pull_group_2, pull_direction,
        auxiliary_cvs=auxiliary_cvs,
    )


# L-30 Step 2: n_workers parameter for parallel SMD campaign execution.
def run_smd_campaign(
    equilibrated_state_path: Path,
    system_xml_path: Path,
    config: SMDConfig,
    pull_group_1: list[int],
    pull_group_2: list[int],
    output_dir: Path,
    auxiliary_cvs: list[CollectiveVariableSpec] | None = None,
    topology_pdb_path: Path | None = None,
    platform_name: str | None = None,
    n_workers: int = 1,
) -> list[dict[str, Any]]:
    """Run a set of SMD pulling replicates from a shared equilibrated starting state."""

    if n_workers < 1:
        raise ValueError("n_workers must be at least 1")

    _validate_config(config)
    system_xml = Path(system_xml_path).read_text(encoding="utf-8")
    state_xml = Path(equilibrated_state_path).read_text(encoding="utf-8")
    system = XmlSerializer.deserialize(system_xml)
    state = XmlSerializer.deserialize(state_xml)
    topology = _load_topology(topology_pdb_path, system.getNumParticles())
    production_defaults = ProductionConfig()
    # Derive the initial pull direction from the equilibrated COM axis.
    direction_positions = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)
    box_vectors = state.getPeriodicBoxVectors()
    box_lengths = np.array([
        box_vectors[0][0].value_in_unit(unit.nanometer),
        box_vectors[1][1].value_in_unit(unit.nanometer),
        box_vectors[2][2].value_in_unit(unit.nanometer),
    ])
    masses = _particle_masses_amu(system)
    com_a = com_vector(direction_positions, masses, np.asarray(pull_group_1, dtype=int), box_lengths=box_lengths)
    com_b = com_vector(direction_positions, masses, np.asarray(pull_group_2, dtype=int), box_lengths=box_lengths)
    pull_direction = com_a - com_b
    pull_direction = pull_direction - box_lengths * np.round(pull_direction / box_lengths)
    pull_direction = pull_direction / np.linalg.norm(pull_direction)

    base_seed = _resolve_seed(config.random_seed, "SMD campaign")

    if n_workers == 1:
        # Sequential path — zero overhead, backward compatible.
        results: list[dict[str, Any]] = []
        for replicate_id in range(1, config.n_replicates + 1):
            replicate_system = XmlSerializer.deserialize(XmlSerializer.serialize(system))
            integrator = openmm.LangevinMiddleIntegrator(
                production_defaults.temperature_k * unit.kelvin,
                production_defaults.friction_per_ps / unit.picosecond,
                production_defaults.timestep_ps * unit.picoseconds,
            )
            simulation = Simulation(topology, replicate_system, integrator, select_platform(platform_name))
            logger.info("Simulation platform: %s", simulation.context.getPlatform().getName())
            simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
            simulation.context.setPositions(state.getPositions())
            _configure_replicate_velocities(simulation, replicate_id, base_seed)
            results.append(
                run_smd_replicate(
                    simulation,
                    config,
                    replicate_id,
                    Path(output_dir),
                    pull_group_1,
                    pull_group_2,
                    pull_direction,
                    auxiliary_cvs=auxiliary_cvs,
                )
            )
    else:
        # Parallel path — dispatch replicates to a process pool.
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(
                    _run_smd_worker,
                    system_xml=system_xml,
                    state_xml=state_xml,
                    config=config,
                    replicate_id=replicate_id,
                    output_dir=output_dir,
                    pull_group_1=pull_group_1,
                    pull_group_2=pull_group_2,
                    pull_direction=pull_direction,
                    base_seed=base_seed,
                    temperature_k=production_defaults.temperature_k,
                    friction_per_ps=production_defaults.friction_per_ps,
                    timestep_ps=production_defaults.timestep_ps,
                    platform_name=platform_name,
                    topology_pdb_path=topology_pdb_path,
                    auxiliary_cvs=auxiliary_cvs,
                ): replicate_id
                for replicate_id in range(1, config.n_replicates + 1)
            }
            results_map: dict[int, dict[str, Any]] = {}
            for future in as_completed(futures):
                rid = futures[future]
                results_map[rid] = future.result()
            results = [results_map[rid] for rid in sorted(results_map)]

    _validate_work_unimodality(results)
    if topology_pdb_path is not None:
        _validate_trajectory_chain_integrity(topology)
    return results