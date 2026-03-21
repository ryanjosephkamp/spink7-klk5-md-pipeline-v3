"""Umbrella-sampling utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import json
import logging
import math
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import openmm
from openmm import XmlSerializer, unit
from openmm.app import DCDFile, Element, Simulation, Topology

from src import PhysicalValidityError
from src.config import ProductionConfig, UmbrellaConfig
from src.physics.collective_variables import com_distance, com_vector, CollectiveVariableSpec, compute_cv
from src.simulate.platform import select_platform


logger = logging.getLogger(__name__)


def _validate_config(config: UmbrellaConfig) -> None:
    """Validate public-boundary umbrella-sampling parameters."""

    if config.xi_min_nm <= 0.0:
        raise ValueError("config.xi_min_nm must be positive")
    if config.xi_max_nm <= config.xi_min_nm:
        raise ValueError("config.xi_max_nm must be greater than config.xi_min_nm")
    if config.window_spacing_nm <= 0.0:
        raise ValueError("config.window_spacing_nm must be positive")
    if config.spring_constant_kj_mol_nm2 <= 0.0:
        raise ValueError("config.spring_constant_kj_mol_nm2 must be positive")
    if config.per_window_duration_ns <= 0.0:
        raise ValueError("config.per_window_duration_ns must be positive")
    if config.save_interval_ps <= 0.0:
        raise ValueError("config.save_interval_ps must be positive")


_MANIFEST_FILENAME = "umbrella_manifest.json"


def _load_manifest(output_dir: Path) -> dict:
    """Load the umbrella campaign manifest, or return an empty manifest."""

    manifest_path = output_dir / _MANIFEST_FILENAME
    if manifest_path.exists():
        with manifest_path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    return {"completed_windows": {}, "version": 1}


def _save_manifest(output_dir: Path, manifest: dict) -> None:
    """Persist the umbrella campaign manifest to disk."""

    manifest_path = output_dir / _MANIFEST_FILENAME
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


def generate_window_centers(config: UmbrellaConfig) -> np.ndarray:
    """Generate umbrella window centers spanning the configured reaction-coordinate range.

    Invariants: None.

    Args:
        config: Umbrella sampling parameters.

    Returns:
        np.ndarray: Window-center positions in nm. Shape: [M].
    """

    _validate_config(config)
    n_windows = int(round((config.xi_max_nm - config.xi_min_nm) / config.window_spacing_nm)) + 1
    return config.xi_min_nm + np.arange(n_windows, dtype=float) * config.window_spacing_nm


def _validate_group_indices(group_indices: list[int], n_particles: int, parameter_name: str) -> list[int]:
    """Validate COM-group selections against the system size."""

    validated = [int(index) for index in group_indices]
    if not validated:
        raise ValueError(f"{parameter_name} must be non-empty")
    if any(index < 0 or index >= n_particles for index in validated):
        raise ValueError(f"{parameter_name} must be within [0, N_atoms)")
    return validated


def _particle_masses_amu(system: openmm.System) -> np.ndarray:
    """Extract particle masses in amu."""

    return np.asarray(
        [system.getParticleMass(index).value_in_unit(unit.dalton) for index in range(system.getNumParticles())],
        dtype=float,
    )


def _positions_nm(simulation: Simulation) -> np.ndarray:
    """Extract the current Cartesian coordinates in nm."""

    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    return np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)


def _infer_pull_groups_from_topology(simulation: Simulation) -> tuple[list[int], list[int]]:
    """Infer COM groups from the first two non-empty topology chains."""

    chain_groups: list[list[int]] = []
    for chain in simulation.topology.chains():
        atom_indices = [atom.index for atom in chain.atoms()]
        if atom_indices:
            chain_groups.append(atom_indices)
    if len(chain_groups) < 2:
        raise ValueError("simulation topology must contain at least two non-empty chains to infer umbrella groups")
    return chain_groups[0], chain_groups[1]


def _create_umbrella_force(
    simulation: Simulation,
    group_a_indices: list[int],
    group_b_indices: list[int],
    window_center_nm: float,
    spring_constant_kj_mol_nm2: float,
) -> openmm.CustomCentroidBondForce:
    """Attach a harmonic centroid-distance umbrella restraint to the system."""

    umbrella_force = openmm.CustomCentroidBondForce(2, "0.5*umbrella_k*(distance(g1, g2)-umbrella_r_ref)^2")
    umbrella_force.addGlobalParameter("umbrella_k", float(spring_constant_kj_mol_nm2))
    umbrella_force.addGlobalParameter("umbrella_r_ref", float(window_center_nm))
    group_a = umbrella_force.addGroup(group_a_indices)
    group_b = umbrella_force.addGroup(group_b_indices)
    umbrella_force.addBond([group_a, group_b], [])
    simulation.system.addForce(umbrella_force)
    simulation.context.reinitialize(preserveState=True)
    return simulation.system.getForce(simulation.system.getNumForces() - 1)


def _window_output_paths(output_dir: Path, window_id: int, window_center_nm: float) -> dict[str, Path]:
    """Construct deterministic file paths for a single umbrella window."""

    output_dir.mkdir(parents=True, exist_ok=True)
    center_label = f"{window_center_nm:.2f}".replace("-", "m")
    return {
        "trajectory_path": output_dir / f"umbrella_window_{window_id:03d}_xi{center_label}nm.dcd",
        "xi_timeseries_path": output_dir / f"umbrella_window_{window_id:03d}_xi{center_label}nm.npy",
    }


def _dummy_topology(n_particles: int) -> Topology:
    """Build a minimal two-chain topology for the Simulation / DCD writer."""

    topology = Topology()
    chain_a = topology.addChain()
    chain_b = topology.addChain()
    residue_a = topology.addResidue("UMB_A", chain_a)
    residue_b = topology.addResidue("UMB_B", chain_b)
    carbon = Element.getByAtomicNumber(6)
    midpoint = n_particles // 2
    for atom_index in range(n_particles):
        residue = residue_a if atom_index < midpoint else residue_b
        topology.addAtom(f"A{atom_index}", carbon, residue)
    return topology


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

    return _dummy_topology(n_particles)


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
            "Chain assignments and pull groups may be incorrect. "
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


def _pre_position_to_target(
    simulation: Simulation,
    group_a: list[int],
    group_b: list[int],
    target_xi_nm: float,
    masses: np.ndarray,
    config: UmbrellaConfig,
) -> float:
    """Steer the system COM distance to the target window center.

    Uses a time-varying harmonic bias to smoothly pull the system from its
    current COM distance to *target_xi_nm*.  The bias force is removed after
    pre-positioning so that production sampling is unaffected.

    Returns the final COM distance after pre-positioning (nm).
    """

    state = simulation.context.getState(getPositions=True)
    positions = np.asarray(
        state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float,
    )
    current_xi = com_distance(
        positions, masses, np.asarray(group_a, dtype=int), np.asarray(group_b, dtype=int),
    )

    if abs(current_xi - target_xi_nm) < 0.01:
        return current_xi

    # Attach a temporary centroid-bond pulling force.
    pre_force = openmm.CustomCentroidBondForce(
        2, "0.5*pre_k*(distance(g1,g2)-pre_r_ref)^2",
    )
    pre_force.addGlobalParameter("pre_k", float(config.pre_position_spring_constant_kj_mol_nm2))
    pre_force.addGlobalParameter("pre_r_ref", float(current_xi))
    ga = pre_force.addGroup(group_a)
    gb = pre_force.addGroup(group_b)
    pre_force.addBond([ga, gb], [])
    force_index = simulation.system.addForce(pre_force)
    simulation.context.reinitialize(preserveState=True)

    timestep_ps = simulation.integrator.getStepSize().value_in_unit(unit.picoseconds)
    velocity = config.pre_position_velocity_nm_per_ps  # nm/ps
    distance_to_cover = abs(target_xi_nm - current_xi)
    n_steps = max(1, math.ceil(distance_to_cover / (velocity * timestep_ps)))
    direction = 1.0 if target_xi_nm > current_xi else -1.0

    for step_i in range(1, n_steps + 1):
        r_ref = current_xi + direction * velocity * timestep_ps * step_i
        # Clamp so we don't overshoot.
        if direction > 0:
            r_ref = min(r_ref, target_xi_nm)
        else:
            r_ref = max(r_ref, target_xi_nm)
        simulation.context.setParameter("pre_r_ref", float(r_ref))
        simulation.step(1)

    # Remove the pre-positioning force.
    simulation.system.removeForce(force_index)
    simulation.context.reinitialize(preserveState=True)

    final_state = simulation.context.getState(getPositions=True)
    final_positions = np.asarray(
        final_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float,
    )
    final_xi = com_distance(
        final_positions, masses, np.asarray(group_a, dtype=int), np.asarray(group_b, dtype=int),
    )
    logger.info(
        "Pre-positioning: xi %.4f -> %.4f nm (target %.4f nm, %d steps)",
        current_xi, final_xi, target_xi_nm, n_steps,
    )
    return float(final_xi)


def _integrated_autocorrelation_time(series: np.ndarray) -> float:
    """Estimate integrated autocorrelation time using the initial positive sequence estimator.

    Uses Geyer's (1992) method: sum consecutive pairs of autocovariance values
    and truncate the sum at the first non-positive pair to prevent noise from
    inflating the estimate.

    Args:
        series: 1-D stationary timeseries. Must contain at least 3 elements.

    Returns:
        Integrated autocorrelation time (in units of the sampling interval).
    """

    n = len(series)
    if n < 3:
        return 0.0
    mean = np.mean(series)
    var = np.var(series, ddof=0)
    if var < 1e-30:
        return 0.0

    # Normalized autocovariance function.
    centered = series - mean
    tau_int = 0.0
    # Process lag pairs (1,2), (3,4), ... per Geyer's initial positive sequence.
    lag = 1
    while lag < n - 1:
        c_lag = float(np.mean(centered[: n - lag] * centered[lag:]))
        c_lag1 = float(np.mean(centered[: n - lag - 1] * centered[lag + 1 :])) if lag + 1 < n else 0.0
        pair_sum = (c_lag + c_lag1) / var
        if pair_sum <= 0.0:
            break
        tau_int += pair_sum
        lag += 2

    return tau_int


def _detect_and_trim_equilibration(xi_timeseries: np.ndarray) -> tuple[np.ndarray, int]:
    """Detect and discard the equilibration transient from a 1D timeseries.

    Uses the maximum effective sample size criterion (Chodera, JCTC 2016).

    If ``pymbar`` is installed, delegates to ``pymbar.timeseries.detect_equilibration``.
    Otherwise, uses a built-in implementation of the same algorithm.

    Args:
        xi_timeseries: Full production timeseries. Shape: [T].

    Returns:
        Tuple of (equilibrated timeseries, number of discarded frames).
    """

    try:
        from pymbar.timeseries import detect_equilibration  # type: ignore[import-untyped]
        t0, g, n_eff = detect_equilibration(xi_timeseries)
        return xi_timeseries[t0:], int(t0)
    except ImportError:
        pass

    n = len(xi_timeseries)
    min_tail = 10
    if n <= min_tail:
        return xi_timeseries, 0

    best_t0 = 0
    best_n_eff = 0.0

    for t0 in range(n - min_tail):
        tail = xi_timeseries[t0:]
        tau = _integrated_autocorrelation_time(tail)
        g = 1.0 + 2.0 * tau
        n_eff = len(tail) / g
        if n_eff > best_n_eff:
            best_n_eff = n_eff
            best_t0 = t0

    return xi_timeseries[best_t0:], best_t0


def _run_umbrella_window_with_groups(
    simulation: Simulation,
    window_center_nm: float,
    config: UmbrellaConfig,
    window_id: int,
    output_dir: Path,
    pull_group_1: list[int],
    pull_group_2: list[int],
    auxiliary_cvs: list[CollectiveVariableSpec] | None = None,
) -> dict[str, Any]:
    """Run a single restrained umbrella window with explicit COM groups.

    Invariants: None directly. Campaign-level overlap is checked separately under IV-8.

    Args:
        simulation: OpenMM simulation object.
        window_center_nm: Target COM distance in nm.
        config: Umbrella sampling parameters.
        window_id: Integer window identifier.
        output_dir: Directory for on-disk outputs.
        pull_group_1: First COM group. Shape: [N_a].
        pull_group_2: Second COM group. Shape: [N_b].

    Returns:
        dict[str, Any]: Window metadata, xi timeseries, and streamed trajectory path.
    """

    _validate_config(config)
    if window_center_nm <= 0.0:
        raise ValueError("window_center_nm must be positive")
    n_particles = simulation.system.getNumParticles()
    group_a = _validate_group_indices(pull_group_1, n_particles, "pull_group_1")
    group_b = _validate_group_indices(pull_group_2, n_particles, "pull_group_2")
    if set(group_a).intersection(group_b):
        raise ValueError("pull_group_1 and pull_group_2 must be disjoint")

    masses = _particle_masses_amu(simulation.system)

    # Phase 1: Steered pre-positioning to the target window center.
    _pre_position_to_target(
        simulation, group_a, group_b, window_center_nm, masses, config,
    )

    # Phase 2: Attach umbrella bias at the target window center.
    _create_umbrella_force(simulation, group_a, group_b, window_center_nm, config.spring_constant_kj_mol_nm2)

    # Phase 3: Biased equilibration (no frames recorded).
    timestep_ps = simulation.integrator.getStepSize().value_in_unit(unit.picoseconds)
    n_equil_steps = max(1, math.ceil(config.equilibration_duration_ps / timestep_ps))
    simulation.step(n_equil_steps)
    logger.info(
        "Window %d: completed %d-step equilibration (%.1f ps)",
        window_id, n_equil_steps, config.equilibration_duration_ps,
    )

    # Phase 4: Production sampling.
    output_paths = _window_output_paths(Path(output_dir), window_id, window_center_nm)
    report_interval_steps = max(1, round(config.save_interval_ps / timestep_ps))
    n_samples = max(1, round((config.per_window_duration_ns * 1000.0) / config.save_interval_ps))
    xi_timeseries = np.zeros(n_samples, dtype=float)

    active_aux_cvs = auxiliary_cvs if auxiliary_cvs else []
    aux_timeseries: dict[str, np.ndarray] = {
        spec.name: np.zeros(n_samples, dtype=float) for spec in active_aux_cvs
    }

    logger.info(
        "Starting umbrella window %d at center %.3f nm for %d samples",
        window_id,
        window_center_nm,
        n_samples,
    )

    with output_paths["trajectory_path"].open("wb") as trajectory_handle:
        dcd_file = DCDFile(
            trajectory_handle,
            simulation.topology,
            timestep_ps * unit.picoseconds,
            firstStep=simulation.currentStep,
            interval=report_interval_steps,
            append=False,
        )

        for sample_index in range(n_samples):
            simulation.step(report_interval_steps)
            state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
            positions_nm = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)
            box_vectors = state.getPeriodicBoxVectors()
            box_lengths = np.array([
                box_vectors[0][0].value_in_unit(unit.nanometer),
                box_vectors[1][1].value_in_unit(unit.nanometer),
                box_vectors[2][2].value_in_unit(unit.nanometer),
            ])
            xi_timeseries[sample_index] = com_distance(
                positions_nm,
                masses,
                np.asarray(group_a, dtype=int),
                np.asarray(group_b, dtype=int),
                box_lengths=box_lengths,
            )

            for aux_spec in active_aux_cvs:
                aux_timeseries[aux_spec.name][sample_index] = compute_cv(
                    aux_spec, positions_nm, masses,
                    np.asarray(group_a, dtype=int),
                    np.asarray(group_b, dtype=int),
                    box_lengths=box_lengths,
                )

            dcd_file.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())

    np.save(output_paths["xi_timeseries_path"], xi_timeseries)

    for aux_spec in active_aux_cvs:
        aux_npy_path = output_dir / f"umbrella_window_{window_id:03d}_aux_{aux_spec.name}.npy"
        np.save(aux_npy_path, aux_timeseries[aux_spec.name])

    n_discarded = 0
    if config.detect_equilibration:
        xi_timeseries_eq, n_discarded = _detect_and_trim_equilibration(xi_timeseries)
        logger.info(
            "Window %d: equilibration detection discarded %d/%d frames",
            window_id, n_discarded, len(xi_timeseries),
        )
        xi_timeseries = xi_timeseries_eq
        for aux_spec in active_aux_cvs:
            aux_timeseries[aux_spec.name] = aux_timeseries[aux_spec.name][n_discarded:]

    result: dict[str, Any] = {
        "window_id": int(window_id),
        "window_center_nm": float(window_center_nm),
        "xi_timeseries": xi_timeseries,
        "trajectory_path": output_paths["trajectory_path"],
        "mean_xi_nm": float(np.mean(xi_timeseries)),
        "std_xi_nm": float(np.std(xi_timeseries)),
        "xi_timeseries_path": output_paths["xi_timeseries_path"],
        "n_discarded_frames": n_discarded,
    }
    for aux_spec in active_aux_cvs:
        result[f"aux_cv_{aux_spec.name}_timeseries"] = aux_timeseries[aux_spec.name]
    return result


# L-30 Step 3: Process-safe umbrella worker for multi-process parallelism.
def _run_umbrella_worker(
    system_xml: str,
    state_xml: str,
    config: UmbrellaConfig,
    window_id: int,
    window_center_nm: float,
    output_dir: Path,
    pull_group_1: list[int],
    pull_group_2: list[int],
    temperature_k: float,
    friction_per_ps: float,
    timestep_ps: float,
    platform_name: str | None,
    topology_pdb_path: Path | None = None,
    auxiliary_cvs: list[CollectiveVariableSpec] | None = None,
) -> dict[str, Any]:
    """Process-safe umbrella window worker.

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
    simulation.context.setVelocities(state.getVelocities())
    return _run_umbrella_window_with_groups(
        simulation, float(window_center_nm), config,
        window_id, Path(output_dir), pull_group_1, pull_group_2,
        auxiliary_cvs=auxiliary_cvs,
    )


def run_umbrella_window(
    simulation: Simulation,
    window_center_nm: float,
    config: UmbrellaConfig,
    window_id: int,
    output_dir: Path,
    auxiliary_cvs: list[CollectiveVariableSpec] | None = None,
) -> dict[str, Any]:
    """Run a single umbrella-sampling window using groups inferred from topology chains.

    Invariants: None directly. Campaign-level overlap is checked separately under IV-8.

    Args:
        simulation: OpenMM simulation object whose first two chains define the COM groups.
        window_center_nm: Target COM distance in nm.
        config: Umbrella sampling parameters.
        window_id: Integer window identifier.
        output_dir: Directory for on-disk outputs.
        auxiliary_cvs: Optional list of auxiliary CVs to record (observational only).

    Returns:
        dict[str, Any]: Window metadata and xi timeseries. `xi_timeseries` shape: [N_samples].
    """

    group_a, group_b = _infer_pull_groups_from_topology(simulation)
    return _run_umbrella_window_with_groups(simulation, window_center_nm, config, window_id, output_dir, group_a, group_b, auxiliary_cvs=auxiliary_cvs)


def _histogram_overlap_fraction(first_values: np.ndarray, second_values: np.ndarray) -> float:
    """Compute normalized histogram overlap between two xi distributions."""

    if first_values.size == 0 or second_values.size == 0:
        return 0.0
    lower = float(min(np.min(first_values), np.min(second_values)))
    upper = float(max(np.max(first_values), np.max(second_values)))
    if np.isclose(lower, upper):
        return 1.0
    bins = np.linspace(lower, upper, 32)
    first_hist, _ = np.histogram(first_values, bins=bins, density=True)
    second_hist, _ = np.histogram(second_values, bins=bins, density=True)
    bin_width = float(bins[1] - bins[0])
    return float(np.sum(np.minimum(first_hist, second_hist)) * bin_width)


def diagnose_histogram_coverage(
    xi_timeseries_list: list[np.ndarray],
    window_centers: np.ndarray,
    n_bins: int = 100,
) -> dict[str, Any]:
    """Compute the global histogram envelope and identify coverage holes.

    Args:
        xi_timeseries_list: List of per-window xi timeseries. Each element shape: [N_i].
        window_centers: Target center of each window in nm. Shape: [M].
        n_bins: Number of bins for the global histogram grid.

    Returns:
        dict with keys:
            - 'bin_edges': np.ndarray, shape [n_bins + 1]
            - 'bin_centers': np.ndarray, shape [n_bins]
            - 'envelope_counts': np.ndarray, shape [n_bins] — summed histogram counts
            - 'zero_count_bins': np.ndarray — bin center values where envelope is zero
            - 'has_coverage_holes': bool
            - 'coverage_fraction': float — fraction of bins with nonzero counts
    """

    all_values = np.concatenate(xi_timeseries_list)
    xi_min = float(np.min(all_values))
    xi_max = float(np.max(all_values))
    if np.isclose(xi_min, xi_max):
        xi_min -= 1e-3
        xi_max += 1e-3
    bin_edges = np.linspace(xi_min, xi_max, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    envelope = np.zeros(n_bins, dtype=float)
    for xi in xi_timeseries_list:
        counts, _ = np.histogram(xi, bins=bin_edges)
        envelope += counts

    zero_mask = envelope == 0
    return {
        "bin_edges": bin_edges,
        "bin_centers": bin_centers,
        "envelope_counts": envelope,
        "zero_count_bins": bin_centers[zero_mask],
        "has_coverage_holes": bool(np.any(zero_mask)),
        "coverage_fraction": float(np.sum(~zero_mask) / n_bins),
    }


def compute_overlap_matrix(
    xi_timeseries_list: list[np.ndarray],
) -> np.ndarray:
    """Compute the symmetric pairwise histogram overlap matrix for all window pairs.

    Args:
        xi_timeseries_list: List of per-window xi timeseries. Each element shape: [N_i].

    Returns:
        np.ndarray: Symmetric overlap matrix of shape [M, M]. Entry (i, j) is the
        histogram intersection O_{ij} in [0, 1]. Diagonal entries are 1.0.
    """

    m = len(xi_timeseries_list)
    overlap = np.eye(m, dtype=float)
    for i in range(m):
        for j in range(i + 1, m):
            o_ij = _histogram_overlap_fraction(xi_timeseries_list[i], xi_timeseries_list[j])
            overlap[i, j] = o_ij
            overlap[j, i] = o_ij
    return overlap


def compute_effective_sample_sizes(
    xi_timeseries_list: list[np.ndarray],
) -> np.ndarray:
    """Compute the effective sample size for each window timeseries.

    Uses the initial positive sequence estimator for the integrated
    autocorrelation time (Geyer, 1992).

    Args:
        xi_timeseries_list: List of per-window xi timeseries. Each element shape: [N_i].

    Returns:
        np.ndarray: Effective sample sizes. Shape: [M].
    """

    n_eff = np.empty(len(xi_timeseries_list), dtype=float)
    for i, xi in enumerate(xi_timeseries_list):
        tau = _integrated_autocorrelation_time(xi)
        g = max(1.0, 1.0 + 2.0 * tau)
        n_eff[i] = len(xi) / g
    return n_eff


# L-30 Step 4: n_workers parameter for parallel umbrella campaign execution.
def run_umbrella_campaign(
    equilibrated_state_path: Path,
    system_xml_path: Path,
    config: UmbrellaConfig,
    pull_group_1: list[int] | None = None,
    pull_group_2: list[int] | None = None,
    output_dir: Path = Path("umbrella_output"),
    auxiliary_cvs: list[CollectiveVariableSpec] | None = None,
    pdb_path: Path | None = None,
    topology_pdb_path: Path | None = None,
    platform_name: str | None = None,
    n_workers: int = 1,
) -> list[dict[str, Any]]:
    """Run an independent restrained simulation at each umbrella window center.

    Invariants: Enforces IV-8 by requiring adjacent xi histograms to overlap by at least 10%.

    Args:
        equilibrated_state_path: Serialized OpenMM state path.
        system_xml_path: Serialized OpenMM system path.
        config: Umbrella sampling parameters.
        pull_group_1: First COM group, or None to infer from topology.
        pull_group_2: Second COM group, or None to infer from topology.
        output_dir: Directory for per-window outputs.
        auxiliary_cvs: Optional list of auxiliary CVs to record.
        pdb_path: PDB file path for chain-aware group inference when
            pull_group_1/pull_group_2 are not provided.
        topology_pdb_path: Path to the topology reference PDB from NPT
            equilibration. Used for authentic DCD metadata when explicit
            pull groups are provided.
        n_workers: Number of parallel worker processes. 1 = sequential (default).

    Returns:
        list[dict[str, Any]]: Per-window umbrella results.
    """

    if n_workers < 1:
        raise ValueError("n_workers must be at least 1")

    _validate_config(config)
    system_xml = Path(system_xml_path).read_text(encoding="utf-8")
    state_xml = Path(equilibrated_state_path).read_text(encoding="utf-8")
    system = XmlSerializer.deserialize(system_xml)
    state = XmlSerializer.deserialize(state_xml)

    # Resolve topology_pdb_path from either argument (topology_pdb_path takes precedence).
    effective_topo_path = topology_pdb_path or pdb_path

    # Determine pull groups and topology.
    if pull_group_1 is not None and pull_group_2 is not None:
        topology = _load_topology(effective_topo_path, system.getNumParticles())
        effective_group_1 = pull_group_1
        effective_group_2 = pull_group_2
    elif effective_topo_path is not None:
        topology = _load_topology(effective_topo_path, system.getNumParticles())
        # Infer groups from the PDB topology chains.
        chain_groups: list[list[int]] = []
        for chain in topology.chains():
            atom_indices = [atom.index for atom in chain.atoms()]
            if atom_indices:
                chain_groups.append(atom_indices)
        if len(chain_groups) < 2:
            raise ValueError(
                "PDB topology must contain at least two non-empty chains "
                "to infer umbrella pull groups"
            )
        effective_group_1 = chain_groups[0]
        effective_group_2 = chain_groups[1]
    else:
        # Fallback: midpoint heuristic with deprecation warning.
        topology = _generic_topology(system.getNumParticles())
        n = system.getNumParticles()
        effective_group_1 = list(range(n // 2))
        effective_group_2 = list(range(n // 2, n))

    production_defaults = ProductionConfig()
    window_centers = generate_window_centers(config)

    manifest = _load_manifest(output_dir)
    completed = manifest.get("completed_windows", {})

    # Separate already-completed windows from those needing computation.
    results_map: dict[int, dict[str, Any]] = {}
    pending_windows: list[tuple[int, float]] = []
    for window_id, window_center_nm in enumerate(window_centers, start=1):
        window_key = str(window_id)
        if window_key in completed:
            xi_path = Path(completed[window_key]["xi_timeseries_path"])
            if xi_path.exists():
                logger.info(
                    "Skipping completed umbrella window %d (center %.3f nm)",
                    window_id, float(window_center_nm),
                )
                xi_data = np.load(xi_path)
                results_map[window_id] = {
                    "window_id": window_id,
                    "window_center_nm": float(window_center_nm),
                    "xi_timeseries": xi_data,
                    "trajectory_path": Path(completed[window_key]["trajectory_path"]),
                    "mean_xi_nm": float(np.mean(xi_data)),
                    "std_xi_nm": float(np.std(xi_data)),
                    "xi_timeseries_path": xi_path,
                    "n_discarded_frames": 0,
                }
                continue
        pending_windows.append((window_id, float(window_center_nm)))

    if n_workers == 1 or len(pending_windows) <= 1:
        # Sequential path — zero overhead, backward compatible.
        for window_id, window_center_nm in pending_windows:
            window_system = XmlSerializer.deserialize(XmlSerializer.serialize(system))
            integrator = openmm.LangevinMiddleIntegrator(
                production_defaults.temperature_k * unit.kelvin,
                production_defaults.friction_per_ps / unit.picosecond,
                production_defaults.timestep_ps * unit.picoseconds,
            )
            simulation = Simulation(topology, window_system, integrator, select_platform(platform_name))
            logger.info("Simulation platform: %s", simulation.context.getPlatform().getName())
            simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
            simulation.context.setPositions(state.getPositions())
            simulation.context.setVelocities(state.getVelocities())
            logger.info(
                "Campaign: pre-positioning window %d/%d to xi=%.3f nm",
                window_id, len(window_centers), float(window_center_nm),
            )

            window_result = _run_umbrella_window_with_groups(
                simulation,
                float(window_center_nm),
                config,
                window_id,
                Path(output_dir),
                effective_group_1,
                effective_group_2,
                auxiliary_cvs=auxiliary_cvs,
            )
            results_map[window_id] = window_result

            # Record completion in manifest.
            completed[str(window_id)] = {
                "window_center_nm": float(window_center_nm),
                "trajectory_path": str(window_result["trajectory_path"]),
                "xi_timeseries_path": str(window_result["xi_timeseries_path"]),
            }
            manifest["completed_windows"] = completed
            _save_manifest(output_dir, manifest)
    else:
        # Parallel path — dispatch windows to a process pool.
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(
                    _run_umbrella_worker,
                    system_xml=system_xml,
                    state_xml=state_xml,
                    config=config,
                    window_id=window_id,
                    window_center_nm=window_center_nm,
                    output_dir=output_dir,
                    pull_group_1=effective_group_1,
                    pull_group_2=effective_group_2,
                    temperature_k=production_defaults.temperature_k,
                    friction_per_ps=production_defaults.friction_per_ps,
                    timestep_ps=production_defaults.timestep_ps,
                    platform_name=platform_name,
                    topology_pdb_path=effective_topo_path,
                    auxiliary_cvs=auxiliary_cvs,
                ): window_id
                for window_id, window_center_nm in pending_windows
            }
            for future in as_completed(futures):
                wid = futures[future]
                window_result = future.result()
                results_map[wid] = window_result
                # Update manifest after each completed window.
                completed[str(wid)] = {
                    "window_center_nm": window_result["window_center_nm"],
                    "trajectory_path": str(window_result["trajectory_path"]),
                    "xi_timeseries_path": str(window_result["xi_timeseries_path"]),
                }
        manifest["completed_windows"] = completed
        _save_manifest(output_dir, manifest)

    # Sort results by window_id for correct adjacency in overlap checks.
    results = [results_map[wid] for wid in sorted(results_map)]

    # --- Enhanced overlap / coverage diagnostics ---
    xi_timeseries_list = [r["xi_timeseries"] for r in results]
    window_centers_array = np.array([r["window_center_nm"] for r in results])

    # 1. Global coverage-hole detection.
    total_samples = sum(ts.size for ts in xi_timeseries_list)
    coverage_bins = min(100, max(10, total_samples // 5))
    coverage = diagnose_histogram_coverage(xi_timeseries_list, window_centers_array, n_bins=coverage_bins)
    if coverage["has_coverage_holes"]:
        # Only flag interior holes (between first and last non-zero bins).
        envelope = coverage["envelope_counts"]
        nonzero_indices = np.flatnonzero(envelope > 0)
        if nonzero_indices.size >= 2:
            interior_mask = np.zeros(len(envelope), dtype=bool)
            interior_mask[nonzero_indices[0]:nonzero_indices[-1] + 1] = True
            interior_holes = coverage["bin_centers"][(envelope == 0) & interior_mask]
            interior_span = nonzero_indices[-1] - nonzero_indices[0] + 1
            interior_hole_fraction = len(interior_holes) / interior_span if interior_span > 0 else 0.0
            if interior_hole_fraction > 0.10:
                raise PhysicalValidityError(
                    f"IV-8 violated: coverage holes detected at xi = "
                    f"{', '.join(f'{h:.3f}' for h in interior_holes[:5])} nm "
                    f"({len(interior_holes)} zero-count bins total). "
                    f"Add umbrella windows in these regions."
                )

    # 2. Adjacent overlap check (retained from original IV-8 enforcement).
    for left_result, right_result in zip(results[:-1], results[1:], strict=False):
        overlap_fraction = _histogram_overlap_fraction(left_result["xi_timeseries"], right_result["xi_timeseries"])
        if overlap_fraction < 0.10:
            raise PhysicalValidityError(
                "IV-8 violated: adjacent umbrella histograms must overlap by at least 10%"
            )

    # 3. Effective sample size warning.
    n_eff = compute_effective_sample_sizes(xi_timeseries_list)
    for i, (eff, result) in enumerate(zip(n_eff, results)):
        if eff < 50:
            logger.warning(
                "Window %d (xi=%.3f nm): N_eff=%.0f < 50 — consider extending sampling",
                result["window_id"], result["window_center_nm"], eff,
            )

    # 4. Log campaign-level diagnostics.
    overlap_matrix = compute_overlap_matrix(xi_timeseries_list)
    min_adjacent_overlap = min(
        overlap_matrix[i, i + 1] for i in range(len(results) - 1)
    ) if len(results) > 1 else 1.0
    logger.info(
        "Umbrella campaign: coverage=%.1f%%, min adjacent overlap=%.3f, "
        "min N_eff=%.0f",
        coverage["coverage_fraction"] * 100, min_adjacent_overlap, float(np.min(n_eff)),
    )

    # 5. Validate chain integrity when authentic topology was provided.
    if effective_topo_path is not None:
        _validate_trajectory_chain_integrity(topology)

    return results