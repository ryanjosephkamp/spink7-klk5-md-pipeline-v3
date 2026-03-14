"""Umbrella-sampling utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import openmm
from openmm import XmlSerializer, unit
from openmm.app import DCDFile, Element, Simulation, Topology

from src import PhysicalValidityError
from src.config import ProductionConfig, UmbrellaConfig
from src.physics.collective_variables import com_distance, com_vector


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


def _generic_topology(n_particles: int) -> Topology:
    """Construct a generic topology when only serialized state and system are available."""

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


def _run_umbrella_window_with_groups(
    simulation: Simulation,
    window_center_nm: float,
    config: UmbrellaConfig,
    window_id: int,
    output_dir: Path,
    pull_group_1: list[int],
    pull_group_2: list[int],
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

    _create_umbrella_force(simulation, group_a, group_b, window_center_nm, config.spring_constant_kj_mol_nm2)
    output_paths = _window_output_paths(Path(output_dir), window_id, window_center_nm)
    timestep_ps = simulation.integrator.getStepSize().value_in_unit(unit.picoseconds)
    report_interval_steps = max(1, round(config.save_interval_ps / timestep_ps))
    n_samples = max(1, round((config.per_window_duration_ns * 1000.0) / config.save_interval_ps))
    masses = _particle_masses_amu(simulation.system)
    xi_timeseries = np.zeros(n_samples, dtype=float)

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
            xi_timeseries[sample_index] = com_distance(
                positions_nm,
                masses,
                np.asarray(group_a, dtype=int),
                np.asarray(group_b, dtype=int),
            )
            dcd_file.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())

    np.save(output_paths["xi_timeseries_path"], xi_timeseries)
    return {
        "window_id": int(window_id),
        "window_center_nm": float(window_center_nm),
        "xi_timeseries": xi_timeseries,
        "trajectory_path": output_paths["trajectory_path"],
        "mean_xi_nm": float(np.mean(xi_timeseries)),
        "std_xi_nm": float(np.std(xi_timeseries)),
        "xi_timeseries_path": output_paths["xi_timeseries_path"],
    }


def run_umbrella_window(
    simulation: Simulation,
    window_center_nm: float,
    config: UmbrellaConfig,
    window_id: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Run a single umbrella-sampling window using groups inferred from topology chains.

    Invariants: None directly. Campaign-level overlap is checked separately under IV-8.

    Args:
        simulation: OpenMM simulation object whose first two chains define the COM groups.
        window_center_nm: Target COM distance in nm.
        config: Umbrella sampling parameters.
        window_id: Integer window identifier.
        output_dir: Directory for on-disk outputs.

    Returns:
        dict[str, Any]: Window metadata and xi timeseries. `xi_timeseries` shape: [N_samples].
    """

    group_a, group_b = _infer_pull_groups_from_topology(simulation)
    return _run_umbrella_window_with_groups(simulation, window_center_nm, config, window_id, output_dir, group_a, group_b)


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


def run_umbrella_campaign(
    equilibrated_state_path: Path,
    system_xml_path: Path,
    config: UmbrellaConfig,
    pull_group_1: list[int],
    pull_group_2: list[int],
    output_dir: Path,
) -> list[dict[str, Any]]:
    """Run an independent restrained simulation at each umbrella window center.

    Invariants: Enforces IV-8 by requiring adjacent xi histograms to overlap by at least 10%.

    Args:
        equilibrated_state_path: Serialized OpenMM state path.
        system_xml_path: Serialized OpenMM system path.
        config: Umbrella sampling parameters.
        pull_group_1: First COM group. Shape: [N_a].
        pull_group_2: Second COM group. Shape: [N_b].
        output_dir: Directory for per-window outputs.

    Returns:
        list[dict[str, Any]]: Per-window umbrella results.
    """

    _validate_config(config)
    system = XmlSerializer.deserialize(Path(system_xml_path).read_text(encoding="utf-8"))
    state = XmlSerializer.deserialize(Path(equilibrated_state_path).read_text(encoding="utf-8"))
    topology = _generic_topology(system.getNumParticles())
    production_defaults = ProductionConfig()
    window_centers = generate_window_centers(config)

    results: list[dict[str, Any]] = []
    for window_id, window_center_nm in enumerate(window_centers, start=1):
        window_system = XmlSerializer.deserialize(XmlSerializer.serialize(system))
        integrator = openmm.LangevinMiddleIntegrator(
            production_defaults.temperature_k * unit.kelvin,
            production_defaults.friction_per_ps / unit.picosecond,
            production_defaults.timestep_ps * unit.picoseconds,
        )
        simulation = Simulation(topology, window_system, integrator, openmm.Platform.getPlatformByName("CPU"))
        simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        simulation.context.setPositions(state.getPositions())
        simulation.context.setVelocities(state.getVelocities())
        results.append(
            _run_umbrella_window_with_groups(
                simulation,
                float(window_center_nm),
                config,
                window_id,
                Path(output_dir),
                pull_group_1,
                pull_group_2,
            )
        )

    for left_result, right_result in zip(results[:-1], results[1:], strict=False):
        overlap_fraction = _histogram_overlap_fraction(left_result["xi_timeseries"], right_result["xi_timeseries"])
        if overlap_fraction < 0.10:
            raise PhysicalValidityError(
                "IV-8 violated: adjacent umbrella histograms must overlap by at least 10%"
            )

    return results