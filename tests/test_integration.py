"""End-to-end integration test for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import openmm
import pytest
from openmm import XmlSerializer, unit
from openmm.app import ForceField, HBonds, PME, PDBFile, Simulation

from src.analyze.structural import compute_rmsd
from src.analyze.trajectory import load_trajectory
from src.analyze.wham import solve_wham
from src.config import EquilibrationConfig, MinimizationConfig, SMDConfig, SystemConfig, UmbrellaConfig, WHAMConfig
from src.physics.collective_variables import com_distance, com_vector
from src.physics.restraints import create_positional_restraints
from src.prep.solvate import solvate_system
from src.prep.topology import build_topology
from src.simulate.equilibrate import run_npt, run_nvt
from src.simulate.minimizer import minimize_energy
from src.simulate.smd import run_smd_campaign
from src.simulate.umbrella import generate_window_centers, run_umbrella_campaign


_SOLVENT_RESIDUES = {"HOH", "WAT", "NA", "CL"}


def _build_restrained_equilibration_simulation(
    alanine_dipeptide_pdb: Path,
) -> tuple[Simulation, object, ForceField, int]:
    """Build a solvated alanine-dipeptide simulation with equilibration restraints."""

    system_config = SystemConfig(box_padding_nm=1.2, ionic_strength_molar=0.0)
    _, _, modeller = build_topology(alanine_dipeptide_pdb, system_config)
    solute_atom_count = modeller.topology.getNumAtoms()
    solvated_modeller, _, _, _ = solvate_system(modeller, system_config)

    force_field = ForceField(system_config.force_field, system_config.water_model)
    system = force_field.createSystem(
        solvated_modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )

    restrained_indices: list[int] = []
    reference_positions: list[list[float]] = []
    for atom_index, atom in enumerate(solvated_modeller.topology.atoms()):
        if atom.residue.name.upper() in _SOLVENT_RESIDUES:
            continue
        if atom.element is None or atom.element.symbol == "H":
            continue
        restrained_indices.append(atom_index)
        position_nm = solvated_modeller.positions[atom_index].value_in_unit(unit.nanometer)
        reference_positions.append([float(position_nm[0]), float(position_nm[1]), float(position_nm[2])])

    create_positional_restraints(
        system,
        restrained_indices,
        reference_positions=np.asarray(reference_positions, dtype=float),
        force_constant_kj_mol_nm2=500.0,
    )

    integrator = openmm.LangevinMiddleIntegrator(310.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(solvated_modeller.topology, system, integrator, platform)
    simulation.context.setPositions(solvated_modeller.positions)
    return simulation, solvated_modeller, force_field, solute_atom_count


def _terminal_residue_groups(topology, solute_atom_count: int) -> tuple[list[int], list[int]]:
    """Split the solvated peptide into terminal COM groups using the first and last solute residues."""

    solute_residues = [residue for residue in topology.residues() if residue.name.upper() not in _SOLVENT_RESIDUES]
    if len(solute_residues) < 2:
        raise ValueError("integration test requires at least two solute residues")

    group_a = [atom.index for atom in solute_residues[-1].atoms() if atom.index < solute_atom_count]
    group_b = [atom.index for atom in solute_residues[0].atoms() if atom.index < solute_atom_count]
    return group_a, group_b


def _write_topology_pdb(topology, positions, output_path: Path) -> Path:
    """Write a solvated topology PDB for MDTraj loading."""

    with output_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology, positions, handle, keepIds=True)
    return output_path


def _serialize_sampling_inputs(
    force_field: ForceField,
    topology,
    state_xml_path: Path,
    output_dir: Path,
    system_config: SystemConfig,
) -> Path:
    """Serialize an unrestrained sampling system aligned with the equilibrated state."""

    sampling_system = force_field.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )
    system_xml_path = output_dir / "sampling_system.xml"
    system_xml_path.write_text(XmlSerializer.serialize(sampling_system), encoding="utf-8")
    return system_xml_path


def _reaction_coordinate_state(
    final_state_path: Path,
    system_xml_path: Path,
    pull_group_1: list[int],
    pull_group_2: list[int],
) -> tuple[np.ndarray, float]:
    """Extract a deterministic pull direction and initial COM distance from serialized inputs."""

    state = XmlSerializer.deserialize(final_state_path.read_text(encoding="utf-8"))
    system = XmlSerializer.deserialize(system_xml_path.read_text(encoding="utf-8"))
    positions_nm = np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float)
    masses = np.asarray(
        [system.getParticleMass(index).value_in_unit(unit.dalton) for index in range(system.getNumParticles())],
        dtype=float,
    )
    com_a = com_vector(positions_nm, masses, np.asarray(pull_group_1, dtype=int))
    com_b = com_vector(positions_nm, masses, np.asarray(pull_group_2, dtype=int))
    pull_direction = com_a - com_b
    pull_direction /= np.linalg.norm(pull_direction)
    initial_distance_nm = com_distance(
        positions_nm,
        masses,
        np.asarray(pull_group_1, dtype=int),
        np.asarray(pull_group_2, dtype=int),
    )
    return pull_direction, float(initial_distance_nm)


def test_end_to_end_alanine_dipeptide_pipeline(tmp_output_dir: Path, alanine_dipeptide_pdb: Path) -> None:
    """Run the Chunk 24 end-to-end alanine-dipeptide pipeline in under five minutes on CPU."""

    start_time = time.perf_counter()
    system_config = SystemConfig(box_padding_nm=1.2, ionic_strength_molar=0.0)
    simulation, solvated_modeller, force_field, solute_atom_count = _build_restrained_equilibration_simulation(
        alanine_dipeptide_pdb
    )
    pull_group_1, pull_group_2 = _terminal_residue_groups(solvated_modeller.topology, solute_atom_count)

    minimization_result = minimize_energy(
        simulation,
        MinimizationConfig(max_iterations=100, tolerance_kj_mol_nm=10.0),
    )
    assert minimization_result["final_energy_kj_mol"] < minimization_result["initial_energy_kj_mol"]

    nvt_result = run_nvt(
        simulation,
        EquilibrationConfig(nvt_duration_ps=100.0, friction_per_ps=10.0, save_interval_ps=0.5),
        tmp_output_dir / "nvt",
    )
    assert abs(nvt_result["avg_temperature_k"] - 310.0) < 5.0

    npt_result = run_npt(
        simulation,
        EquilibrationConfig(npt_duration_ps=100.0, friction_per_ps=10.0, barostat_interval=5, save_interval_ps=0.5),
        tmp_output_dir / "npt",
    )
    assert abs(npt_result["avg_temperature_k"] - 310.0) < 5.0
    assert 0.95 < npt_result["avg_density_g_cm3"] < 1.05

    sampling_state_path = tmp_output_dir / "sampling_state.xml"
    sampling_state_path.write_text(Path(npt_result["final_state_path"]).read_text(encoding="utf-8"), encoding="utf-8")
    sampling_system_path = _serialize_sampling_inputs(
        force_field,
        solvated_modeller.topology,
        sampling_state_path,
        tmp_output_dir,
        system_config,
    )

    _, initial_distance_nm = _reaction_coordinate_state(
        sampling_state_path,
        sampling_system_path,
        pull_group_1,
        pull_group_2,
    )

    smd_results = run_smd_campaign(
        sampling_state_path,
        sampling_system_path,
        SMDConfig(
            spring_constant_kj_mol_nm2=2000.0,
            pulling_velocity_nm_per_ps=0.02,
            pull_distance_nm=0.02,
            n_replicates=2,
            save_interval_ps=0.25,
        ),
        pull_group_1,
        pull_group_2,
        tmp_output_dir / "smd",
    )
    assert len(smd_results) == 2
    assert all(result["trajectory_path"].exists() for result in smd_results)
    assert all(np.isfinite(result["total_work_kj_mol"]) for result in smd_results)

    umbrella_config = UmbrellaConfig(
        xi_min_nm=max(0.05, initial_distance_nm - 0.01),
        xi_max_nm=initial_distance_nm + 0.01,
        window_spacing_nm=0.01,
        spring_constant_kj_mol_nm2=1200.0,
        per_window_duration_ns=0.01,
        save_interval_ps=0.5,
    )
    umbrella_results = run_umbrella_campaign(
        sampling_state_path,
        sampling_system_path,
        umbrella_config,
        pull_group_1,
        pull_group_2,
        tmp_output_dir / "umbrella",
    )
    assert len(umbrella_results) == 3
    assert all(result["trajectory_path"].exists() for result in umbrella_results)

    wham_result = solve_wham(
        [result["xi_timeseries"] for result in umbrella_results],
        generate_window_centers(umbrella_config),
        np.full(3, umbrella_config.spring_constant_kj_mol_nm2, dtype=float),
        310.0,
        WHAMConfig(tolerance=1e-6, max_iterations=25_000, n_bootstrap=4, histogram_bins=24),
    )
    assert wham_result["converged"] is True
    assert wham_result["n_iterations"] > 0

    topology_path = _write_topology_pdb(
        solvated_modeller.topology,
        solvated_modeller.positions,
        tmp_output_dir / "solvated_topology.pdb",
    )
    nvt_trajectory = load_trajectory(nvt_result["trajectory_path"], topology_path)
    rmsd_nm = compute_rmsd(nvt_trajectory, nvt_trajectory[0], atom_selection="backbone")
    assert rmsd_nm.shape == (nvt_trajectory.n_frames,)
    assert rmsd_nm[0] == pytest.approx(0.0, abs=1e-8)

    elapsed_seconds = time.perf_counter() - start_time
    assert elapsed_seconds < 300.0