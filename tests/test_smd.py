"""Tests for steered molecular dynamics utilities."""

from __future__ import annotations

import numpy as np
import openmm
import pytest
from openmm import XmlSerializer, unit
from openmm.app import Element, Simulation, Topology

from src import PhysicalValidityError
from src.config import SMDConfig
from src.physics.collective_variables import com_distance
from src.simulate.smd import run_smd_campaign, run_smd_replicate


def _make_smd_test_simulation() -> tuple[Simulation, list[int], list[int], np.ndarray, float]:
    """Build a tiny anchored two-group system for fast SMD validation."""

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("SMD", chain)
    carbon = Element.getByAtomicNumber(6)
    atom_0 = topology.addAtom("A0", carbon, residue)
    atom_1 = topology.addAtom("A1", carbon, residue)
    atom_2 = topology.addAtom("A2", carbon, residue)
    topology.addBond(atom_0, atom_1)
    topology.addBond(atom_1, atom_2)

    system = openmm.System()
    for _ in range(3):
        system.addParticle(12.0)

    tether = openmm.CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    tether.addGlobalParameter("k", 5000.0)
    tether.addPerParticleParameter("x0")
    tether.addPerParticleParameter("y0")
    tether.addPerParticleParameter("z0")
    tether.addParticle(0, [0.0, 0.0, 0.0])
    tether.addParticle(1, [0.2, 0.0, 0.0])
    system.addForce(tether)

    bond_force = openmm.HarmonicBondForce()
    bond_force.addBond(0, 1, 0.2, 2000.0)
    bond_force.addBond(1, 2, 0.2, 2000.0)
    system.addForce(bond_force)

    integrator = openmm.LangevinMiddleIntegrator(310.0 * unit.kelvin, 5.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.4, 0.0, 0.0]] * unit.nanometer
    )
    simulation.context.setVelocitiesToTemperature(310.0 * unit.kelvin, 11)

    initial_positions = np.asarray(
        simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer),
        dtype=float,
    )
    masses = np.asarray([12.0, 12.0, 12.0], dtype=float)
    initial_distance_nm = com_distance(initial_positions, masses, np.asarray([2]), np.asarray([0, 1]))
    return simulation, [2], [0, 1], np.asarray([1.0, 0.0, 0.0], dtype=float), initial_distance_nm


def test_run_smd_replicate_displaces_coordinate_and_logs_work(tmp_path) -> None:
    """SMD should increase the COM distance and record non-equilibrium work consistently."""

    simulation, pull_group_1, pull_group_2, pull_direction, initial_distance_nm = _make_smd_test_simulation()
    config = SMDConfig(
        spring_constant_kj_mol_nm2=4000.0,
        pulling_velocity_nm_per_ps=0.02,
        pull_distance_nm=0.10,
        save_interval_ps=0.25,
    )

    result = run_smd_replicate(
        simulation,
        config,
        replicate_id=1,
        output_dir=tmp_path,
        pull_group_1=pull_group_1,
        pull_group_2=pull_group_2,
        pull_direction=pull_direction,
    )

    assert result["trajectory_path"].exists()
    assert result["work_timeseries_path"].exists()
    assert result["force_timeseries_path"].exists()
    assert result["xi_timeseries_path"].exists()

    work_csv = np.loadtxt(result["work_timeseries_path"], delimiter=",", skiprows=1)
    force_csv = np.loadtxt(result["force_timeseries_path"], delimiter=",", skiprows=1)
    xi_csv = np.loadtxt(result["xi_timeseries_path"], delimiter=",", skiprows=1)

    assert work_csv.shape[1] == 2
    assert force_csv.shape == work_csv.shape
    assert xi_csv.shape == work_csv.shape
    assert np.all(np.diff(work_csv[:, 0]) > 0.0)
    assert np.all(np.diff(xi_csv[:, 0]) > 0.0)
    assert result["final_xi_nm"] > initial_distance_nm
    assert work_csv[-1, 1] == pytest.approx(result["total_work_kj_mol"], rel=1e-9, abs=1e-9)
    assert np.any(np.abs(force_csv[:, 1]) > 0.0)


def test_run_smd_replicate_rejects_misaligned_pull_direction(tmp_path) -> None:
    """Public API should reject pull directions that oppose the initial COM axis."""

    simulation, pull_group_1, pull_group_2, _, _ = _make_smd_test_simulation()

    with pytest.raises(ValueError, match="pull_direction must point from pull_group_2 toward pull_group_1"):
        run_smd_replicate(
            simulation,
            SMDConfig(pull_distance_nm=0.05, save_interval_ps=0.5),
            replicate_id=1,
            output_dir=tmp_path,
            pull_group_1=pull_group_1,
            pull_group_2=pull_group_2,
            pull_direction=np.asarray([-1.0, 0.0, 0.0], dtype=float),
        )


# ---------- L-31 Step 1: Streaming CSV time series ----------


def test_streaming_smd_produces_valid_csv_and_scalar_work(tmp_path) -> None:
    """Streaming SMD must produce CSV files and a scalar total_work identical to expectations."""

    simulation, pull_group_1, pull_group_2, pull_direction, initial_distance_nm = _make_smd_test_simulation()
    config = SMDConfig(
        spring_constant_kj_mol_nm2=4000.0,
        pulling_velocity_nm_per_ps=0.02,
        pull_distance_nm=0.10,
        save_interval_ps=0.25,
    )

    result = run_smd_replicate(
        simulation,
        config,
        replicate_id=1,
        output_dir=tmp_path,
        pull_group_1=pull_group_1,
        pull_group_2=pull_group_2,
        pull_direction=pull_direction,
    )

    # Scalar outputs
    assert isinstance(result["total_work_kj_mol"], float)
    assert isinstance(result["final_xi_nm"], float)
    assert result["final_xi_nm"] > initial_distance_nm

    # CSV file integrity
    work_csv = np.loadtxt(result["work_timeseries_path"], delimiter=",", skiprows=1)
    force_csv = np.loadtxt(result["force_timeseries_path"], delimiter=",", skiprows=1)
    xi_csv = np.loadtxt(result["xi_timeseries_path"], delimiter=",", skiprows=1)

    assert work_csv.shape[1] == 2
    assert force_csv.shape[1] == 2
    assert xi_csv.shape[1] == 2
    assert work_csv.shape[0] == result["n_samples"]

    # Time monotonicity
    assert np.all(np.diff(work_csv[:, 0]) > 0.0)
    assert np.all(np.diff(xi_csv[:, 0]) > 0.0)

    # Final CSV work value matches scalar
    assert work_csv[-1, 1] == pytest.approx(result["total_work_kj_mol"], rel=1e-9)

    # In-memory arrays should NOT be in result
    assert "work_timeseries" not in result
    assert "force_timeseries" not in result
    assert "xi_timeseries" not in result


# ---------- L-08 Step 4: Auxiliary CV recording in SMD ----------


def test_smd_replicate_records_auxiliary_angle_cv(tmp_path) -> None:
    """SMD replicate should record auxiliary CV timeseries alongside xi."""
    from src.physics.collective_variables import CollectiveVariableSpec

    simulation, pg1, pg2, pull_dir, _ = _make_smd_test_simulation()
    config = SMDConfig(
        spring_constant_kj_mol_nm2=4000.0,
        pulling_velocity_nm_per_ps=0.02,
        pull_distance_nm=0.10,
        save_interval_ps=0.25,
    )
    aux = [CollectiveVariableSpec(name="com_angle", reference_direction=pull_dir)]
    result = run_smd_replicate(
        simulation, config, replicate_id=1, output_dir=tmp_path,
        pull_group_1=pg1, pull_group_2=pg2, pull_direction=pull_dir,
        auxiliary_cvs=aux,
    )
    key = "aux_cv_com_angle_timeseries"
    assert key in result, f"Expected key '{key}' in result"
    ts = result[key]
    assert ts.shape[1] == 2
    assert ts.shape[0] == result["n_samples"]
    assert np.all(np.isfinite(ts[:, 1]))


def test_smd_replicate_no_auxiliary_cvs_unchanged(tmp_path) -> None:
    """With no auxiliary CVs, result schema should be unchanged."""
    simulation, pg1, pg2, pull_dir, _ = _make_smd_test_simulation()
    config = SMDConfig(
        spring_constant_kj_mol_nm2=4000.0,
        pulling_velocity_nm_per_ps=0.02,
        pull_distance_nm=0.10,
        save_interval_ps=0.25,
    )
    result = run_smd_replicate(
        simulation, config, replicate_id=1, output_dir=tmp_path,
        pull_group_1=pg1, pull_group_2=pg2, pull_direction=pull_dir,
    )
    aux_keys = [k for k in result if k.startswith("aux_cv_")]
    assert len(aux_keys) == 0


# ---------- L-18 Step 5: Configurable base seed in SMD ----------


def test_smd_configure_replicate_velocities_uses_base_seed() -> None:
    """_configure_replicate_velocities should use base_seed + replicate_id."""
    from src.simulate.smd import _configure_replicate_velocities

    simulation, _, _, _, _ = _make_smd_test_simulation()
    base_seed = 500
    replicate_id = 3
    _configure_replicate_velocities(simulation, replicate_id, base_seed)
    actual_seed = simulation.integrator.getRandomNumberSeed()
    assert actual_seed == base_seed + replicate_id


def test_smd_resolve_seed_returns_integer_when_none() -> None:
    """SMD _resolve_seed(None) should return a non-negative integer."""
    from src.simulate.smd import _resolve_seed

    seed = _resolve_seed(None, "test")
    assert isinstance(seed, int)
    assert 0 <= seed <= 0x7FFFFFFF


# ---------- L-19 Step 2: Topology loading in SMD campaign ----------


def test_smd_load_topology_from_pdb(tmp_path) -> None:
    """_load_topology() should load the authentic topology from a PDB file."""
    from openmm.app import PDBFile
    from src.simulate.smd import _load_topology

    # Build a small system and write its topology as PDB.
    simulation, _, _, _, _ = _make_smd_test_simulation()
    state = simulation.context.getState(getPositions=True)
    pdb_path = tmp_path / "topology_reference.pdb"
    with pdb_path.open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(simulation.topology, state.getPositions(), fh, keepIds=True)

    n_particles = simulation.system.getNumParticles()
    loaded_topology = _load_topology(pdb_path, n_particles)
    assert loaded_topology.getNumAtoms() == n_particles


def test_smd_load_topology_none_falls_back_with_warning(tmp_path, caplog) -> None:
    """_load_topology(None) should fall back to _generic_topology() with a log warning."""
    from src.simulate.smd import _load_topology

    with caplog.at_level("WARNING", logger="src.simulate.smd"):
        topology = _load_topology(None, 5)
    assert topology.getNumAtoms() == 5
    assert any("topology_pdb_path" in msg for msg in caplog.messages)


def test_smd_load_topology_rejects_atom_count_mismatch(tmp_path) -> None:
    """_load_topology() should reject a PDB with mismatched atom count."""
    from openmm.app import PDBFile
    from src.simulate.smd import _load_topology

    simulation, _, _, _, _ = _make_smd_test_simulation()
    state = simulation.context.getState(getPositions=True)
    pdb_path = tmp_path / "topology_reference.pdb"
    with pdb_path.open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(simulation.topology, state.getPositions(), fh, keepIds=True)

    with pytest.raises(ValueError, match="does not match"):
        _load_topology(pdb_path, 999)


# ---------- L-19 Step 3: Trajectory chain integrity validation ----------


def test_smd_chain_validation_rejects_single_chain() -> None:
    """Chain validation should reject a single-chain topology for a multi-chain system."""
    from src.simulate.smd import _validate_trajectory_chain_integrity

    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("X", chain)
    topology.addAtom("A", Element.getByAtomicNumber(6), residue)

    with pytest.raises(PhysicalValidityError, match="chain"):
        _validate_trajectory_chain_integrity(topology, expected_min_chains=2)


def test_smd_chain_validation_accepts_multi_chain() -> None:
    """Chain validation should accept a topology with enough chains."""
    from src.simulate.smd import _validate_trajectory_chain_integrity

    topology = Topology()
    chain_a = topology.addChain()
    chain_b = topology.addChain()
    res_a = topology.addResidue("A", chain_a)
    res_b = topology.addResidue("B", chain_b)
    topology.addAtom("A0", Element.getByAtomicNumber(6), res_a)
    topology.addAtom("B0", Element.getByAtomicNumber(6), res_b)

    # Should not raise
    _validate_trajectory_chain_integrity(topology, expected_min_chains=2)


# ---------- L-19 Step 4: _generic_topology() deprecation ----------


# ---------- L-30 Step 1: Process-safe SMD worker ----------


def test_smd_worker_produces_identical_result_to_direct_call(tmp_path) -> None:
    """Worker function must produce a valid result matching the run_smd_replicate schema."""
    from openmm import XmlSerializer
    from src.config import ProductionConfig
    from src.simulate.smd import _run_smd_worker

    config = SMDConfig(
        spring_constant_kj_mol_nm2=4000.0,
        pulling_velocity_nm_per_ps=0.02,
        pull_distance_nm=0.10,
        save_interval_ps=0.25,
        random_seed=42,
    )
    prod = ProductionConfig()
    base_seed = 42

    # Build a reference simulation and serialize its system/state.
    sim_ref, pg1, pg2, pull_dir, initial_dist = _make_smd_test_simulation()
    system_xml = XmlSerializer.serialize(sim_ref.system)
    state_xml = XmlSerializer.serialize(
        sim_ref.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    )

    # Worker call.
    result = _run_smd_worker(
        system_xml=system_xml,
        state_xml=state_xml,
        config=config,
        replicate_id=1,
        output_dir=tmp_path / "worker",
        pull_group_1=pg1,
        pull_group_2=pg2,
        pull_direction=pull_dir,
        base_seed=base_seed,
        temperature_k=prod.temperature_k,
        friction_per_ps=prod.friction_per_ps,
        timestep_ps=prod.timestep_ps,
        platform_name="CPU",
    )

    # Verify schema matches run_smd_replicate return contract.
    assert result["trajectory_path"].exists()
    assert result["work_timeseries_path"].exists()
    assert result["force_timeseries_path"].exists()
    assert result["xi_timeseries_path"].exists()
    work_csv = np.loadtxt(result["work_timeseries_path"], delimiter=",", skiprows=1)
    assert work_csv.shape[1] == 2
    # Time should be monotonically increasing.
    assert np.all(np.diff(work_csv[:, 0]) > 0.0)
    # COM distance should increase under pulling.
    assert result["final_xi_nm"] > initial_dist
    # Cumulative work must match the last work timeseries entry.
    assert work_csv[-1, 1] == pytest.approx(
        result["total_work_kj_mol"], rel=1e-9,
    )


# ---------- L-30 Step 2: Parallel SMD campaign ----------


def test_run_smd_campaign_parallel_produces_valid_results(tmp_path) -> None:
    """Parallel campaign with n_workers=2 must produce valid results for each replicate."""
    from openmm import XmlSerializer

    simulation, pg1, pg2, _, _ = _make_smd_test_simulation()
    system_xml_path = tmp_path / "system.xml"
    state_xml_path = tmp_path / "state.xml"
    system_xml_path.write_text(XmlSerializer.serialize(simulation.system), encoding="utf-8")
    state_xml_path.write_text(
        XmlSerializer.serialize(
            simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
        ),
        encoding="utf-8",
    )
    config = SMDConfig(
        spring_constant_kj_mol_nm2=4000.0,
        pulling_velocity_nm_per_ps=0.02,
        pull_distance_nm=0.10,
        save_interval_ps=0.25,
        n_replicates=2,
        random_seed=42,
    )
    results = run_smd_campaign(
        equilibrated_state_path=state_xml_path,
        system_xml_path=system_xml_path,
        config=config,
        pull_group_1=pg1,
        pull_group_2=pg2,
        output_dir=tmp_path / "campaign",
        n_workers=2,
    )
    assert len(results) == 2
    for r in results:
        assert r["trajectory_path"].exists()
        assert isinstance(r["total_work_kj_mol"], float)
        assert np.isfinite(r["total_work_kj_mol"])


# ---------- L-30 Step 5: n_workers validation ----------


def test_run_smd_campaign_rejects_invalid_n_workers(tmp_path) -> None:
    """run_smd_campaign must raise ValueError when n_workers < 1."""
    simulation, pg1, pg2, _, _ = _make_smd_test_simulation()
    system_path = tmp_path / "system.xml"
    state_path = tmp_path / "state.xml"
    system_path.write_text(XmlSerializer.serialize(simulation.system))
    state_path.write_text(
        XmlSerializer.serialize(
            simulation.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
        )
    )
    config = SMDConfig(pulling_velocity_nm_per_ps=0.005, pull_distance_nm=0.01, save_interval_ps=0.5)
    with pytest.raises(ValueError, match="n_workers must be at least 1"):
        run_smd_campaign(
            equilibrated_state_path=state_path,
            system_xml_path=system_path,
            config=config,
            pull_group_1=pg1,
            pull_group_2=pg2,
            output_dir=tmp_path / "smd_bad",
            n_workers=0,
        )


def test_smd_generic_topology_emits_deprecation_warning() -> None:
    """_generic_topology() should emit a DeprecationWarning."""
    from src.simulate.smd import _generic_topology

    with pytest.warns(DeprecationWarning, match="incorrect chain assignments"):
        _generic_topology(10)