"""Tests for production-MD utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import openmm
import pytest
from openmm import LangevinMiddleIntegrator, unit
from openmm.app import ForceField, HBonds, PME, Simulation

from src.config import EquilibrationConfig, MinimizationConfig, ProductionConfig, SystemConfig
from src.prep.solvate import solvate_system
from src.prep.topology import build_topology
from src.simulate.equilibrate import run_npt, run_nvt
from src.simulate.minimizer import minimize_energy
from src.simulate.production import run_production


def _write_neutral_protonated_peptide_pdb(pdb_path: Path) -> Path:
    """Write a compact neutral peptide suitable for fast production tests."""

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        """ATOM      1  N   ALA A   1      -0.648   1.210   0.000  1.00 20.00           N
ATOM      2  H   ALA A   1      -1.300   1.500   0.730  1.00 20.00           H
ATOM      3  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      4  C   ALA A   1       1.526   0.000   0.000  1.00 20.00           C
ATOM      5  O   ALA A   1       2.046  -1.143   0.000  1.00 20.00           O
ATOM      6  CB  ALA A   1      -0.507  -0.774  -1.206  1.00 20.00           C
ATOM      7  N   GLY A   2       2.260   1.105   0.000  1.00 20.00           N
ATOM      8  H   GLY A   2       1.836   2.019   0.000  1.00 20.00           H
ATOM      9  CA  GLY A   2       3.714   1.084   0.000  1.00 20.00           C
ATOM     10  C   GLY A   2       4.234  -0.337   0.000  1.00 20.00           C
ATOM     11  O   GLY A   2       5.430  -0.596   0.000  1.00 20.00           O
ATOM     12  OXT GLY A   2       3.451  -1.314   0.000  1.00 20.00           O
TER
END
""",
        encoding="utf-8",
    )
    return pdb_path


def _make_equilibrated_simulation(tmp_path: Path) -> Simulation:
    """Create a minimized, briefly equilibrated solvated peptide simulation."""

    system_config = SystemConfig(box_padding_nm=1.2, ionic_strength_molar=0.0)
    protonated_path = _write_neutral_protonated_peptide_pdb(
        tmp_path / "data" / "pdb" / "prepared" / "production_peptide.pdb"
    )
    _, _, modeller = build_topology(protonated_path, system_config)
    solvated_modeller, _, _, _ = solvate_system(modeller, system_config)

    force_field = ForceField(system_config.force_field, system_config.water_model)
    system = force_field.createSystem(
        solvated_modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )

    integrator = LangevinMiddleIntegrator(310.0 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    simulation = Simulation(solvated_modeller.topology, system, integrator, platform)
    simulation.context.setPositions(solvated_modeller.positions)
    minimize_energy(simulation, MinimizationConfig(max_iterations=1000, tolerance_kj_mol_nm=10.0))
    run_nvt(
        simulation,
        EquilibrationConfig(nvt_duration_ps=100.0, friction_per_ps=10.0, save_interval_ps=0.5),
        tmp_path / "preproduction_nvt",
    )
    run_npt(
        simulation,
        EquilibrationConfig(npt_duration_ps=100.0, friction_per_ps=10.0, barostat_interval=5, save_interval_ps=1.0),
        tmp_path / "preproduction_npt",
    )
    return simulation


def test_run_production_writes_streamed_outputs_and_respects_iv5(tmp_path: Path) -> None:
    """IV-5: production MD should pass its internal NVE-reference drift diagnostic."""

    simulation = _make_equilibrated_simulation(tmp_path)
    config = ProductionConfig(duration_ns=0.02, save_interval_ps=1.0, checkpoint_interval_ps=5.0)

    result = run_production(simulation, config, tmp_path / "production_output")

    assert result["trajectory_path"].exists()
    assert result["energy_timeseries_path"].exists()
    assert result["n_frames"] == 20
    assert result["total_time_ns"] == pytest.approx(0.02, abs=1e-9)
    assert len(result["checkpoint_paths"]) >= 4
    assert all(checkpoint_path.exists() for checkpoint_path in result["checkpoint_paths"])

    energy_timeseries = np.loadtxt(result["energy_timeseries_path"], delimiter=",", skiprows=1)
    assert energy_timeseries.shape == (20, 3)
    assert np.all(np.diff(energy_timeseries[:, 0]) > 0.0)
    assert np.all(np.isfinite(energy_timeseries))
    assert simulation.currentStep > 0


def test_run_production_rejects_non_positive_duration(tmp_path: Path) -> None:
    """Public API should reject non-physical production durations."""

    simulation = _make_equilibrated_simulation(tmp_path)

    with pytest.raises(ValueError, match="config.duration_ns must be positive"):
        run_production(simulation, ProductionConfig(duration_ns=0.0), tmp_path / "invalid_production")


# ---------- L-18 Step 4: Configurable seed in production ----------


def test_run_production_with_none_seed_returns_seed(tmp_path: Path) -> None:
    """Production MD with random_seed=None should return the auto-generated seed."""

    simulation = _make_equilibrated_simulation(tmp_path)
    config = ProductionConfig(duration_ns=0.001, save_interval_ps=0.5, checkpoint_interval_ps=1.0)

    result = run_production(simulation, config, tmp_path / "production")

    assert isinstance(result["random_seed"], int)
    assert result["random_seed"] >= 0


# ---------- L-23 Step 1: Checkpoint resume ----------


def test_resume_production_continues_from_checkpoint(tmp_path: Path) -> None:
    """resume_production() should continue from the latest checkpoint."""

    from src.simulate.production import resume_production

    simulation = _make_equilibrated_simulation(tmp_path)
    config = ProductionConfig(
        duration_ns=0.001,  # 1 ps total
        save_interval_ps=0.1,
        checkpoint_interval_ps=0.5,
    )

    # Run initial production
    result_initial = run_production(simulation, config, tmp_path / "prod_output")
    assert result_initial["n_frames"] > 0

    # Resume — since the initial run already completed the full duration,
    # this should return with n_frames == 0 (nothing left to do).
    result_resumed = resume_production(simulation, config, tmp_path / "prod_output")

    assert result_resumed["trajectory_path"].exists()
    assert result_resumed["total_time_ns"] > 0


def test_resume_production_raises_when_no_checkpoint(tmp_path: Path) -> None:
    """resume_production() should raise FileNotFoundError when no checkpoint exists."""

    from unittest.mock import MagicMock

    from src.simulate.production import resume_production

    # resume_production raises FileNotFoundError before using the simulation,
    # so a mock is sufficient and avoids the expensive equilibration setup.
    simulation = MagicMock()
    config = ProductionConfig(duration_ns=0.001)
    output_dir = tmp_path / "empty_output"
    output_dir.mkdir()

    with pytest.raises(FileNotFoundError, match="No production checkpoint found"):
        resume_production(simulation, config, output_dir)


# ---------- L-23 Step 3: CLI --resume flag ----------


def test_production_cli_accepts_resume_flag() -> None:
    """Production CLI should accept --resume without error."""

    from scripts.run_production import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "--topology-pdb", "fake.pdb",
        "--system-xml", "fake.xml",
        "--state-xml", "fake.xml",
        "--resume",
    ])
    assert args.resume is True


def test_umbrella_cli_accepts_resume_flag() -> None:
    """Umbrella CLI should accept --resume without error."""

    from scripts.run_umbrella import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "--state-xml", "fake.xml",
        "--system-xml", "fake.xml",
        "--pull-group-1", "0,1",
        "--pull-group-2", "2,3",
        "--resume",
    ])
    assert args.resume is True


# ---------- L-25 Step 1: Platform selection utility ----------


def test_select_platform_returns_valid_platform() -> None:
    """select_platform() must return a valid OpenMM platform."""

    import openmm as _openmm
    from src.simulate.platform import select_platform

    platform = select_platform()
    assert isinstance(platform, _openmm.Platform)
    assert platform.getName() in ("CUDA", "OpenCL", "CPU")


def test_select_platform_explicit_cpu() -> None:
    """select_platform('CPU') must always return the CPU platform."""

    from src.simulate.platform import select_platform

    platform = select_platform("CPU")
    assert platform.getName() == "CPU"


def test_select_platform_invalid_raises() -> None:
    """select_platform('NonExistent') must raise RuntimeError."""

    from src.simulate.platform import select_platform

    with pytest.raises(RuntimeError, match="not available"):
        select_platform("NonExistent")


def test_platform_performance_summary_returns_dict() -> None:
    """platform_performance_summary must return a dict for any platform."""

    from src.simulate.platform import platform_performance_summary, select_platform

    platform = select_platform("CPU")
    summary = platform_performance_summary(platform)
    assert isinstance(summary, dict)


# ---------- L-25 Step 2: Source module platform_name parameter ----------


def test_smd_campaign_accepts_platform_name() -> None:
    """run_smd_campaign must accept a platform_name argument."""

    import inspect
    from src.simulate.smd import run_smd_campaign

    sig = inspect.signature(run_smd_campaign)
    assert "platform_name" in sig.parameters


def test_umbrella_campaign_accepts_platform_name() -> None:
    """run_umbrella_campaign must accept a platform_name argument."""

    import inspect
    from src.simulate.umbrella import run_umbrella_campaign

    sig = inspect.signature(run_umbrella_campaign)
    assert "platform_name" in sig.parameters


# ---------- L-25 Step 3: CLI --platform flag ----------


def test_equilibration_help_includes_platform_flag() -> None:
    """run_equilibration.py --help must show the --platform option."""

    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-m", "scripts.run_equilibration", "--help"],
        capture_output=True, text=True, timeout=30,
    )
    assert "--platform" in result.stdout


def test_production_help_includes_platform_flag() -> None:
    """run_production.py --help must show the --platform option."""

    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-m", "scripts.run_production", "--help"],
        capture_output=True, text=True, timeout=30,
    )
    assert "--platform" in result.stdout


def test_smd_help_includes_platform_flag() -> None:
    """run_smd.py --help must show the --platform option."""

    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-m", "scripts.run_smd", "--help"],
        capture_output=True, text=True, timeout=30,
    )
    assert "--platform" in result.stdout


def test_umbrella_help_includes_platform_flag() -> None:
    """run_umbrella.py --help must show the --platform option."""

    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-m", "scripts.run_umbrella", "--help"],
        capture_output=True, text=True, timeout=30,
    )
    assert "--platform" in result.stdout


# ---------- L-25 Step 4: Integration test uses CPU platform ----------


def test_integration_test_uses_cpu_platform() -> None:
    """Integration test source must pass platform_name='CPU' to campaign functions."""

    test_source = (
        Path(__file__).resolve().parents[0] / "test_integration.py"
    ).read_text(encoding="utf-8")
    assert 'platform_name="CPU"' in test_source or "platform_name='CPU'" in test_source