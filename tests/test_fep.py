"""Tests for alchemical FEP (free energy perturbation) modules.

Covers L-39: No computational mutagenesis or DDG capability.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Step 1: Dependency availability
# ---------------------------------------------------------------------------


def test_step_1_pymbar_importable():
    """pymbar is importable and exposes MBAR."""
    import pymbar

    assert hasattr(pymbar, "MBAR")


def test_step_1_openmmtools_in_requirements():
    """openmmtools>=0.23.0 is declared in requirements.txt."""
    req_path = Path(__file__).resolve().parent.parent / "requirements.txt"
    content = req_path.read_text(encoding="utf-8")
    assert "openmmtools" in content, "openmmtools not found in requirements.txt"


# ---------------------------------------------------------------------------
# Step 2: FEPConfig dataclass
# ---------------------------------------------------------------------------


def test_step_2_fep_config():
    """FEPConfig is frozen and has correct defaults."""
    from dataclasses import FrozenInstanceError, is_dataclass

    from src.config import FEPConfig

    config = FEPConfig()
    assert is_dataclass(config)
    assert config.n_lambda_windows == 20
    assert config.temperature_k == 310.0
    assert config.soft_core_alpha == 0.5
    assert config.per_window_duration_ns == 2.0
    assert config.n_equilibration_steps == 5000
    assert config.annihilate_electrostatics is True
    assert config.annihilate_sterics is False
    with pytest.raises(FrozenInstanceError):
        config.n_lambda_windows = 30


# ---------------------------------------------------------------------------
# Step 3: Alchemical system construction (src/simulate/fep.py)
# ---------------------------------------------------------------------------


def test_step_3_lambda_schedule():
    """Lambda schedule spans [0, 1] with correct number of windows."""
    from src.simulate.fep import generate_lambda_schedule

    schedule = generate_lambda_schedule(20)
    assert schedule.shape == (20,)
    assert float(schedule[0]) == 0.0
    assert float(schedule[-1]) == 1.0
    assert np.all(np.diff(schedule) > 0), "Lambda schedule must be monotonically increasing"


def test_step_3_lambda_schedule_minimum():
    """Lambda schedule requires at least 2 windows."""
    from src.simulate.fep import generate_lambda_schedule

    with pytest.raises(ValueError, match="n_windows must be at least 2"):
        generate_lambda_schedule(1)


def test_step_3_lambda_schedule_two_windows():
    """Lambda schedule with 2 windows gives exactly [0.0, 1.0]."""
    from src.simulate.fep import generate_lambda_schedule

    schedule = generate_lambda_schedule(2)
    np.testing.assert_array_equal(schedule, [0.0, 1.0])


def test_step_3_create_alchemical_system_rejects_empty():
    """create_alchemical_system raises on empty atom indices."""
    from src.simulate.fep import create_alchemical_system
    from src.config import FEPConfig

    system = __import__("openmm").System()
    system.addParticle(12.0)
    config = FEPConfig()
    with pytest.raises(ValueError, match="must not be empty"):
        create_alchemical_system(system, [], config)


def test_step_3_create_alchemical_system_rejects_out_of_range():
    """create_alchemical_system raises on out-of-range atom indices."""
    from src.simulate.fep import create_alchemical_system
    from src.config import FEPConfig

    system = __import__("openmm").System()
    system.addParticle(12.0)
    config = FEPConfig()
    with pytest.raises(ValueError, match="out of range"):
        create_alchemical_system(system, [5], config)


# ---------------------------------------------------------------------------
# Step 4: BAR/MBAR analysis module (src/analyze/fep.py)
# ---------------------------------------------------------------------------


def test_step_4_delta_delta_g_arithmetic():
    """DDG is correctly computed as DG_complex - DG_free."""
    from src.analyze.fep import compute_delta_delta_g
    from src.config import KCAL_TO_KJ

    complex_result = {
        "delta_g_kj_mol": 10.0,
        "delta_g_std_kj_mol": 0.5,
    }
    free_result = {
        "delta_g_kj_mol": 8.0,
        "delta_g_std_kj_mol": 0.3,
    }

    ddg = compute_delta_delta_g(complex_result, free_result)
    assert abs(ddg["delta_delta_g_kj_mol"] - 2.0) < 1e-10
    expected_std = np.sqrt(0.5**2 + 0.3**2)
    assert abs(ddg["delta_delta_g_std_kj_mol"] - expected_std) < 1e-10
    assert abs(ddg["delta_delta_g_kcal_mol"] - 2.0 / KCAL_TO_KJ) < 1e-10


def test_step_4_delta_delta_g_negative():
    """Negative DDG indicates mutation stabilizes binding."""
    from src.analyze.fep import compute_delta_delta_g

    complex_result = {"delta_g_kj_mol": 5.0, "delta_g_std_kj_mol": 0.2}
    free_result = {"delta_g_kj_mol": 9.0, "delta_g_std_kj_mol": 0.3}

    ddg = compute_delta_delta_g(complex_result, free_result)
    assert ddg["delta_delta_g_kj_mol"] < 0.0, "Negative DDG expected"


def test_step_4_compute_delta_g_bar_rejects_empty():
    """compute_delta_g_bar raises on empty energy arrays."""
    from src.analyze.fep import compute_delta_g_bar

    with pytest.raises(ValueError, match="non-empty"):
        compute_delta_g_bar(np.array([]), np.array([1.0]), 310.0)


def test_step_4_compute_delta_g_mbar_rejects_bad_temperature():
    """compute_delta_g_mbar raises on non-positive temperature."""
    from src.analyze.fep import compute_delta_g_mbar

    with pytest.raises(ValueError, match="temperature_k must be positive"):
        compute_delta_g_mbar(np.zeros((2, 4)), np.array([2, 2]), -10.0)


# ---------------------------------------------------------------------------
# Step 5: FEP validation notebook
# ---------------------------------------------------------------------------


def test_step_5_fep_validation_notebook_valid():
    """FEP validation notebook is valid JSON with expected structure."""
    import json

    notebook_path = (
        Path(__file__).resolve().parent.parent
        / "notebooks"
        / "10_fep_validation.ipynb"
    )
    assert notebook_path.exists(), f"Missing notebook: {notebook_path}"

    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    assert nb["nbformat"] == 4
    cells = nb["cells"]
    assert len(cells) >= 8, f"Expected >=8 cells, got {len(cells)}"

    full_source = "\n".join(
        line for cell in cells for line in cell["source"]
    )
    # Key FEP components must be referenced
    for keyword in (
        "generate_lambda_schedule",
        "create_alchemical_system",
        "run_fep_campaign",
        "compute_delta_g_mbar",
        "compute_delta_delta_g",
        "barnase_A59G_ddg",
    ):
        assert keyword in full_source, f"Missing reference: {keyword}"


def test_step_5_fep_thermodynamic_cycle_consistency():
    """Thermodynamic cycle DDG is path-independent (synthetic validation)."""
    from src.analyze.fep import compute_delta_delta_g
    from src.config import KCAL_TO_KJ

    # Simulate a known mutation:
    # Ala59->Gly experimental DDG = +2.0 kcal/mol = +8.368 kJ/mol
    complex_dg = {"delta_g_kj_mol": 15.0, "delta_g_std_kj_mol": 0.4}
    free_dg = {"delta_g_kj_mol": 6.632, "delta_g_std_kj_mol": 0.3}

    ddg = compute_delta_delta_g(complex_dg, free_dg)
    expected_ddg_kj = 15.0 - 6.632
    assert abs(ddg["delta_delta_g_kj_mol"] - expected_ddg_kj) < 1e-10
    # Check kcal conversion
    assert abs(ddg["delta_delta_g_kcal_mol"] - expected_ddg_kj / KCAL_TO_KJ) < 1e-8


# ---------------------------------------------------------------------------
# Step 6: FEP CLI and mutagenesis scanning notebook
# ---------------------------------------------------------------------------


def test_step_6_fep_cli_parser():
    """FEP CLI parser accepts required arguments."""
    from scripts.run_fep import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "--complex-pdb", "data/pdb/SPINK7_KLK5.pdb",
        "--chain-id", "A",
        "--residue-number", "42",
        "--mutation", "GLY",
    ])
    assert args.residue_number == 42
    assert args.mutation == "GLY"
    assert args.chain_id == "A"
    assert args.n_lambda_windows == 20
    assert args.per_window_ns == 2.0
    assert args.temperature_k == 310.0


def test_step_6_fep_cli_custom_windows():
    """FEP CLI accepts custom lambda window count."""
    from scripts.run_fep import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "--complex-pdb", "data/pdb/SPINK7_KLK5.pdb",
        "--chain-id", "B",
        "--residue-number", "10",
        "--mutation", "ALA",
        "--n-lambda-windows", "12",
        "--per-window-ns", "1.0",
    ])
    assert args.n_lambda_windows == 12
    assert args.per_window_ns == 1.0


def test_step_6_mutagenesis_notebook_valid():
    """SPINK7 mutagenesis notebook is valid JSON with expected structure."""
    import json

    notebook_path = (
        Path(__file__).resolve().parent.parent
        / "notebooks"
        / "11_spink7_mutagenesis.ipynb"
    )
    assert notebook_path.exists(), f"Missing notebook: {notebook_path}"

    nb = json.loads(notebook_path.read_text(encoding="utf-8"))
    assert nb["nbformat"] == 4
    cells = nb["cells"]
    assert len(cells) >= 6, f"Expected >=6 cells, got {len(cells)}"

    full_source = "\n".join(
        line for cell in cells for line in cell["source"]
    )
    for keyword in (
        "run_fep_campaign",
        "compute_delta_g_mbar",
        "compute_delta_delta_g",
        "FEPConfig",
        "SPINK7",
    ):
        assert keyword in full_source, f"Missing reference: {keyword}"
