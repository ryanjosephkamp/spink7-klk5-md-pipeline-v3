"""Tests for inter-module data flow contracts (L-37)."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

# Import the fixed loader from scripts/run_analysis.py.
# The scripts directory is importable as a side effect of the project
# structure; we access the private function directly for unit testing.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from run_analysis import _load_work_values


class TestLoadWorkValuesCSV:
    """Step 1: _load_work_values correctly handles SMD work CSV files."""

    def test_extracts_work_column_from_two_column_csv(self, tmp_path: Path) -> None:
        """Work column is extracted without interleaving time values."""

        csv_path = tmp_path / "smd_work.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["time_ps", "work_kj_mol"])
            writer.writerow([1.0, 10.5])
            writer.writerow([2.0, 21.3])
            writer.writerow([3.0, 32.7])

        work_values = _load_work_values(csv_path)

        assert work_values.shape == (3,), f"Expected shape (3,), got {work_values.shape}"
        np.testing.assert_allclose(work_values, [10.5, 21.3, 32.7])
        assert not np.any(np.isclose(work_values, [1.0, 2.0, 3.0])), \
            "Time values leaked into work array"

    def test_single_row_csv(self, tmp_path: Path) -> None:
        """Single-row CSV files are handled correctly."""

        csv_path = tmp_path / "smd_work_single.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["time_ps", "work_kj_mol"])
            writer.writerow([5.0, 42.0])

        work_values = _load_work_values(csv_path)

        assert work_values.shape == (1,)
        np.testing.assert_allclose(work_values, [42.0])

    def test_npy_format_still_works(self, tmp_path: Path) -> None:
        """NPY loading path is not broken by the CSV fix."""

        npy_path = tmp_path / "work.npy"
        expected = np.array([10.0, 20.0, 30.0])
        np.save(npy_path, expected)

        work_values = _load_work_values(npy_path)

        np.testing.assert_allclose(work_values, expected)

    def test_five_row_csv_matches_expected(self, tmp_path: Path) -> None:
        """Realistic multi-row CSV matches expected work values exactly."""

        csv_path = tmp_path / "smd_work_5.csv"
        expected_work = [10.5, 21.3, 32.7, 44.1, 55.9]
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["time_ps", "work_kj_mol"])
            for i, w in enumerate(expected_work, 1):
                writer.writerow([float(i), w])

        work_values = _load_work_values(csv_path)

        assert work_values.shape == (5,)
        np.testing.assert_allclose(work_values, expected_work)


class TestTopologyCoSaving:
    """Step 2: Production MD co-saves topology PDB alongside DCD trajectory."""

    def test_save_topology_pdb_creates_valid_file(self, tmp_path: Path) -> None:
        """The shared save_topology_pdb utility creates a valid PDB file."""

        import openmm
        from openmm import unit
        from openmm.app import Element, Simulation, Topology

        from src.simulate._topology_io import save_topology_pdb

        topology = Topology()
        chain = topology.addChain()
        residue = topology.addResidue("ALA", chain)
        carbon = Element.getByAtomicNumber(6)
        topology.addAtom("CA", carbon, residue)

        system = openmm.System()
        system.addParticle(12.0)
        system.setDefaultPeriodicBoxVectors(
            openmm.Vec3(3, 0, 0) * unit.nanometer,
            openmm.Vec3(0, 3, 0) * unit.nanometer,
            openmm.Vec3(0, 0, 3) * unit.nanometer,
        )

        integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
        platform = openmm.Platform.getPlatformByName("CPU")
        sim = Simulation(topology, system, integrator, platform)
        sim.context.setPositions([[0.5, 0.5, 0.5]] * unit.nanometer)

        pdb_path = tmp_path / "test_topology.pdb"
        result = save_topology_pdb(sim, pdb_path)

        assert result == pdb_path
        assert pdb_path.exists()
        assert pdb_path.stat().st_size > 0
        content = pdb_path.read_text()
        assert "ATOM" in content or "HETATM" in content
