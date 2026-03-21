"""Shared topology I/O utilities for simulation modules."""

from __future__ import annotations

import logging
from pathlib import Path

from openmm.app import PDBFile, Simulation


logger = logging.getLogger(__name__)


def save_topology_pdb(simulation: Simulation, output_path: Path) -> Path:
    """Persist the simulation topology as a PDB snapshot.

    Args:
        simulation: OpenMM simulation with initialized context.
        output_path: Destination PDB file path.

    Returns:
        The written topology file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with output_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(simulation.topology, state.getPositions(), handle)
    logger.info("Topology saved to %s", output_path)
    return output_path
