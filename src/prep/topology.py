"""OpenMM topology construction utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
from pathlib import Path

import openmm
from openmm.app import ForceField, HBonds, Modeller, NoCutoff, PDBFile, Topology

from src import PhysicalValidityError
from src.config import SystemConfig


logger = logging.getLogger(__name__)


def _validate_inputs(pdb_path: Path, system_config: SystemConfig) -> Path:
    """Validate public API inputs for OpenMM topology construction."""

    validated_path = Path(pdb_path)
    if not validated_path.exists():
        raise FileNotFoundError(f"Input structure does not exist: {validated_path}")
    if validated_path.suffix.lower() != ".pdb":
        raise ValueError("pdb_path must have a .pdb extension")
    if not system_config.force_field.strip():
        raise ValueError("system_config.force_field must be a non-empty string")
    if not system_config.water_model.strip():
        raise ValueError("system_config.water_model must be a non-empty string")
    return validated_path


def _count_hydrogens(topology: Topology) -> int:
    """Count hydrogen atoms in an OpenMM topology."""

    return sum(1 for atom in topology.atoms() if atom.element is not None and atom.element.symbol == "H")


def build_topology(
    pdb_path: Path,
    system_config: SystemConfig,
) -> tuple[Topology, openmm.System, Modeller]:
    """Build an OpenMM topology, modeller, and AMBER-mapped system.

    Invariants: Enforces the constitution's explicit-hydrogen requirement for the
    all-atom system definition and fails fast if AMBER ff14SB parameter mapping
    cannot be completed without inventing parameters.

    Args:
        pdb_path: Protonated PDB path.
        system_config: Immutable system-preparation parameters.

    Returns:
        tuple[Topology, openmm.System, Modeller]: Parameterized OpenMM topology,
        system, and modeller. Positions in the returned modeller have shape
        [N_atoms, 3] in nm, and system particle masses are defined for [N_atoms].
    """

    validated_path = _validate_inputs(pdb_path, system_config)
    logger.info("Building OpenMM topology from %s", validated_path)

    try:
        pdb = PDBFile(str(validated_path))
    except Exception as exc:
        logger.error("Topology generation failed while parsing %s: %s", validated_path, exc)
        raise PhysicalValidityError("Input structure is not a valid non-empty PDB") from exc

    initial_topology = pdb.topology
    if initial_topology.getNumAtoms() == 0:
        logger.error("Topology generation failed: input structure contains zero atoms")
        raise PhysicalValidityError("Input structure contains zero atoms")

    force_field = ForceField(system_config.force_field, system_config.water_model)
    modeller = Modeller(initial_topology, pdb.positions)

    try:
        modeller.addHydrogens(force_field, pH=system_config.ph)
        system = force_field.createSystem(
            modeller.topology,
            constraints=HBonds,
            nonbondedMethod=NoCutoff,
            rigidWater=True,
        )
    except Exception as exc:  # pragma: no cover - exercised by integration with OpenMM internals
        logger.error("Topology parameterization failed for %s: %s", validated_path, exc)
        raise PhysicalValidityError(
            "AMBER topology mapping failed; structure contains unsupported or unmatched atoms"
        ) from exc

    topology = modeller.topology
    atom_count = topology.getNumAtoms()
    particle_count = system.getNumParticles()
    hydrogen_count = _count_hydrogens(topology)

    if particle_count != atom_count:
        logger.error(
            "Topology particle mismatch: topology atoms=%d, system particles=%d",
            atom_count,
            particle_count,
        )
        raise PhysicalValidityError(
            "Topology atom count does not match the parameterized system particle count"
        )

    if hydrogen_count == 0:
        logger.error("Explicit-hydrogen requirement violated for %s", validated_path)
        raise PhysicalValidityError("Parameterized topology must contain explicit hydrogens")

    logger.info(
        "Built topology with %d atoms, %d hydrogens, and %d system particles",
        atom_count,
        hydrogen_count,
        particle_count,
    )
    return topology, system, modeller