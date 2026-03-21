"""OpenMM topology construction utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import openmm
from openmm.app import ForceField, HBonds, Modeller, NoCutoff, PDBFile, PME, Topology

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


def _is_pre_protonated(pdb_path: Path) -> bool:
    """Detect whether the input PDB was produced by assign_protonation().

    Detection criteria:
    1. Filename ends with '_protonated.pdb' (the naming convention from
       assign_protonation()).
    2. The PDB contains at least one hydrogen atom (element symbol 'H' in
       the PDB ATOM records).

    Both criteria must be met to return True.
    """

    if not pdb_path.stem.endswith("_protonated"):
        return False

    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line[0:6].strip() in {"ATOM", "HETATM"}:
                element = line[76:78].strip().upper()
                if element == "H":
                    return True
    return False


def build_topology(
    pdb_path: Path,
    system_config: SystemConfig,
    nonbonded_method: object = NoCutoff,
    nonbonded_cutoff_nm: float = 1.0,
    skip_hydrogens: bool = False,
) -> tuple[Topology, openmm.System, Modeller]:
    """Build an OpenMM topology, modeller, and AMBER-mapped system.

    Invariants: Enforces the constitution's explicit-hydrogen requirement for the
    all-atom system definition and fails fast if AMBER ff14SB parameter mapping
    cannot be completed without inventing parameters.

    Args:
        pdb_path: Protonated PDB path.
        system_config: Immutable system-preparation parameters.
        nonbonded_method: OpenMM nonbonded method for the returned System.
            Defaults to NoCutoff because the pre-solvation topology lacks a
            periodic box.  The system must be re-parameterized with PME after
            solvation for production dynamics.  A UserWarning is emitted when
            NoCutoff is used on systems with more than 500 atoms.
        nonbonded_cutoff_nm: Cutoff distance in nm for cutoff-based nonbonded
            methods (PME, CutoffPeriodic, CutoffNonPeriodic). Ignored when
            nonbonded_method is NoCutoff. Defaults to 1.0 nm per
            global_constitution.md section 2.3.
        skip_hydrogens: If True, skip the modeller.addHydrogens() call,
            preserving existing protonation states from upstream preparation.
            Set to True when the input PDB has already been protonated by
            assign_protonation() or an equivalent tool.

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

    if not skip_hydrogens and _is_pre_protonated(validated_path):
        logger.info(
            "Detected pre-protonated input (%s); automatically skipping addHydrogens",
            validated_path.name,
        )
        skip_hydrogens = True

    try:
        initial_h_count = _count_hydrogens(modeller.topology)

        if skip_hydrogens:
            logger.info("Skipping addHydrogens — input structure is pre-protonated")
        else:
            modeller.addHydrogens(force_field, pH=system_config.ph)
            if initial_h_count > 0:
                post_h_count = _count_hydrogens(modeller.topology)
                warnings.warn(
                    f"addHydrogens() was called on a structure that already "
                    f"contained {initial_h_count} hydrogen atoms. "
                    f"Post-addHydrogens count: {post_h_count}. If the input "
                    f"was pre-protonated, set skip_hydrogens=True to avoid "
                    f"double protonation.",
                    UserWarning,
                    stacklevel=2,
                )

        create_system_kwargs: dict[str, object] = {
            "constraints": HBonds,
            "nonbondedMethod": nonbonded_method,
            "rigidWater": True,
        }
        if nonbonded_method is not NoCutoff:
            from openmm import unit as omm_unit
            create_system_kwargs["nonbondedCutoff"] = (
                nonbonded_cutoff_nm * omm_unit.nanometer
            )

        system = force_field.createSystem(modeller.topology, **create_system_kwargs)
    except Exception as exc:  # pragma: no cover - exercised by integration with OpenMM internals
        logger.error("Topology parameterization failed for %s: %s", validated_path, exc)
        raise PhysicalValidityError(
            "AMBER topology mapping failed; structure contains unsupported or unmatched atoms"
        ) from exc

    topology = modeller.topology
    atom_count = topology.getNumAtoms()

    logger.info("Nonbonded method for returned system: %s", type(nonbonded_method).__name__)

    if nonbonded_method is NoCutoff and atom_count > 500:
        warnings.warn(
            f"build_topology() created a system with {atom_count} atoms "
            f"using NoCutoff nonbonded method. For systems of this size, "
            f"NoCutoff produces physically inaccurate long-range "
            f"electrostatics. Ensure the system is re-parameterized with "
            f"PME after solvation, or pass nonbonded_method=PME if a "
            f"periodic box is already defined.",
            UserWarning,
            stacklevel=2,
        )
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