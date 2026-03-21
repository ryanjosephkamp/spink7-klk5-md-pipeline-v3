"""Solvation and ion-placement utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
from pathlib import Path

from openmm import unit
from openmm.app import ForceField, Modeller

from src import PhysicalValidityError
from src.config import SUPPORTED_BOX_SHAPES, SystemConfig


logger = logging.getLogger(__name__)

_WATER_RESIDUE_NAMES = {"HOH", "WAT"}


def _validate_inputs(modeller: Modeller, system_config: SystemConfig) -> None:
    """Validate public API inputs for solvation."""

    if modeller.topology.getNumAtoms() == 0:
        raise PhysicalValidityError("Cannot solvate an empty topology")
    if system_config.box_padding_nm <= 0.0:
        raise ValueError("system_config.box_padding_nm must be positive")
    if system_config.ionic_strength_molar < 0.0:
        raise ValueError("system_config.ionic_strength_molar must be non-negative")
    if not system_config.positive_ion.strip():
        raise ValueError("system_config.positive_ion must be a non-empty string")
    if not system_config.negative_ion.strip():
        raise ValueError("system_config.negative_ion must be a non-empty string")
    if system_config.box_shape not in SUPPORTED_BOX_SHAPES:
        raise ValueError(
            f"system_config.box_shape must be one of {SUPPORTED_BOX_SHAPES}, "
            f"got '{system_config.box_shape}'"
        )


def _water_model_name(water_model: str) -> str:
    """Map the configured water-model XML path to the OpenMM solvent model name.

    Consults SUPPORTED_WATER_MODELS for validated mappings.
    Falls back to stem extraction for unregistered models.
    """
    from src.config import SUPPORTED_WATER_MODELS

    for _key, (xml_path, solvent_name) in SUPPORTED_WATER_MODELS.items():
        if water_model == xml_path:
            return solvent_name
    return Path(water_model).stem


def _ion_residue_name(ion_name: str) -> str:
    """Convert an ion label like Na+ or Cl- to an OpenMM residue name."""

    return "".join(character for character in ion_name.upper() if character.isalpha())


def _count_solvent_species(
    modeller: Modeller,
    positive_ion_residue: str,
    negative_ion_residue: str,
) -> tuple[int, int, int]:
    """Count waters and monoatomic ions in the solvated topology."""

    water_count = 0
    positive_ion_count = 0
    negative_ion_count = 0

    for residue in modeller.topology.residues():
        residue_name = residue.name.upper()
        if residue_name in _WATER_RESIDUE_NAMES:
            water_count += 1
        elif residue_name == positive_ion_residue:
            positive_ion_count += 1
        elif residue_name == negative_ion_residue:
            negative_ion_count += 1

    return water_count, positive_ion_count, negative_ion_count


def solvate_system(
    modeller: Modeller,
    system_config: SystemConfig,
) -> tuple[Modeller, int, int, int]:
    """Add explicit solvent and counter-ions to an OpenMM modeller.

    Invariants: Enforces a periodic box definition for downstream PME dynamics and
    fails fast if no explicit water is added or if the solvated topology lacks box
    vectors, both prerequisites for satisfying IV-7 in later simulation stages.

    Args:
        modeller: Protein topology and positions. Positions have shape [N_atoms, 3] in nm.
        system_config: Immutable system-preparation parameters.

    Returns:
        tuple[Modeller, int, int, int]: Solvated modeller plus counts of waters,
        positive ions, and negative ions. Final positions have shape [N_total_atoms, 3] in nm.
    """

    _validate_inputs(modeller, system_config)

    force_field = ForceField(system_config.force_field, system_config.water_model)
    solvent_model = _water_model_name(system_config.water_model)
    positive_ion_residue = _ion_residue_name(system_config.positive_ion)
    negative_ion_residue = _ion_residue_name(system_config.negative_ion)
    openmm_box_shape = "cube" if system_config.box_shape == "cubic" else system_config.box_shape

    logger.info(
        "Adding solvent with model=%s, padding=%.3f nm, ionic strength=%.3f M, box_shape=%s",
        solvent_model,
        system_config.box_padding_nm,
        system_config.ionic_strength_molar,
        system_config.box_shape,
    )

    try:
        modeller.addSolvent(
            force_field,
            model=solvent_model,
            padding=system_config.box_padding_nm * unit.nanometer,
            ionicStrength=system_config.ionic_strength_molar * unit.molar,
            positiveIon=system_config.positive_ion,
            negativeIon=system_config.negative_ion,
            neutralize=True,
            boxShape=openmm_box_shape,
        )
    except Exception as exc:  # pragma: no cover - exercised by OpenMM internals
        logger.error("Solvation failed: %s", exc)
        raise PhysicalValidityError("OpenMM solvation failed for the provided topology") from exc

    water_count, positive_ion_count, negative_ion_count = _count_solvent_species(
        modeller,
        positive_ion_residue,
        negative_ion_residue,
    )
    box_vectors = modeller.topology.getPeriodicBoxVectors()

    if box_vectors is None:
        logger.error("Solvation failed to define periodic box vectors")
        raise PhysicalValidityError("Solvated topology must define periodic box vectors")

    if water_count == 0:
        logger.error("Solvation failed to add any explicit water molecules")
        raise PhysicalValidityError("Solvated topology must contain explicit water molecules")

    logger.info(
        "Solvated topology contains %d waters, %d %s ions, and %d %s ions",
        water_count,
        positive_ion_count,
        positive_ion_residue,
        negative_ion_count,
        negative_ion_residue,
    )
    return modeller, water_count, positive_ion_count, negative_ion_count