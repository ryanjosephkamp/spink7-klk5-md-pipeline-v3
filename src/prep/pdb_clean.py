"""PDB cleaning utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import List

try:
    import gemmi
    _HAS_GEMMI = True
except ImportError:
    _HAS_GEMMI = False

from src import PipelineError


logger = logging.getLogger(__name__)

_WATER_RESIDUES = {"HOH", "WAT", "TIP", "TIP3", "TIP3P"}


def _parse_cif_to_pdb_lines(cif_path: Path) -> list[str]:
    """Parse a CIF/mmCIF file and convert atoms to PDB-format lines.

    Uses Gemmi to read the CIF file and produce PDB-compatible fixed-width
    lines for each coordinate record. This allows the existing PDB filtering
    logic to operate on CIF-sourced data without modification.

    Args:
        cif_path: Path to a .cif or .mmcif file.

    Returns:
        PDB-format lines including ATOM/HETATM/TER/END records.

    Raises:
        PipelineError: If Gemmi is not installed or the CIF file cannot be parsed.
    """
    if not _HAS_GEMMI:
        raise PipelineError(
            "Gemmi is required for CIF/mmCIF support. "
            "Install it with: pip install gemmi"
        )

    try:
        structure = gemmi.read_structure(str(cif_path))
    except Exception as exc:
        raise PipelineError(f"Failed to parse CIF file {cif_path}: {exc}") from exc

    if len(structure) == 0:
        raise PipelineError(f"CIF file {cif_path} contains no models")

    pdb_string = structure.make_pdb_string()
    return pdb_string.splitlines(keepends=True)


def _validate_pdb_path(pdb_path: Path) -> Path:
    """Validate the public PDB path input."""

    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"Input structure does not exist: {path}")
    if path.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
        raise ValueError("pdb_path must have a .pdb, .cif, or .mmcif extension")
    return path


def _prepared_output_dir(pdb_path: Path) -> Path:
    """Infer the prepared output directory from the raw input path."""

    if pdb_path.parent.name == "raw":
        output_dir = pdb_path.parent.parent / "prepared"
    else:
        output_dir = pdb_path.parent / "prepared"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _should_keep_atom_line(
    line: str,
    chains_to_keep: set[str],
    remove_heteroatoms: bool,
    remove_waters: bool,
) -> bool:
    """Decide whether an ATOM/HETATM record should be preserved."""

    record = line[0:6].strip()
    chain_id = line[21].strip()
    residue_name = line[17:20].strip().upper()

    if chain_id not in chains_to_keep:
        return False
    if record == "HETATM" and remove_heteroatoms:
        return False
    if residue_name in _WATER_RESIDUES and remove_waters:
        return False
    return True


def clean_structure(
    pdb_path: Path,
    chains_to_keep: List[str],
    remove_heteroatoms: bool = True,
    remove_waters: bool = True,
    model_index: int | None = 1,
) -> Path:
    """Clean a PDB or CIF structure by selecting chains and removing optional exclusions.

    CIF/mmCIF inputs are converted to PDB-format lines via Gemmi before
    applying the same filtering logic.  Output is always a cleaned ``.pdb``
    file.

    Invariants: None.

    Args:
        pdb_path: Input PDB, CIF, or mmCIF path.
        chains_to_keep: Chain identifiers to retain.
        remove_heteroatoms: Whether to remove non-protein HETATM records.
        remove_waters: Whether to remove water molecules.
        model_index: 1-based model index to select from multi-model PDB files
            (e.g., NMR ensembles). Defaults to 1 (first model). Set to None
            to retain all models (not recommended for MD starting structures).

    Returns:
        Path: Path to the cleaned PDB file.
    """

    validated_path = _validate_pdb_path(pdb_path)
    normalized_chains = {chain.strip() for chain in chains_to_keep if chain.strip()}
    if not normalized_chains:
        raise ValueError("chains_to_keep must contain at least one non-empty chain identifier")

    output_dir = _prepared_output_dir(validated_path)
    output_path = output_dir / f"{validated_path.stem}_cleaned.pdb"

    is_cif = validated_path.suffix.lower() in {".cif", ".mmcif"}

    if is_cif:
        logger.info("Detected CIF/mmCIF input; converting via Gemmi for cleaning")
        source_lines = _parse_cif_to_pdb_lines(validated_path)
    else:
        with validated_path.open("r", encoding="utf-8") as handle:
            source_lines = handle.readlines()

    cleaned_lines: list[str] = []
    # Track whether the last ATOM/HETATM was retained so that
    # TER records are only emitted after kept coordinate blocks.
    previous_atom_kept = False

    # Multi-model (NMR ensemble) tracking state.
    current_model: int | None = None
    model_count = 0
    is_multi_model = False
    # If model_index is None, accept all models; otherwise select the target.
    target_model = model_index

    for line in source_lines:
            record = line[0:6].strip()

            # --- MODEL/ENDMDL: track multi-model structure and select target model ---
            if record == "MODEL":
                model_count += 1
                is_multi_model = True
                current_model = model_count
                if target_model is None or current_model == target_model:
                    cleaned_lines.append(line)
                previous_atom_kept = False
                continue

            if record == "ENDMDL":
                if target_model is None or current_model == target_model:
                    cleaned_lines.append(line)
                previous_atom_kept = False
                continue

            # --- Coordinate records: filter by chain, heteroatom, and water rules ---
            if record in {"ATOM", "HETATM"}:
                # Skip atoms from non-selected models in multi-model files.
                if is_multi_model and target_model is not None and current_model != target_model:
                    previous_atom_kept = False
                    continue
                keep_line = _should_keep_atom_line(
                    line,
                    normalized_chains,
                    remove_heteroatoms=remove_heteroatoms,
                    remove_waters=remove_waters,
                )
                if keep_line:
                    cleaned_lines.append(line)
                previous_atom_kept = keep_line
                continue

            # --- TER records: emit only when they follow a kept atom block ---
            if record == "TER":
                if is_multi_model and target_model is not None and current_model != target_model:
                    previous_atom_kept = False
                    continue
                if previous_atom_kept:
                    cleaned_lines.append(line)
                previous_atom_kept = False
                continue

    if is_multi_model and target_model is not None and target_model > model_count:
        raise PipelineError(
            f"Requested model_index={target_model} but PDB contains only {model_count} model(s)"
        )

    if is_multi_model and model_count > 1 and target_model is not None:
        warnings.warn(
            f"Multi-model PDB detected ({model_count} models). "
            f"Selected model {target_model}; discarded {model_count - 1} other model(s). "
            f"This is standard practice for MD starting structures from NMR ensembles.",
            stacklevel=2,
        )
        logger.info(
            "Multi-model PDB: selected model %d of %d", target_model, model_count
        )

    if not cleaned_lines:
        raise PipelineError("Cleaning removed all structural records from the input PDB")

    if not cleaned_lines[-1].startswith("END"):
        cleaned_lines.append("END\n")

    with output_path.open("w", encoding="utf-8") as handle:
        handle.writelines(cleaned_lines)

    logger.info("Cleaned PDB written to %s", output_path)
    return output_path
