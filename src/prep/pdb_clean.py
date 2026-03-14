"""PDB cleaning utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List

from src import PipelineError


logger = logging.getLogger(__name__)

_WATER_RESIDUES = {"HOH", "WAT", "TIP", "TIP3", "TIP3P"}


def _validate_pdb_path(pdb_path: Path) -> Path:
    """Validate the public PDB path input."""

    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"Input structure does not exist: {path}")
    if path.suffix.lower() not in {".pdb", ".cif"}:
        raise ValueError("pdb_path must have a .pdb or .cif extension")
    if path.suffix.lower() == ".cif":
        raise PipelineError("clean_structure currently supports .pdb inputs only")
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
) -> Path:
    """Clean a PDB by selecting chains and removing optional exclusions.

    Invariants: None.

    Args:
        pdb_path: Input PDB path.
        chains_to_keep: Chain identifiers to retain.
        remove_heteroatoms: Whether to remove non-protein HETATM records.
        remove_waters: Whether to remove water molecules.

    Returns:
        Path: Path to the cleaned PDB file.
    """

    validated_path = _validate_pdb_path(pdb_path)
    normalized_chains = {chain.strip() for chain in chains_to_keep if chain.strip()}
    if not normalized_chains:
        raise ValueError("chains_to_keep must contain at least one non-empty chain identifier")

    output_dir = _prepared_output_dir(validated_path)
    output_path = output_dir / f"{validated_path.stem}_cleaned.pdb"

    cleaned_lines: list[str] = []
    # Track whether the last ATOM/HETATM was retained so that
    # TER records are only emitted after kept coordinate blocks.
    previous_atom_kept = False

    with validated_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = line[0:6].strip()

            # --- Coordinate records: filter by chain, heteroatom, and water rules ---
            if record in {"ATOM", "HETATM"}:
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
                if previous_atom_kept:
                    cleaned_lines.append(line)
                previous_atom_kept = False
                continue

            # --- MODEL/ENDMDL: always preserve multi-model structure ---
            if record in {"MODEL", "ENDMDL"}:
                cleaned_lines.append(line)
                previous_atom_kept = False

    if not cleaned_lines:
        raise PipelineError("Cleaning removed all structural records from the input PDB")

    if not cleaned_lines[-1].startswith("END"):
        cleaned_lines.append("END\n")

    with output_path.open("w", encoding="utf-8") as handle:
        handle.writelines(cleaned_lines)

    logger.info("Cleaned PDB written to %s", output_path)
    return output_path
