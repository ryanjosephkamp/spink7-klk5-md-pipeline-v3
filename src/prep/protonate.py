"""Protonation-state assignment utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
import warnings
from collections import OrderedDict
from pathlib import Path

from propka.molecular_container import MolecularContainer
from propka.run import single as propka_single

from src import PipelineError


logger = logging.getLogger(__name__)


_PROTONATION_RULES: dict[str, tuple[float, str, str, list[str], list[str]]] = {
    "ASP": (4.0, "ASH", "ASP", ["HD2"], []),
    "GLU": (4.4, "GLH", "GLU", ["HE2"], []),
    "HIS": (6.5, "HIP", "HIE", ["HD1", "HE2"], ["HE2"]),
    "LYS": (10.5, "LYS", "LYN", ["HZ1", "HZ2", "HZ3"], []),
    "CYS": (8.3, "CYS", "CYM", ["HG"], []),
    "TYR": (10.1, "TYR", "TYM", ["HH"], []),
    "ARG": (12.0, "ARG", "ARG", ["HH11", "HH12", "HH21", "HH22"], ["HH11", "HH12", "HH21", "HH22"]),
    # Phosphorylated residues — pKa values for the phosphate second dissociation
    "SEP": (5.8, "SEP", "SEP", [], []),   # Phosphoserine
    "TPO": (6.3, "TPO", "TPO", [], []),   # Phosphothreonine
    "PTR": (5.8, "PTR", "PTR", [], []),   # Phosphotyrosine
}

# PROPKA residue type names that map to _PROTONATION_RULES keys.
_PROPKA_TYPE_MAP: dict[str, str] = {
    "ASP": "ASP",
    "GLU": "GLU",
    "HIS": "HIS",
    "LYS": "LYS",
    "CYS": "CYS",
    "TYR": "TYR",
    "ARG": "ARG",
}


def _run_propka(pdb_path: Path) -> dict[tuple[str, str, str], float]:
    """Run PROPKA to compute environment-dependent pKa values.

    Returns a mapping of (chain_id, residue_number, residue_name) to the
    PROPKA-predicted effective pKa.  Returns an empty dict if PROPKA fails,
    allowing the caller to fall back to reference pKa values.
    """
    try:
        molecule: MolecularContainer = propka_single(str(pdb_path), write_pka=False)
    except Exception:
        logger.warning(
            "PROPKA failed to process %s; falling back to reference pKa values",
            pdb_path,
        )
        return {}

    pka_dict: dict[tuple[str, str, str], float] = {}

    # Use the averaged conformation ("AVR") when available; fall back to first.
    if "AVR" in molecule.conformations:
        conformation = molecule.conformations["AVR"]
    else:
        conformations = list(molecule.conformations.values())
        if not conformations:
            return {}
        conformation = conformations[0]

    for group in conformation.groups:
        residue_type = group.residue_type
        if residue_type not in _PROPKA_TYPE_MAP:
            continue
        canonical_name = _PROPKA_TYPE_MAP[residue_type]
        chain_id = str(group.atom.chain_id).strip()
        res_num = str(group.atom.res_num).strip()
        pka_dict[(chain_id, res_num, canonical_name)] = group.pka_value

    return pka_dict


def _validate_inputs(pdb_path: Path, ph: float, force_field: str) -> Path:
    """Validate public API inputs for protonation assignment."""

    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"Input structure does not exist: {path}")
    if path.suffix.lower() != ".pdb":
        raise ValueError("pdb_path must have a .pdb extension")
    if ph <= 0.0:
        raise ValueError("ph must be positive")
    if force_field.strip().upper() != "AMBER":
        raise ValueError("force_field must be 'AMBER'")
    return path


def _prepared_output_dir(pdb_path: Path) -> Path:
    """Infer the prepared output directory from the cleaned input path."""

    if pdb_path.parent.name == "prepared":
        output_dir = pdb_path.parent
    elif pdb_path.parent.name == "raw":
        output_dir = pdb_path.parent.parent / "prepared"
    else:
        output_dir = pdb_path.parent / "prepared"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _parse_atom_line(line: str) -> dict[str, object]:
    """Parse a PDB ATOM/HETATM line into a structured dictionary."""

    return {
        "record": line[0:6].strip(),
        "atom_name": line[12:16].strip(),
        "residue_name": line[17:20].strip(),
        "chain_id": line[21],
        "residue_seq": line[22:26],
        "insertion_code": line[26],
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "occupancy": line[54:60].strip() or "1.00",
        "temp_factor": line[60:66].strip() or "20.00",
        "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
    }


def _format_atom_line(serial: int, atom: dict[str, object]) -> str:
    """Format a structured atom dictionary as a PDB ATOM line."""

    return (
        f"{str(atom['record']).ljust(6)}{serial:5d} "
        f"{str(atom['atom_name']).rjust(4)} "
        f"{str(atom['residue_name']).rjust(3)} {str(atom['chain_id'])}"
        f"{str(atom['residue_seq']).rjust(4)}{str(atom['insertion_code'])}   "
        f"{float(atom['x']):8.3f}{float(atom['y']):8.3f}{float(atom['z']):8.3f}"
        f"{float(atom['occupancy']):6.2f}{float(atom['temp_factor']):6.2f}          "
        f"{str(atom['element']).rjust(2)}\n"
    )


def _protonation_decision(
    residue_name: str,
    ph: float,
    pka_override: float | None = None,
) -> tuple[str, list[str], str] | None:
    """Determine the protonation-state decision for a residue if titratable.

    When *pka_override* is provided (e.g., from PROPKA), it replaces the
    tabulated reference pKa for the Henderson-Hasselbalch comparison.
    When *pka_override* is ``None``, the tabulated reference pKa is used
    (backward-compatible behavior).
    """

    residue_key = residue_name.upper()
    if residue_key not in _PROTONATION_RULES:
        return None

    pka_ref, protonated_name, deprotonated_name, protonated_hydrogens, deprotonated_hydrogens = _PROTONATION_RULES[residue_key]

    if pka_override is not None:
        pka_used = pka_override
        source = "PROPKA"
    else:
        pka_used = pka_ref
        source = "ref"

    # Non-standard titratable residues (e.g., phosphorylated) have empty
    # hydrogen lists and require no hydrogen modification.
    if not protonated_hydrogens and not deprotonated_hydrogens:
        target_name = protonated_name if ph < pka_used else deprotonated_name
        return target_name, [], f"non-standard titratable residue (pKa {pka_used:.2f}); no hydrogen modification"

    if ph < pka_used:
        return protonated_name, protonated_hydrogens, f"pH {ph:.2f} < pKa({source}) {pka_used:.2f}"
    return deprotonated_name, deprotonated_hydrogens, f"pH {ph:.2f} >= pKa({source}) {pka_used:.2f}"


def _centroid(atom_records: list[dict[str, object]]) -> tuple[float, float, float]:
    """Compute a residue centroid for synthetic hydrogen placement."""

    heavy_atoms = [atom for atom in atom_records if str(atom["element"]).upper() != "H"]
    atoms = heavy_atoms or atom_records
    n_atoms = float(len(atoms))
    return (
        sum(float(atom["x"]) for atom in atoms) / n_atoms,
        sum(float(atom["y"]) for atom in atoms) / n_atoms,
        sum(float(atom["z"]) for atom in atoms) / n_atoms,
    )


def assign_protonation(
    pdb_path: Path,
    ph: float = 7.4,
    force_field: str = "AMBER",
    use_propka: bool = True,
) -> Path:
    """Assign heuristic AMBER-compatible protonation states and add missing hydrogens.

    When *use_propka* is ``True`` (the default), PROPKA is run on the input
    structure to obtain environment-corrected pKa values.  When ``False`` or
    when PROPKA fails, the tabulated reference Henderson-Hasselbalch pKa
    values are used instead.

    Invariants: None.

    Args:
        pdb_path: Cleaned PDB path.
        ph: Target pH for protonation-state assignment.
        force_field: Force-field family used for residue naming.
        use_propka: If ``True``, use PROPKA for environment-dependent pKa
            prediction.  If ``False``, use reference pKa values only.

    Returns:
        Path: Path to the protonated PDB file.
    """

    validated_path = _validate_inputs(pdb_path, ph, force_field)
    output_dir = _prepared_output_dir(validated_path)
    output_path = output_dir / f"{validated_path.stem}_protonated.pdb"
    log_path = output_dir / f"{validated_path.stem}_protonation.log"

    # Group atoms by residue (chain, seq, icode) to apply per-residue
    # protonation-state decisions while preserving input ordering.
    ordered_residues: OrderedDict[tuple[str, str, str], list[dict[str, object]]] = OrderedDict()
    passthrough_lines: list[str] = []

    with validated_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = line[0:6].strip()
            if record in {"ATOM", "HETATM"}:
                atom = _parse_atom_line(line)
                residue_key = (
                    str(atom["chain_id"]),
                    str(atom["residue_seq"]),
                    str(atom["insertion_code"]),
                )
                ordered_residues.setdefault(residue_key, []).append(atom)
            elif record in {"MODEL", "ENDMDL", "TER"}:
                passthrough_lines.append(line)

    if not ordered_residues:
        raise PipelineError("No ATOM or HETATM records were found for protonation")

    # Optionally compute environment-dependent pKa values via PROPKA.
    propka_pkas: dict[tuple[str, str, str], float] = {}
    if use_propka:
        propka_pkas = _run_propka(validated_path)
        if propka_pkas:
            logger.info("PROPKA computed pKa values for %d titratable groups", len(propka_pkas))
        else:
            logger.warning("PROPKA returned no pKa values; falling back to reference pKa values")

    _PHOSPHORYLATED = {"SEP", "TPO", "PTR"}

    protonated_lines: list[str] = []
    decision_lines: list[str] = []
    phosphorylated_found: list[str] = []
    serial = 1
    passthrough_iter = iter(passthrough_lines)

    for residue_index, ((chain_id, residue_seq, insertion_code), atom_records) in enumerate(ordered_residues.items()):
        residue_name = str(atom_records[0]["residue_name"]).upper()
        # Look up whether PROPKA computed a pKa for this residue.
        propka_key = (chain_id.strip(), residue_seq.strip(), residue_name)
        pka_override = propka_pkas.get(propka_key)
        decision = _protonation_decision(residue_name, ph, pka_override=pka_override)
        target_name = residue_name
        hydrogen_names: list[str] = []

        if decision is not None:
            target_name, hydrogen_names, rationale = decision
            decision_lines.append(
                f"{residue_name} {chain_id.strip() or '-'} {residue_seq.strip()}{insertion_code.strip()} -> {target_name}: {rationale}"
            )
            logger.info(
                "Protonation decision for %s chain %s residue %s%s: %s",
                residue_name,
                chain_id.strip() or "-",
                residue_seq.strip(),
                insertion_code.strip(),
                target_name,
            )
            if residue_name in _PHOSPHORYLATED:
                phosphorylated_found.append(
                    f"{residue_name} {chain_id.strip()}{residue_seq.strip()}"
                )

        # Identify hydrogens already present so we only add missing ones.
        existing_hydrogens = {str(atom["atom_name"]).upper() for atom in atom_records if str(atom["element"]).upper() == "H"}
        centroid_x, centroid_y, centroid_z = _centroid(atom_records)

        # Re-emit all original atoms with the (possibly renamed) residue.
        for atom in atom_records:
            updated_atom = dict(atom)
            updated_atom["residue_name"] = target_name
            protonated_lines.append(_format_atom_line(serial, updated_atom))
            serial += 1

        # Place any missing protonation-state hydrogens at small offsets
        # from the heavy-atom centroid.  Exact coordinates are refined
        # downstream by the force-field minimizer.
        for hydrogen_index, hydrogen_name in enumerate(hydrogen_names):
            if hydrogen_name.upper() in existing_hydrogens:
                continue
            hydrogen_atom = {
                "record": "ATOM",
                "atom_name": hydrogen_name,
                "residue_name": target_name,
                "chain_id": chain_id,
                "residue_seq": residue_seq,
                "insertion_code": insertion_code,
                "x": centroid_x + 0.1 * (hydrogen_index + 1),
                "y": centroid_y + 0.05 * (hydrogen_index + 1),
                "z": centroid_z + 0.03 * (hydrogen_index + 1),
                "occupancy": "1.00",
                "temp_factor": "20.00",
                "element": "H",
            }
            protonated_lines.append(_format_atom_line(serial, hydrogen_atom))
            serial += 1

        if residue_index < len(ordered_residues) - 1:
            protonated_lines.append("TER\n")

    protonated_lines.append("END\n")

    if phosphorylated_found:
        warnings.warn(
            f"Phosphorylated residues detected ({', '.join(phosphorylated_found)}); "
            "these may require specialized AMBER parameters (e.g., phosaa14SB).",
            UserWarning,
            stacklevel=2,
        )

    output_path.write_text("".join(protonated_lines), encoding="utf-8")
    log_path.write_text("\n".join(decision_lines) + ("\n" if decision_lines else ""), encoding="utf-8")

    logger.info("Protonated PDB written to %s", output_path)
    logger.info("Protonation decision log written to %s", log_path)
    return output_path
