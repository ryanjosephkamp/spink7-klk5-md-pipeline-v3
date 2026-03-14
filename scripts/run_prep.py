"""CLI entry point for system preparation."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from openmm import XmlSerializer, unit
from openmm.app import ForceField, HBonds, PME, PDBFile

from src.config import DATA_DIR, SystemConfig
from src.prep.pdb_clean import clean_structure
from src.prep.pdb_fetch import fetch_alphafold, fetch_pdb
from src.prep.protonate import assign_protonation
from src.prep.solvate import solvate_system
from src.prep.topology import build_topology


logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for the preparation workflow."""

    defaults = SystemConfig()
    parser = argparse.ArgumentParser(
        description="Fetch or load a structure, clean it, protonate it, solvate it, and serialize the solvated system.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--input-pdb", type=Path, help="Existing local PDB file to prepare.")
    source_group.add_argument("--pdb-id", help="Four-character RCSB PDB identifier to download.")
    source_group.add_argument("--alphafold-id", help="UniProt accession to retrieve from AlphaFold.")
    parser.add_argument("--chains", nargs="+", required=True, help="Chain identifiers to retain during cleaning.")
    parser.add_argument("--output-root", type=Path, default=DATA_DIR, help="Root data directory for preparation outputs.")
    parser.add_argument("--ph", type=float, default=defaults.ph, help="Target protonation pH.")
    parser.add_argument("--box-padding-nm", type=float, default=defaults.box_padding_nm, help="Solvent box padding in nm.")
    parser.add_argument(
        "--ionic-strength-molar",
        type=float,
        default=defaults.ionic_strength_molar,
        help="Bulk ionic strength in molar units.",
    )
    parser.add_argument("--positive-ion", default=defaults.positive_ion, help="Positive ion name passed to OpenMM.")
    parser.add_argument("--negative-ion", default=defaults.negative_ion, help="Negative ion name passed to OpenMM.")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Python logging level.",
    )
    return parser


def _resolve_input_structure(args: argparse.Namespace, raw_dir: Path) -> Path:
    """Fetch or accept a local structure path for preparation."""

    if args.input_pdb is not None:
        return Path(args.input_pdb)
    if args.pdb_id is not None:
        return fetch_pdb(args.pdb_id, raw_dir)
    return fetch_alphafold(args.alphafold_id, raw_dir)


def main(argv: list[str] | None = None) -> int:
    """Run the preparation CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(name)s: %(message)s")

    system_config = SystemConfig(
        ph=args.ph,
        box_padding_nm=args.box_padding_nm,
        ionic_strength_molar=args.ionic_strength_molar,
        positive_ion=args.positive_ion,
        negative_ion=args.negative_ion,
    )

    output_root = Path(args.output_root)
    raw_dir = output_root / "pdb" / "raw"
    prepared_dir = output_root / "pdb" / "prepared"
    topology_dir = output_root / "topologies"
    raw_dir.mkdir(parents=True, exist_ok=True)
    prepared_dir.mkdir(parents=True, exist_ok=True)
    topology_dir.mkdir(parents=True, exist_ok=True)

    input_path = _resolve_input_structure(args, raw_dir)
    cleaned_path = clean_structure(input_path, chains_to_keep=args.chains)
    protonated_path = assign_protonation(cleaned_path, ph=system_config.ph, force_field="AMBER")
    _, _, modeller = build_topology(protonated_path, system_config)
    solvated_modeller, n_water, n_positive_ions, n_negative_ions = solvate_system(modeller, system_config)

    force_field = ForceField(system_config.force_field, system_config.water_model)
    solvated_system = force_field.createSystem(
        solvated_modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )

    solvated_pdb_path = prepared_dir / f"{protonated_path.stem}_solvated.pdb"
    with solvated_pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(solvated_modeller.topology, solvated_modeller.positions, handle, keepIds=True)

    system_xml_path = topology_dir / f"{protonated_path.stem}_system.xml"
    system_xml_path.write_text(XmlSerializer.serialize(solvated_system), encoding="utf-8")

    logger.info("Prepared solvated structure written to %s", solvated_pdb_path)
    logger.info("Serialized OpenMM system written to %s", system_xml_path)
    logger.info(
        "Solvation summary: %d waters, %d positive ions, %d negative ions",
        n_water,
        n_positive_ions,
        n_negative_ions,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())