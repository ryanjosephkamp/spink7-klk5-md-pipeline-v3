"""CLI entry point for alchemical FEP mutagenesis campaigns."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for the FEP mutagenesis workflow."""
    parser = argparse.ArgumentParser(
        description="Run alchemical FEP for computational mutagenesis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--complex-pdb", type=Path, required=True,
        help="Path to the prepared complex PDB file.",
    )
    parser.add_argument(
        "--chain-id", type=str, required=True,
        help="Chain identifier containing the residue to mutate.",
    )
    parser.add_argument(
        "--residue-number", type=int, required=True,
        help="Residue number to mutate.",
    )
    parser.add_argument(
        "--mutation", type=str, required=True,
        help="Three-letter code of the target residue (e.g., GLY).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/analysis/fep"),
        help="Output directory for FEP results.",
    )
    parser.add_argument(
        "--n-lambda-windows", type=int, default=20,
        help="Number of lambda windows.",
    )
    parser.add_argument(
        "--per-window-ns", type=float, default=2.0,
        help="Production duration per lambda window in ns.",
    )
    parser.add_argument(
        "--temperature-k", type=float, default=310.0,
        help="Simulation temperature in Kelvin.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the FEP mutagenesis campaign.

    Args:
        argv: Command-line arguments (defaults to sys.argv[1:]).

    Returns:
        Exit code: 0 on success, 1 on failure.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.complex_pdb.exists():
        logger.error("PDB file not found: %s", args.complex_pdb)
        return 1

    from openmm.app import ForceField, PDBFile, PME, HBonds
    from openmm import unit

    from src.config import FEPConfig
    from src.simulate.fep import run_fep_campaign
    from src.analyze.fep import compute_delta_g_mbar

    config = FEPConfig(
        n_lambda_windows=args.n_lambda_windows,
        per_window_duration_ns=args.per_window_ns,
        temperature_k=args.temperature_k,
    )

    logger.info(
        "Running FEP: %s chain %s residue %d -> %s",
        args.complex_pdb.name, args.chain_id, args.residue_number, args.mutation,
    )
    logger.info("Lambda windows: %d, per-window: %.1f ns", config.n_lambda_windows, config.per_window_duration_ns)

    # Load structure
    pdb = PDBFile(str(args.complex_pdb))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=HBonds,
    )

    # Identify mutant residue atoms
    mutant_atom_indices = []
    for atom in pdb.topology.atoms():
        if (atom.residue.chain.id == args.chain_id
                and atom.residue.id == str(args.residue_number)):
            mutant_atom_indices.append(atom.index)

    if not mutant_atom_indices:
        logger.error(
            "No atoms found for chain %s residue %d",
            args.chain_id, args.residue_number,
        )
        return 1

    logger.info("Mutant residue: %d atoms at indices %s", len(mutant_atom_indices), mutant_atom_indices)

    # Run FEP campaign
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = run_fep_campaign(
        system=system,
        positions=pdb.positions,
        mutant_atom_indices=mutant_atom_indices,
        config=config,
        output_dir=args.output_dir,
    )

    # Analyze with MBAR
    dg_result = compute_delta_g_mbar(
        results["energy_matrix"],
        results["n_samples_per_state"],
        config.temperature_k,
    )

    logger.info(
        "DG_alch = %.2f +/- %.2f kJ/mol (%.2f +/- %.2f kcal/mol)",
        dg_result["delta_g_kj_mol"],
        dg_result["delta_g_std_kj_mol"],
        dg_result["delta_g_kcal_mol"],
        dg_result["delta_g_std_kcal_mol"],
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
