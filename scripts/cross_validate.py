"""Cross-validation of Jarzynski and WHAM binding free energy estimates.

Compares ΔG_bind from the Jarzynski equality (SMD) against the WHAM
potential of mean force (umbrella sampling) and enforces the internal
consistency criterion:

    |ΔG_Jarzynski − ΔG_WHAM| < 2 k_B T

At T = 310 K this threshold is ~5.15 kJ/mol (~1.23 kcal/mol).
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np

from src.config import BOLTZMANN_KJ, KCAL_TO_KJ


logger = logging.getLogger(__name__)


def cross_validate(
    jarzynski_npz_path: Path,
    wham_npz_path: Path,
    temperature_k: float = 310.0,
    *,
    threshold_kbt: float = 2.0,
) -> dict[str, float]:
    """Compare Jarzynski and WHAM ΔG estimates.

    Parameters
    ----------
    jarzynski_npz_path : Path
        NPZ file produced by ``scripts/run_analysis.py jarzynski``.
        Must contain key ``delta_g_kj_mol``.
    wham_npz_path : Path
        NPZ file produced by ``scripts/run_analysis.py wham``.
        Must contain keys ``pmf_kj_mol`` and ``xi_bins``.
    temperature_k : float
        Simulation temperature in Kelvin.
    threshold_kbt : float
        Discrepancy threshold in units of k_B T (default 2.0).

    Returns
    -------
    dict[str, float]
        ``delta_g_jarzynski_kj_mol``, ``delta_g_wham_kj_mol``,
        ``discrepancy_kj_mol``, ``discrepancy_kcal_mol``,
        ``threshold_kj_mol``, ``passed``.
    """
    kbt = BOLTZMANN_KJ * temperature_k  # kJ/mol
    threshold_kj = threshold_kbt * kbt

    # Jarzynski estimate.
    jarz = np.load(jarzynski_npz_path)
    dg_jarz = float(jarz["delta_g_kj_mol"])

    # WHAM PMF: ΔG = PMF(bound) − PMF(dissociated).
    wham = np.load(wham_npz_path)
    pmf_kj = np.asarray(wham["pmf_kj_mol"], dtype=float)
    dg_wham = float(np.min(pmf_kj) - pmf_kj[-1])

    discrepancy = abs(dg_jarz - dg_wham)
    passed = discrepancy < threshold_kj

    result = {
        "delta_g_jarzynski_kj_mol": dg_jarz,
        "delta_g_jarzynski_kcal_mol": dg_jarz / KCAL_TO_KJ,
        "delta_g_wham_kj_mol": dg_wham,
        "delta_g_wham_kcal_mol": dg_wham / KCAL_TO_KJ,
        "discrepancy_kj_mol": discrepancy,
        "discrepancy_kcal_mol": discrepancy / KCAL_TO_KJ,
        "threshold_kj_mol": threshold_kj,
        "passed": float(passed),
    }

    if passed:
        logger.info(
            "Cross-validation PASSED: |ΔG_Jarzynski − ΔG_WHAM| = %.2f kJ/mol < %.2f kJ/mol (2 k_BT)",
            discrepancy,
            threshold_kj,
        )
    else:
        logger.warning(
            "Cross-validation FAILED: |ΔG_Jarzynski − ΔG_WHAM| = %.2f kJ/mol ≥ %.2f kJ/mol (2 k_BT)",
            discrepancy,
            threshold_kj,
        )

    return result


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser."""
    parser = argparse.ArgumentParser(
        description="Cross-validate Jarzynski and WHAM ΔG estimates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--jarzynski-npz",
        type=Path,
        required=True,
        help="Path to the Jarzynski summary NPZ file.",
    )
    parser.add_argument(
        "--wham-npz",
        type=Path,
        required=True,
        help="Path to the WHAM PMF NPZ file.",
    )
    parser.add_argument(
        "--temperature-k",
        type=float,
        default=310.0,
        help="Simulation temperature in Kelvin.",
    )
    parser.add_argument(
        "--threshold-kbt",
        type=float,
        default=2.0,
        help="Discrepancy threshold in units of k_B T.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output NPZ path to save cross-validation results.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Python logging level.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the cross-validation CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s %(name)s: %(message)s",
    )

    result = cross_validate(
        args.jarzynski_npz,
        args.wham_npz,
        temperature_k=args.temperature_k,
        threshold_kbt=args.threshold_kbt,
    )

    print(f"Jarzynski ΔG: {result['delta_g_jarzynski_kj_mol']:.2f} kJ/mol "
          f"({result['delta_g_jarzynski_kcal_mol']:.2f} kcal/mol)")
    print(f"WHAM ΔG:      {result['delta_g_wham_kj_mol']:.2f} kJ/mol "
          f"({result['delta_g_wham_kcal_mol']:.2f} kcal/mol)")
    print(f"Discrepancy:  {result['discrepancy_kj_mol']:.2f} kJ/mol "
          f"({result['discrepancy_kcal_mol']:.2f} kcal/mol)")
    print(f"Threshold:    {result['threshold_kj_mol']:.2f} kJ/mol (2 k_BT)")
    print(f"Result:       {'PASSED' if result['passed'] else 'FAILED'}")

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        np.savez(args.output, **result)
        print(f"Results saved to {args.output}")

    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
