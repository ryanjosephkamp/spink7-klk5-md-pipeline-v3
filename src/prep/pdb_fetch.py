"""Structure retrieval utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
from pathlib import Path

import requests

from src import PipelineError


logger = logging.getLogger(__name__)

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"


def _validate_output_dir(output_dir: Path) -> Path:
    """Ensure that the output directory exists and is usable."""

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def _download_to_path(url: str, output_path: Path) -> Path:
    """Download a structure file to disk and return the saved path."""

    try:
        response = requests.get(url, timeout=30)
    except requests.RequestException as exc:
        raise PipelineError(f"Failed to download structure from {url}") from exc

    if response.status_code == 404:
        raise FileNotFoundError(f"Structure not found at {url}")

    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        raise PipelineError(f"Unexpected response while downloading structure from {url}") from exc

    if not response.content:
        raise PipelineError(f"Downloaded structure from {url} is empty")

    output_path.write_bytes(response.content)
    logger.info("Downloaded structure to %s", output_path)
    return output_path


def fetch_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download a PDB structure from RCSB.

    Invariants: None.

    Args:
        pdb_id: Four-character RCSB identifier.
        output_dir: Directory where the downloaded PDB file is written.

    Returns:
        Path: Path to the downloaded PDB file.

    Raises:
        ValueError: If the PDB identifier is invalid.
        FileNotFoundError: If the structure does not exist on RCSB.
        PipelineError: If the download fails for a non-404 reason.
    """

    normalized_id = pdb_id.strip().upper()
    if len(normalized_id) != 4 or not normalized_id.isalnum():
        raise ValueError("pdb_id must be a 4-character alphanumeric identifier")

    output_directory = _validate_output_dir(output_dir)
    output_path = output_directory / f"{normalized_id}.pdb"
    download_url = RCSB_PDB_URL.format(pdb_id=normalized_id)

    logger.info("Fetching PDB structure %s", normalized_id)
    return _download_to_path(download_url, output_path)


def fetch_alphafold(uniprot_id: str, output_dir: Path) -> Path:
    """Download an AlphaFold model using a UniProt accession.

    Invariants: None.

    Args:
        uniprot_id: UniProt accession identifier.
        output_dir: Directory where the downloaded model file is written.

    Returns:
        Path: Path to the downloaded AlphaFold PDB file.

    Raises:
        ValueError: If the UniProt identifier is empty.
        FileNotFoundError: If the AlphaFold model is unavailable.
        PipelineError: If the download fails for a non-404 reason.
    """

    normalized_id = uniprot_id.strip().upper()
    if not normalized_id:
        raise ValueError("uniprot_id must be a non-empty identifier")

    output_directory = _validate_output_dir(output_dir)
    output_path = output_directory / f"AF-{normalized_id}.pdb"
    download_url = ALPHAFOLD_PDB_URL.format(uniprot_id=normalized_id)

    logger.info("Fetching AlphaFold model %s", normalized_id)
    return _download_to_path(download_url, output_path)
