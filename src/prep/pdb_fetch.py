"""Structure retrieval utilities for the SPINK7-KLK5 MD pipeline."""

from __future__ import annotations

import logging
import time
from pathlib import Path

import requests

from src import PipelineError


logger = logging.getLogger(__name__)

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

_MAX_RETRIES = 3
_BACKOFF_BASE_SECONDS = 1.0
_BACKOFF_MULTIPLIER = 4.0
_RETRYABLE_STATUS_CODES = {500, 502, 503, 504}


_CACHE_DIR_NAME = "cache"


def _cache_dir() -> Path:
    """Return the local PDB cache directory, creating it if necessary."""
    from src.config import DATA_DIR

    cache_path = DATA_DIR / "pdb" / _CACHE_DIR_NAME
    cache_path.mkdir(parents=True, exist_ok=True)
    return cache_path


def _validate_output_dir(output_dir: Path) -> Path:
    """Ensure that the output directory exists and is usable."""

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def _download_to_path(url: str, output_path: Path) -> Path:
    """Download a structure file to disk with exponential backoff retry.

    Retries up to ``_MAX_RETRIES`` times for transient network errors and
    server-side 5xx responses.  Permanent errors (404, invalid content)
    are raised immediately without retry.

    Args:
        url: The URL to download from.
        output_path: The local path to write the downloaded content to.

    Returns:
        Path to the downloaded file.

    Raises:
        FileNotFoundError: If the server returns HTTP 404.
        PipelineError: If all retry attempts are exhausted or a
            non-retryable error occurs.
    """
    last_exception: Exception | None = None

    for attempt in range(1, _MAX_RETRIES + 1):
        try:
            response = requests.get(url, timeout=30)
        except requests.RequestException as exc:
            last_exception = exc
            if attempt < _MAX_RETRIES:
                delay = _BACKOFF_BASE_SECONDS * (_BACKOFF_MULTIPLIER ** (attempt - 1))
                logger.warning(
                    "Download attempt %d/%d failed for %s: %s. Retrying in %.0f s.",
                    attempt, _MAX_RETRIES, url, exc, delay,
                )
                time.sleep(delay)
                continue
            raise PipelineError(
                f"Failed to download structure from {url} after {_MAX_RETRIES} attempts"
            ) from exc

        if response.status_code == 404:
            raise FileNotFoundError(f"Structure not found at {url}")

        if response.status_code in _RETRYABLE_STATUS_CODES:
            last_exception = requests.HTTPError(f"HTTP {response.status_code}")
            if attempt < _MAX_RETRIES:
                delay = _BACKOFF_BASE_SECONDS * (_BACKOFF_MULTIPLIER ** (attempt - 1))
                logger.warning(
                    "Download attempt %d/%d got HTTP %d for %s. Retrying in %.0f s.",
                    attempt, _MAX_RETRIES, response.status_code, url, delay,
                )
                time.sleep(delay)
                continue
            raise PipelineError(
                f"Failed to download structure from {url} after {_MAX_RETRIES} attempts"
            ) from last_exception

        try:
            response.raise_for_status()
        except requests.HTTPError as exc:
            raise PipelineError(
                f"Unexpected response while downloading structure from {url}"
            ) from exc

        if not response.content:
            raise PipelineError(f"Downloaded structure from {url} is empty")

        output_path.write_bytes(response.content)
        logger.info("Downloaded structure to %s (attempt %d/%d)", output_path, attempt, _MAX_RETRIES)
        return output_path

    raise PipelineError(
        f"Failed to download structure from {url} after {_MAX_RETRIES} attempts"
    )


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

    cached_path = _cache_dir() / f"{normalized_id}.pdb"
    output_directory = _validate_output_dir(output_dir)
    output_path = output_directory / f"{normalized_id}.pdb"

    if cached_path.exists() and cached_path.stat().st_size > 0:
        logger.info("Cache hit for PDB %s at %s", normalized_id, cached_path)
        output_path.write_bytes(cached_path.read_bytes())
        return output_path

    download_url = RCSB_PDB_URL.format(pdb_id=normalized_id)
    logger.info("Cache miss for PDB %s; fetching from RCSB", normalized_id)
    result = _download_to_path(download_url, output_path)

    cached_path.write_bytes(output_path.read_bytes())
    logger.info("Cached PDB %s at %s", normalized_id, cached_path)
    return result


def fetch_cif(pdb_id: str, output_dir: Path) -> Path:
    """Download an mmCIF structure from RCSB.

    Invariants: None.

    Args:
        pdb_id: Four-character RCSB identifier.
        output_dir: Directory where the downloaded CIF file is written.

    Returns:
        Path to the downloaded CIF file.

    Raises:
        ValueError: If the PDB identifier is invalid.
        FileNotFoundError: If the structure does not exist on RCSB.
        PipelineError: If the download fails for a non-404 reason.
    """
    normalized_id = pdb_id.strip().upper()
    if len(normalized_id) != 4 or not normalized_id.isalnum():
        raise ValueError("pdb_id must be a 4-character alphanumeric identifier")

    output_directory = _validate_output_dir(output_dir)
    output_path = output_directory / f"{normalized_id}.cif"
    download_url = RCSB_CIF_URL.format(pdb_id=normalized_id)

    logger.info("Fetching CIF structure %s", normalized_id)
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

    cached_path = _cache_dir() / f"AF-{normalized_id}.pdb"
    output_directory = _validate_output_dir(output_dir)
    output_path = output_directory / f"AF-{normalized_id}.pdb"

    if cached_path.exists() and cached_path.stat().st_size > 0:
        logger.info("Cache hit for AlphaFold %s at %s", normalized_id, cached_path)
        output_path.write_bytes(cached_path.read_bytes())
        return output_path

    download_url = ALPHAFOLD_PDB_URL.format(uniprot_id=normalized_id)
    logger.info("Cache miss for AlphaFold %s; fetching from EBI", normalized_id)
    result = _download_to_path(download_url, output_path)

    cached_path.write_bytes(output_path.read_bytes())
    logger.info("Cached AlphaFold %s at %s", normalized_id, cached_path)
    return result
