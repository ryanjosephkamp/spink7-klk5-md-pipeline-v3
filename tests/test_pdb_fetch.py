"""Tests for structure retrieval utilities."""

from __future__ import annotations

from pathlib import Path

import pytest
import requests

from src import PipelineError
from src.prep.pdb_fetch import fetch_alphafold, fetch_pdb


class _MockResponse:
    """Minimal mock response for requests-based downloader tests."""

    def __init__(self, status_code: int, content: bytes) -> None:
        self.status_code = status_code
        self.content = content

    def raise_for_status(self) -> None:
        if self.status_code >= 400 and self.status_code != 404:
            raise requests.HTTPError(f"HTTP {self.status_code}")


def test_fetch_pdb_downloads_valid_pdb(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Chunk 4 gate: fetch_pdb should write a valid PDB payload to disk."""

    pdb_payload = b"HEADER    TEST PDB\nATOM      1  N   ALA A   1      11.104  13.207   9.331\nEND\n"

    def mock_get(url: str, timeout: int) -> _MockResponse:
        assert url.endswith("/1ABC.pdb")
        assert timeout == 30
        return _MockResponse(status_code=200, content=pdb_payload)

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)

    output_path = fetch_pdb("1abc", tmp_path / "raw")

    assert output_path == tmp_path / "raw" / "1ABC.pdb"
    assert output_path.exists()
    assert output_path.read_bytes() == pdb_payload


def test_fetch_pdb_handles_404_gracefully(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Chunk 4 gate: a missing PDB should raise FileNotFoundError."""

    def mock_get(url: str, timeout: int) -> _MockResponse:
        assert url.endswith("/9ZZZ.pdb")
        assert timeout == 30
        return _MockResponse(status_code=404, content=b"")

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)

    with pytest.raises(FileNotFoundError, match="Structure not found"):
        fetch_pdb("9ZZZ", tmp_path)


def test_fetch_alphafold_downloads_valid_model(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """AlphaFold retrieval should write the model payload to disk."""

    model_payload = b"HEADER    ALPHAFOLD MODEL\nATOM      1  CA  GLY A   1       0.000   0.000   0.000\nEND\n"

    def mock_get(url: str, timeout: int) -> _MockResponse:
        assert url.endswith("AF-Q6UWN8-F1-model_v4.pdb")
        assert timeout == 30
        return _MockResponse(status_code=200, content=model_payload)

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)

    output_path = fetch_alphafold("Q6UWN8", tmp_path / "alphafold")

    assert output_path == tmp_path / "alphafold" / "AF-Q6UWN8.pdb"
    assert output_path.exists()
    assert output_path.read_bytes() == model_payload


def test_fetch_pdb_rejects_invalid_identifier(tmp_path: Path) -> None:
    """Public API should reject malformed RCSB identifiers before any network call."""

    with pytest.raises(ValueError, match="pdb_id must be a 4-character alphanumeric identifier"):
        fetch_pdb("ABC", tmp_path)


def test_fetch_pdb_raises_pipeline_error_on_transport_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Transport-level request failures should raise PipelineError."""

    def mock_get(url: str, timeout: int) -> _MockResponse:
        raise requests.ConnectionError("network down")

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)

    with pytest.raises(PipelineError, match="Failed to download structure"):
        fetch_pdb("1ABC", tmp_path)
