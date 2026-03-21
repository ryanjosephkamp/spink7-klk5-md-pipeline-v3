"""Tests for structure retrieval utilities."""

from __future__ import annotations

from pathlib import Path

import pytest
import requests

from src import PipelineError
from src.prep.pdb_fetch import fetch_alphafold, fetch_cif, fetch_pdb


class _MockResponse:
    """Minimal mock response for requests-based downloader tests."""

    def __init__(self, status_code: int, content: bytes) -> None:
        self.status_code = status_code
        self.content = content

    def raise_for_status(self) -> None:
        if self.status_code >= 400 and self.status_code != 404:
            raise requests.HTTPError(f"HTTP {self.status_code}")


@pytest.fixture(autouse=True)
def _mock_cache_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Redirect the PDB cache to a temp directory for every test."""
    cache_path = tmp_path / "test_cache"
    cache_path.mkdir()
    monkeypatch.setattr("src.prep.pdb_fetch._cache_dir", lambda: cache_path)


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


def test_fetch_cif_rejects_invalid_pdb_id(tmp_path: Path) -> None:
    """fetch_cif should validate the PDB ID format identically to fetch_pdb."""
    with pytest.raises(ValueError, match="4-character alphanumeric"):
        fetch_cif("XY", tmp_path)

    with pytest.raises(ValueError, match="4-character alphanumeric"):
        fetch_cif("toolong", tmp_path)


def test_fetch_cif_downloads_valid_cif(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """fetch_cif should write a valid CIF payload to disk."""
    cif_payload = b"data_1ABC\n_cell.length_a  50.0\n"

    def mock_get(url: str, timeout: int) -> _MockResponse:
        assert url.endswith("/1ABC.cif")
        assert timeout == 30
        return _MockResponse(status_code=200, content=cif_payload)

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)

    output_path = fetch_cif("1abc", tmp_path / "raw")
    assert output_path == tmp_path / "raw" / "1ABC.cif"
    assert output_path.exists()
    assert output_path.read_bytes() == cif_payload


def test_fetch_pdb_raises_pipeline_error_on_transport_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Transport-level request failures should raise PipelineError."""

    def mock_get(url: str, timeout: int) -> _MockResponse:
        raise requests.ConnectionError("network down")

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)
    monkeypatch.setattr("src.prep.pdb_fetch.time.sleep", lambda _: None)

    with pytest.raises(PipelineError, match="Failed to download structure"):
        fetch_pdb("1ABC", tmp_path)


# ---------------------------------------------------------------------------
# L-22 Step 1 — Exponential backoff retry verification tests
# ---------------------------------------------------------------------------


def test_fetch_pdb_retries_on_transient_503(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """fetch_pdb should retry on HTTP 503 and succeed on the third attempt."""
    attempt_count = {"value": 0}
    pdb_payload = b"HEADER    TEST PDB\nEND\n"

    def mock_get(url: str, timeout: int) -> _MockResponse:
        attempt_count["value"] += 1
        if attempt_count["value"] < 3:
            return _MockResponse(status_code=503, content=b"")
        return _MockResponse(status_code=200, content=pdb_payload)

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)
    monkeypatch.setattr("src.prep.pdb_fetch.time.sleep", lambda _: None)

    output_path = fetch_pdb("1ABC", tmp_path / "retry_test")
    assert output_path.read_bytes() == pdb_payload
    assert attempt_count["value"] == 3


def test_fetch_pdb_raises_after_exhausting_retries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """fetch_pdb should raise PipelineError after all retries are exhausted."""

    def mock_get(url: str, timeout: int) -> _MockResponse:
        return _MockResponse(status_code=503, content=b"")

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)
    monkeypatch.setattr("src.prep.pdb_fetch.time.sleep", lambda _: None)

    with pytest.raises(PipelineError, match="after 3 attempts"):
        fetch_pdb("1ABC", tmp_path)


# ---------------------------------------------------------------------------
# L-22 Step 2 — Local disk cache verification test
# ---------------------------------------------------------------------------


def test_fetch_pdb_uses_cache_on_second_call(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Second fetch_pdb call should use the cache — no network request."""
    pdb_payload = b"HEADER    CACHED PDB\nEND\n"
    call_count = {"value": 0}

    def mock_get(url: str, timeout: int) -> _MockResponse:
        call_count["value"] += 1
        return _MockResponse(status_code=200, content=pdb_payload)

    monkeypatch.setattr("src.prep.pdb_fetch.requests.get", mock_get)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    monkeypatch.setattr("src.prep.pdb_fetch._cache_dir", lambda: cache_dir)

    fetch_pdb("1ABC", tmp_path / "run1")
    assert call_count["value"] == 1

    fetch_pdb("1ABC", tmp_path / "run2")
    assert call_count["value"] == 1  # no second network call
    assert (tmp_path / "run2" / "1ABC.pdb").read_bytes() == pdb_payload
