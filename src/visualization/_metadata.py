"""Figure provenance metadata utility."""

from __future__ import annotations

import datetime
from pathlib import Path


def figure_metadata(extra: dict[str, str] | None = None) -> dict[str, str]:
    """Build a metadata dict for embedding into saved figures.

    PNG files support arbitrary text metadata via matplotlib's
    ``savefig(metadata={...})`` parameter. SVG supports ``metadata``
    as well. PDF embeds via the ``metadata`` kwarg to the PDF backend.

    Args:
        extra: Optional additional key-value pairs to include.

    Returns:
        dict[str, str]: Metadata dictionary suitable for ``savefig()``.
    """
    meta: dict[str, str] = {
        "Creator": "SPINK7-KLK5 MD Pipeline",
        "Author": "Ryan Kamp",
        "Source": "University of Cincinnati — Dept. of Computer Science",
        "CreationDate": datetime.datetime.now(tz=datetime.timezone.utc).isoformat(),
    }
    if extra:
        meta.update(extra)
    return meta


# SVG and PDF backends accept only specific metadata keys.
_SVG_VALID_KEYS = {
    "Creator", "Date", "Format", "Identifier", "Keywords", "Language",
    "Relation", "Rights", "Source", "Subject", "Title", "Type",
    "Coverage", "Description", "Contributor", "Publisher",
}
_PDF_VALID_KEYS = {
    "Title", "Author", "Subject", "Keywords", "Creator", "Producer",
    "CreationDate", "ModDate", "Trapped",
}


def metadata_for_format(output_path: Path) -> dict[str, str]:
    """Return a metadata dict filtered to keys valid for *output_path*'s format.

    PNG accepts arbitrary keys; SVG and PDF accept only backend-specific subsets.
    """
    meta = figure_metadata()
    ext = Path(output_path).suffix.lower()
    if ext == ".svg":
        return {k: v for k, v in meta.items() if k in _SVG_VALID_KEYS}
    if ext == ".pdf":
        filtered = {k: v for k, v in meta.items() if k in _PDF_VALID_KEYS}
        # PDF backend expects CreationDate as a datetime object.
        if "CreationDate" in filtered:
            filtered["CreationDate"] = datetime.datetime.now(  # type: ignore[assignment]
                tz=datetime.timezone.utc,
            )
        return filtered
    # PNG and EPS accept arbitrary keys.
    return meta
