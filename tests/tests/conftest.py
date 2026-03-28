"""Shared pytest fixtures for the SPINK7-KLK5 MD pipeline test suite."""

from __future__ import annotations

import pytest


def pytest_configure(config):
    """Register custom pytest markers."""
    config.addinivalue_line(
        "markers",
        "optimized: marks tests that verify behavior under python -O (assert stripped)",
    )
