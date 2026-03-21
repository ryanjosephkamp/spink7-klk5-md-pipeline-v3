"""OpenMM platform auto-detection and selection.

Provides a centralized platform selection utility that probes
available platforms in descending performance order (CUDA -> OpenCL -> CPU)
and returns the first available.  An explicit platform name can be passed
to override auto-detection.
"""

from __future__ import annotations

import logging

import openmm

logger = logging.getLogger(__name__)

_PLATFORM_PRIORITY: tuple[str, ...] = ("CUDA", "OpenCL", "CPU")


def select_platform(requested: str | None = None) -> openmm.Platform:
    """Select the fastest available OpenMM platform.

    Parameters
    ----------
    requested : str or None
        If provided, return the platform with this exact name.
        Raises ``RuntimeError`` if the platform is unavailable.
        If ``None``, auto-detect by probing in priority order.

    Returns
    -------
    openmm.Platform
        The selected platform instance.

    Raises
    ------
    RuntimeError
        If the requested platform is not available, or no platform
        could be found during auto-detection.
    """
    if requested is not None:
        try:
            platform = openmm.Platform.getPlatformByName(requested)
        except Exception as exc:
            raise RuntimeError(
                f"Requested OpenMM platform '{requested}' is not available"
            ) from exc
        logger.info("Using explicitly requested platform: %s", platform.getName())
        return platform

    for name in _PLATFORM_PRIORITY:
        try:
            platform = openmm.Platform.getPlatformByName(name)
            logger.info("Auto-selected platform: %s", platform.getName())
            return platform
        except Exception:
            logger.debug("Platform '%s' not available, trying next", name)
            continue

    raise RuntimeError(
        "No suitable OpenMM platform found. Tried: "
        + ", ".join(_PLATFORM_PRIORITY)
    )


def platform_performance_summary(platform: openmm.Platform) -> dict[str, str]:
    """Return a dict of platform property names and values for logging.

    Parameters
    ----------
    platform : openmm.Platform
        The active OpenMM platform.

    Returns
    -------
    dict[str, str]
        A mapping of property names to their default values.
    """
    properties: dict[str, str] = {}
    for prop_name in platform.getPropertyNames():
        properties[prop_name] = platform.getPropertyDefaultValue(prop_name)
    return properties
