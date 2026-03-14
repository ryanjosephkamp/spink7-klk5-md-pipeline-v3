"""Core package for the SPINK7-KLK5 molecular dynamics pipeline."""


class PipelineError(Exception):
    """Base exception for all pipeline errors."""


class PhysicalValidityError(PipelineError):
    """Raised when a physical invariant is violated."""


class ConvergenceError(PipelineError):
    """Raised when an iterative solver fails to converge."""


class InsufficientSamplingError(PipelineError):
    """Raised when sampling is insufficient for reliable estimates."""
