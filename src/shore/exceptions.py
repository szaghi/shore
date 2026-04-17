"""SHORE-specific exception hierarchy."""

from __future__ import annotations


class ShoreError(Exception):
    """Base exception for all SHORE-specific errors."""


class MeshInputError(ShoreError, ValueError):
    """Raised when the input geometry is invalid or cannot be read."""


class ParameterError(ShoreError, ValueError):
    """Raised when a function argument is out of range or mutually exclusive."""


class MarchingError(ShoreError, RuntimeError):
    """Raised when hyperbolic marching detects an invalid grid state."""
