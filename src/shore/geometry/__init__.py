"""Geometry backends behind one protocol (ADR 0002)."""

from __future__ import annotations

from pathlib import Path

from shore.exceptions import MeshInputError
from shore.geometry.base import BODY, EntityKind, EntityRef, Geometry, Projection
from shore.geometry.stl import StlGeometry

CAD_SUFFIXES = (".step", ".stp", ".iges", ".igs", ".brep", ".brp")

__all__ = [
    "BODY",
    "CAD_SUFFIXES",
    "EntityKind",
    "EntityRef",
    "Geometry",
    "Projection",
    "StlGeometry",
    "load_geometry",
]


def load_geometry(path: str | Path, skip_hygiene: bool = False) -> Geometry:
    """Load a geometry file, choosing the backend from its suffix.

    Parameters
    ----------
    path :
        ``.stl`` (STL backend) or a CAD file (``.step``/``.stp``/``.iges``/
        ``.igs``/``.brep``/``.brp``, OCCT backend).
    skip_hygiene :
        STL only: bypass mesh cleaning (see ``shore.surface.io.load_stl``).

    Raises
    ------
    MeshInputError
        Unknown suffix, unreadable file, or a CAD file while the OCCT backend
        (``shore-mesh[cad]``) is not installed.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".stl":
        return StlGeometry.from_file(path, skip_hygiene=skip_hygiene)
    if suffix in CAD_SUFFIXES:
        try:
            from shore.geometry.occt import OcctGeometry
        except ImportError as err:
            raise MeshInputError(
                f"{path.name}: CAD input needs the OCCT backend: pip install 'shore-mesh[cad]'."
            ) from err
        return OcctGeometry.from_file(path)
    raise MeshInputError(f"{path.name}: unknown geometry format '{suffix}'.")
