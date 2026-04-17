"""STL reader: load a raw triangulated surface into a trimesh.Trimesh."""

from __future__ import annotations

from pathlib import Path

import trimesh

from shore.exceptions import MeshInputError


def load_stl(path: str | Path) -> trimesh.Trimesh:
    """Load an STL file and return a watertight-checked Trimesh.

    Parameters
    ----------
    path :
        Path to the .stl file (ASCII or binary).

    Returns
    -------
    trimesh.Trimesh
        The loaded surface mesh.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the loaded geometry is not a single mesh (e.g. a scene with
        multiple components).
    """
    path = Path(path)
    if not path.exists():
        raise MeshInputError(f"STL file not found: {path}")

    geometry = trimesh.load(path, force="mesh")

    if not isinstance(geometry, trimesh.Trimesh):
        raise MeshInputError(f"Expected a single mesh, got {type(geometry)}")

    return geometry
