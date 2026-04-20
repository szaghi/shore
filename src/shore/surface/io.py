"""STL reader: load a raw triangulated surface into a trimesh.Trimesh."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import trimesh

from shore.exceptions import MeshInputError

log = logging.getLogger(__name__)


def load_stl(path: str | Path, skip_hygiene: bool = False) -> trimesh.Trimesh:
    """Load an STL file and return a trimesh.Trimesh.

    Parameters
    ----------
    path :
        Path to the .stl file (ASCII or binary).
    skip_hygiene :
        When *True*, bypass mesh hygiene steps (merge_vertices,
        remove_degenerate_faces, fix_normals).  Preserves legacy behaviour
        for regression tests.  Default *False*.

    Returns
    -------
    trimesh.Trimesh
        The loaded (and optionally repaired) surface mesh.

    Raises
    ------
    MeshInputError
        If *path* does not exist or the loaded geometry is not a single mesh.
    """
    path = Path(path)
    if not path.exists():
        raise MeshInputError(f"STL file not found: {path}")

    geometry = trimesh.load(path, force="mesh")

    if not isinstance(geometry, trimesh.Trimesh):
        raise MeshInputError(f"Expected a single mesh, got {type(geometry)}")

    if not skip_hygiene:
        n_before = len(geometry.faces)

        geometry.merge_vertices()
        # Remove degenerate (zero-area) faces via Trimesh instance method
        geometry.update_faces(geometry.nondegenerate_faces())
        trimesh.repair.fix_normals(geometry)

        n_after = len(geometry.faces)
        if n_after != n_before:
            log.info(
                "load_stl hygiene: repaired %d triangles (%d → %d)",
                n_before - n_after,
                n_before,
                n_after,
            )

        if not geometry.is_watertight:
            warnings.warn(
                f"Loaded mesh '{path.name}' is not watertight "
                f"({len(geometry.faces)} triangles after hygiene).  "
                "Ray casting may miss some rays or produce incorrect results.  "
                "Consider repairing the STL before use.",
                stacklevel=2,
            )

    return geometry
