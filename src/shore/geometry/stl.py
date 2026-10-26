"""STL backend of the ``Geometry`` protocol (ADR 0002).

Wraps a ``trimesh.Trimesh``. The body is the only entity; there is no
parametrisation. Ray casts go through the existing ``Projector`` (one cached
BVH per anchor), so results are bit-identical to the pre-protocol code.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import trimesh
import trimesh.proximity
import trimesh.triangles

from shore.exceptions import ParameterError
from shore.geometry.base import BODY, EntityKind, EntityRef, Projection
from shore.geometry.keys import EntityKey, make_key
from shore.surface.io import load_stl
from shore.surface.projector import ProjectionResult, Projector


class StlGeometry:
    """Triangulated body surface as a ``Geometry``.

    Parameters
    ----------
    mesh :
        Triangulated body surface. Normals and vertex normals are taken from
        it as-is; ``load_stl`` (used by ``from_file``) cleans the mesh so that
        vertex normals are meaningful.
    """

    def __init__(self, mesh: trimesh.Trimesh) -> None:
        self.mesh = mesh
        self._projectors: dict[bytes, Projector] = {}

    @classmethod
    def from_file(cls, path: str | Path, skip_hygiene: bool = False) -> StlGeometry:
        """Load an STL file (see ``shore.surface.io.load_stl``)."""
        return cls(load_stl(path, skip_hygiene=skip_hygiene))

    # -- Geometry protocol ---------------------------------------------------

    @property
    def bounds(self) -> np.ndarray:
        return np.asarray(self.mesh.bounds, dtype=np.float64).copy()

    def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]:
        return (BODY,) if kind is EntityKind.BODY else ()

    def ray_cast(
        self, anchor: np.ndarray, dirs: np.ndarray, on_failure: str = "raise"
    ) -> ProjectionResult:
        return self.projector(anchor).project(np.asarray(dirs, dtype=np.float64), on_failure)

    def project(
        self,
        points: np.ndarray,
        entity: EntityRef | None = None,
        *,
        exact: bool = False,
    ) -> Projection:
        """Closest point on the triangulation.

        The facets *are* the geometry of an STL, so the answer is exact for
        either value of *exact*. Normals are barycentric blends of the
        vertex normals of the closest triangle, as in
        ``volume.hyperbolic._stl_normals(interpolate=True)``.
        """
        self._check_entity(entity)
        pts = np.asarray(points, dtype=np.float64)
        shape = pts.shape[:-1]
        flat = pts.reshape(-1, 3)
        closest, distances, tri_ids = trimesh.proximity.closest_point(self.mesh, flat)
        bary = trimesh.triangles.points_to_barycentric(self.mesh.triangles[tri_ids], closest)
        vn = self.mesh.vertex_normals[self.mesh.faces[tri_ids]]
        normals = (bary[..., None] * vn).sum(axis=1)
        mag = np.linalg.norm(normals, axis=-1, keepdims=True)
        normals = normals / np.where(mag < 1e-15, 1.0, mag)
        return Projection(
            points=closest.reshape(*shape, 3),
            distances=np.asarray(distances).reshape(shape),
            normals=normals.reshape(*shape, 3),
        )

    def evaluate(self, entity: EntityRef, params: np.ndarray) -> np.ndarray:
        raise ParameterError("STL geometry has no parametric entities; use project().")

    def fingerprint(self, entity: EntityRef) -> EntityKey:
        """Fingerprint of the body: area-weighted facet centroid, area, bbox."""
        self._check_entity(entity)
        areas = self.mesh.area_faces
        centroid = (self.mesh.triangles_center * areas[:, None]).sum(axis=0) / areas.sum()
        return make_key(
            EntityKind.BODY, "mesh", centroid, float(areas.sum()), self.bounds, self.bounds
        )

    # -- STL specifics -------------------------------------------------------

    def projector(self, anchor: np.ndarray) -> Projector:
        """The cached ``Projector`` for *anchor* (one BVH build per anchor)."""
        a = np.asarray(anchor, dtype=np.float64)
        key = a.tobytes()
        if key not in self._projectors:
            self._projectors[key] = Projector(self.mesh, a)
        return self._projectors[key]

    def _check_entity(self, entity: EntityRef | None) -> None:
        if entity is not None and entity != BODY:
            raise ParameterError(f"STL geometry has only the body entity, got {entity}.")
