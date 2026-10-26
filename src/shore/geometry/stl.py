"""STL backend of the ``Geometry`` protocol (ADR 0002).

Wraps a ``trimesh.Trimesh``. The body is the only entity; there is no
parametrisation. Ray casts go through the existing ``Projector`` (one cached
BVH per anchor), so points, distances and ray casts are bit-identical to the
pre-protocol code.

Normals are crease-aware (ADR 0014): each triangle corner gets the
corner-angle-weighted mean of the normals of the faces around its vertex
that lie within ``crease_deg`` of the triangle's own, and the corners of an
open boundary that lies in a plane (a half model's cut) get their normal
projected into that plane. They differ from the legacy
``volume.hyperbolic._stl_normals`` wherever faces meet at a sharp edge.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import trimesh
import trimesh.graph
import trimesh.grouping
import trimesh.proximity
import trimesh.triangles

from shore.exceptions import ParameterError
from shore.geometry.base import BODY, EntityKind, EntityRef, Projection
from shore.geometry.keys import EntityKey, make_key
from shore.surface.io import load_stl
from shore.surface.projector import ProjectionResult, Projector

_CHUNK = 65536  # faces per block of the corner-normal accumulation


class StlGeometry:
    """Triangulated body surface as a ``Geometry``.

    Parameters
    ----------
    mesh :
        Triangulated body surface; ``load_stl`` (used by ``from_file``) cleans
        it so that face orientation and normals are consistent.
    crease_deg :
        Faces whose normals differ by more than this angle are not averaged
        into each other's corner normals (sharp edges). 180 on a closed mesh
        gives trimesh's vertex normals.

    Raises
    ------
    ParameterError
        *crease_deg* outside (0, 180].
    """

    def __init__(self, mesh: trimesh.Trimesh, crease_deg: float = 30.0) -> None:
        crease = float(crease_deg)
        if not 0.0 < crease <= 180.0:  # also rejects NaN
            raise ParameterError(f"crease_deg must be in (0, 180], got {crease_deg!r}.")
        self.mesh = mesh
        self.crease_deg = crease
        self._projectors: dict[bytes, Projector] = {}
        self._corner_normals = _corner_normals(mesh, crease)

    @classmethod
    def from_file(
        cls, path: str | Path, skip_hygiene: bool = False, crease_deg: float = 30.0
    ) -> StlGeometry:
        """Load an STL file (see ``shore.surface.io.load_stl``)."""
        return cls(load_stl(path, skip_hygiene=skip_hygiene), crease_deg=crease_deg)

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
        crease-aware corner normals of the closest triangle (ADR 0014); on a
        crease edge the closest triangle, and so the normal, is either
        side's.
        """
        self._check_entity(entity)
        pts = np.asarray(points, dtype=np.float64)
        shape = pts.shape[:-1]
        flat = pts.reshape(-1, 3)
        closest, distances, tri_ids = trimesh.proximity.closest_point(self.mesh, flat)
        bary = trimesh.triangles.points_to_barycentric(self.mesh.triangles[tri_ids], closest)
        normals = (bary[..., None] * self._corner_normals[tri_ids]).sum(axis=1)
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


def _corner_normals(mesh: trimesh.Trimesh, crease_deg: float) -> np.ndarray:
    """Crease-aware normal of every triangle corner, shape (n_faces, 3, 3).

    Thuerrner & Wuethrich corner-angle weighting (trimesh's vertex normals)
    over the faces around the vertex within *crease_deg* of the corner's
    triangle; then planar open boundaries (see ``_planar_boundaries``).
    Accumulated one vertex-face slot at a time over chunks of faces, so the
    temporaries stay bounded.
    """
    faces = np.asarray(mesh.faces)
    face_normals = np.asarray(mesh.face_normals, dtype=np.float64)
    face_angles = np.asarray(mesh.face_angles, dtype=np.float64)
    vertex_faces = np.asarray(mesh.vertex_faces)  # (n_vertices, k), -1 padded
    cos_crease = np.cos(np.radians(crease_deg)) - 1e-12
    corner = np.empty((len(faces), 3, 3))
    for lo in range(0, len(faces), _CHUNK):
        f = faces[lo : lo + _CHUNK]
        own = face_normals[lo : lo + _CHUNK, None, :]
        summed = np.zeros((len(f), 3, 3))
        for slot in range(vertex_faces.shape[1]):
            other = vertex_faces[f, slot]  # (chunk, 3): a face around each corner's vertex
            valid = other >= 0
            g = np.where(valid, other, 0)
            at = np.argmax(faces[g] == f[:, :, None], axis=-1)  # its corner at that vertex
            n_other = face_normals[g]
            smooth = (n_other * own).sum(axis=-1) >= cos_crease
            weight = np.where(valid & smooth, face_angles[g, at], 0.0)
            summed += weight[..., None] * n_other
        norm = np.linalg.norm(summed, axis=-1, keepdims=True)
        corner[lo : lo + _CHUNK] = np.where(
            norm > 0.0, summed / np.where(norm > 0.0, norm, 1.0), own
        )
    for vertices, normal in _planar_boundaries(mesh):
        on = np.isin(faces, vertices)
        n = corner[on] - (corner[on] @ normal)[:, None] * normal
        size = np.linalg.norm(n, axis=-1, keepdims=True)
        corner[on] = np.where(size > 1e-12, n / np.where(size > 1e-12, size, 1.0), corner[on])
    return corner


def _planar_boundaries(mesh: trimesh.Trimesh) -> list[tuple[np.ndarray, np.ndarray]]:
    """(vertex indices, unit plane normal) of every open boundary that lies in a plane."""
    edges = np.asarray(mesh.edges_sorted)
    open_edges = edges[trimesh.grouping.group_rows(edges, require_count=1)]
    if not len(open_edges):
        return []
    scale = float(np.linalg.norm(np.ptp(np.asarray(mesh.bounds), axis=0)))
    tol = 1e-9 * scale
    out = []
    for component in trimesh.graph.connected_components(open_edges):
        vertices = np.asarray(sorted(component))
        points = np.asarray(mesh.vertices[vertices], dtype=np.float64)
        centred = points - points.mean(axis=0)
        _, sigma, frame = np.linalg.svd(centred, full_matrices=False)
        if len(sigma) < 3 or sigma[1] <= tol:
            continue  # a point or a line: no plane
        normal = frame[2]
        if np.abs(centred @ normal).max() > tol:
            continue  # not planar
        out.append((vertices, normal))
    return out
