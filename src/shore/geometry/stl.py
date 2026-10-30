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
from shore.geometry.base import BODY, Display, EntityKind, EntityRef, Projection
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
        # closest point per query point, by the point's exact bytes (ADR 0032): a point's
        # answer does not depend on the other points of its batch
        self._closest: dict[bytes, tuple[np.ndarray, float, int]] = {}
        self._features: tuple[list[np.ndarray], np.ndarray] | None = None

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

    def display(self) -> Display:
        """The mesh as drawn (ADR 0038): every triangle belongs to the body;
        the feature curves and their junctions (ADR 0049) drawn over it."""
        faces = np.asarray(self.mesh.faces)
        curves, points = self.features()
        return Display(
            points=np.asarray(self.mesh.vertices, dtype=np.float64),
            triangles=faces,
            triangle_entity=np.zeros(len(faces), dtype=np.int64),
            surface_kind=EntityKind.BODY,
            curves=tuple(c.copy() for c in curves),
            vertices=points.copy(),
        )

    def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]:
        if kind is EntityKind.BODY:
            return (BODY,)
        curves, points = self.features()
        if kind is EntityKind.CURVE:
            return tuple(EntityRef(kind, i) for i in range(len(curves)))
        if kind is EntityKind.POINT:
            return tuple(EntityRef(kind, i) for i in range(len(points)))
        return ()

    def features(self) -> tuple[list[np.ndarray], np.ndarray]:
        """Feature curves and points (ADR 0049), computed once.

        Feature edges are the mesh edges whose two facets' normals differ by
        more than ``crease_deg``, and the open-boundary edges. They are
        chained into polylines (the CURVE entities) split at every vertex
        where other than two feature edges meet; those vertices are the POINT
        entities. Chains, and the points, are ordered by their lowest vertex
        index, closed loops starting there.
        """
        if self._features is None:
            self._features = _feature_curves(self.mesh, self.crease_deg)
        return self._features

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
        if entity is not None and entity.kind in (EntityKind.CURVE, EntityKind.POINT):
            curves, points_ = self.features()
            if entity.kind is EntityKind.CURVE:
                closest = _closest_on_polyline(curves[entity.index], flat)
            else:
                closest = np.broadcast_to(points_[entity.index], flat.shape).copy()
            return Projection(
                points=closest.reshape(*shape, 3),
                distances=np.linalg.norm(closest - flat, axis=-1).reshape(shape),
                normals=np.full((*shape, 3), np.nan),
            )
        closest, distances, tri_ids = self._closest_points(flat)
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
        """Fingerprint of the body (area-weighted facet centroid, area, bbox),
        a feature curve (length-weighted centroid, length, bbox) or a point."""
        self._check_entity(entity)
        if entity.kind is EntityKind.CURVE:
            curve = self.features()[0][entity.index]
            seg = np.linalg.norm(np.diff(curve, axis=0), axis=1)
            mid = 0.5 * (curve[1:] + curve[:-1])
            centroid = (mid * seg[:, None]).sum(axis=0) / max(float(seg.sum()), 1e-300)
            box = np.array([curve.min(axis=0), curve.max(axis=0)])
            return make_key(
                EntityKind.CURVE, "polyline", centroid, float(seg.sum()), box, self.bounds
            )
        if entity.kind is EntityKind.POINT:
            point = self.features()[1][entity.index]
            return make_key(
                EntityKind.POINT, "Point", point, 0.0, np.array([point, point]), self.bounds
            )
        areas = self.mesh.area_faces
        centroid = (self.mesh.triangles_center * areas[:, None]).sum(axis=0) / areas.sum()
        return make_key(
            EntityKind.BODY, "mesh", centroid, float(areas.sum()), self.bounds, self.bounds
        )

    # -- STL specifics -------------------------------------------------------

    _CACHE_SIZE = 200_000  # points remembered per geometry (about 20 MB)

    def _closest_points(self, flat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """trimesh's closest points of *flat* (n, 3), each point computed once."""
        keys = [p.tobytes() for p in flat]
        todo = [i for i, k in enumerate(keys) if k not in self._closest]
        if todo:
            new = flat[todo]
            closest, distances, tri_ids = trimesh.proximity.closest_point(self.mesh, new)
            if len(self._closest) + len(todo) > self._CACHE_SIZE:
                self._closest.clear()
            for row, i in enumerate(todo):
                self._closest[keys[i]] = (closest[row], float(distances[row]), int(tri_ids[row]))
        hits = [self._closest[k] for k in keys]
        return (
            np.array([h[0] for h in hits], dtype=np.float64).reshape(-1, 3),
            np.array([h[1] for h in hits], dtype=np.float64),
            np.array([h[2] for h in hits], dtype=np.int64),
        )

    def projector(self, anchor: np.ndarray) -> Projector:
        """The cached ``Projector`` for *anchor* (one BVH build per anchor)."""
        a = np.asarray(anchor, dtype=np.float64)
        key = a.tobytes()
        if key not in self._projectors:
            self._projectors[key] = Projector(self.mesh, a)
        return self._projectors[key]

    def _check_entity(self, entity: EntityRef | None) -> None:
        if entity is None or entity == BODY:
            return
        if entity.kind in (EntityKind.CURVE, EntityKind.POINT):
            curves, points = self.features()
            count = len(curves) if entity.kind is EntityKind.CURVE else len(points)
            if 0 <= entity.index < count:
                return
        raise ParameterError(
            f"STL geometry has the body, {len(self.features()[0])} feature curve(s) and "
            f"{len(self.features()[1])} point(s); got {entity}."
        )


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


def _feature_curves(
    mesh: trimesh.Trimesh, crease_deg: float
) -> tuple[list[np.ndarray], np.ndarray]:
    """Feature polylines and junction points of *mesh* (``StlGeometry.features``)."""
    angles = np.asarray(mesh.face_adjacency_angles)
    sharp = np.sort(np.asarray(mesh.face_adjacency_edges)[angles > np.radians(crease_deg)], axis=1)
    edges, counts = np.unique(np.asarray(mesh.edges_sorted), axis=0, return_counts=True)
    feature = np.unique(np.vstack([sharp.reshape(-1, 2), edges[counts == 1]]), axis=0)
    neighbours: dict[int, list[int]] = {}
    for a, b in feature.tolist():
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)
    for v in neighbours:
        neighbours[v].sort()
    junctions = sorted(v for v, n in neighbours.items() if len(n) != 2)
    used: set[tuple[int, int]] = set()
    chains: list[list[int]] = []

    def walk(start: int, nxt: int) -> list[int]:
        chain, prev, cur = [start], start, nxt
        used.add((min(start, nxt), max(start, nxt)))
        while True:
            chain.append(cur)
            if cur == start or len(neighbours[cur]) != 2:
                return chain
            step = next(v for v in neighbours[cur] if v != prev)
            if (min(cur, step), max(cur, step)) in used:
                return chain
            used.add((min(cur, step), max(cur, step)))
            prev, cur = cur, step

    for j in junctions:
        for n in neighbours[j]:
            if (min(j, n), max(j, n)) not in used:
                chains.append(walk(j, n))
    for v in sorted(neighbours):  # closed loops without junctions
        for n in neighbours[v]:
            if (min(v, n), max(v, n)) not in used:
                chains.append(walk(v, n))
    chains.sort(key=lambda c: (min(c), c))
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    return [vertices[c] for c in chains], vertices[junctions].reshape(-1, 3)


def _closest_on_polyline(curve: np.ndarray, points: np.ndarray) -> np.ndarray:
    """The closest points (n, 3) of the polyline *curve* (m, 3) to *points* (n, 3)."""
    a, b = curve[:-1], curve[1:]
    ab = b - a
    length2 = np.maximum((ab * ab).sum(axis=1), 1e-300)
    t = np.clip(np.einsum("nsj,sj->ns", points[:, None] - a[None], ab) / length2, 0.0, 1.0)
    candidates = a[None] + t[..., None] * ab[None]
    best = np.argmin(((candidates - points[:, None]) ** 2).sum(axis=-1), axis=1)
    return candidates[np.arange(len(points)), best]
