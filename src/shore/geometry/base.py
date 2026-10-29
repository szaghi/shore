"""Geometry protocol: the only geometry API the mesh core sees (ADR 0002).

A ``Geometry`` answers batched numpy queries about a body: anchored ray casts
(the existing templates), closest-point projection onto the body or onto one
of its entities, and parametric evaluation where the backend has parameters.
Backends: ``StlGeometry`` (triangulated STL, via trimesh) and, later,
``OcctGeometry`` (STEP/IGES/BREP, via build123d/OCP).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

from shore.surface.projector import ProjectionResult

if TYPE_CHECKING:
    from shore.geometry.keys import EntityKey


class EntityKind(Enum):
    """Topological kind of a geometry entity."""

    BODY = "body"
    SURFACE = "surface"
    CURVE = "curve"
    POINT = "point"


@dataclass(frozen=True)
class EntityRef:
    """Backend-local handle to one geometry entity.

    ``index`` is only meaningful within the backend instance that issued the
    reference; persistent identity across re-imports is the job of the
    association layer (ADR 0003), not of this handle.
    """

    kind: EntityKind
    index: int = 0


BODY = EntityRef(EntityKind.BODY, 0)
"""The whole body, the target of every pre-association projection."""


@dataclass(frozen=True, eq=False)
class Display:
    """A geometry as it is drawn (ADR 0038): triangles, curves and points.

    Each triangle belongs to entity ``EntityRef(surface_kind,
    triangle_entity[t])``: the body of an STL, a face of a CAD model. Curve
    ``c`` is the polyline of ``EntityRef(CURVE, c)``, point ``p`` the
    position of ``EntityRef(POINT, p)``. Coordinates are float64.
    """

    points: np.ndarray  # (n, 3) triangle vertices
    triangles: np.ndarray  # (m, 3) indices into points
    triangle_entity: np.ndarray  # (m,) entity index of each triangle
    surface_kind: EntityKind
    curves: tuple[np.ndarray, ...] = ()
    vertices: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))


@dataclass
class Projection:
    """Closest-point projection of a batch of points.

    Parameters
    ----------
    points : np.ndarray, shape (..., 3)
        Projected positions.
    distances : np.ndarray, shape (...)
        Distance from each query point to its projection.
    normals : np.ndarray, shape (..., 3)
        Unit outward surface normals at the projected points (for a curve
        entity, NaN).
    params : np.ndarray or None
        Parameters of the projected points on the target entity, shape
        (..., 2) for a surface and (...) for a curve; *None* when the
        backend has no parametrisation (STL), and for body projections
        (parameters only make sense per entity).
    entity_index : np.ndarray or None
        For body projections on a multi-entity backend: index of the SURFACE
        entity each point landed on, shape (...); *None* otherwise.
    """

    points: np.ndarray
    distances: np.ndarray
    normals: np.ndarray
    params: np.ndarray | None = None
    entity_index: np.ndarray | None = None


@runtime_checkable
class Geometry(Protocol):
    """Batched geometry queries used by the mesh core.

    All array arguments and results are numpy arrays; leading dimensions of
    query arrays are preserved in the results.
    """

    @property
    def bounds(self) -> np.ndarray:
        """Axis-aligned bounding box, shape (2, 3): ``[min, max]``."""
        ...

    def entities(self, kind: EntityKind) -> tuple[EntityRef, ...]:
        """All entities of *kind* (empty when the backend has none)."""
        ...

    def ray_cast(
        self, anchor: np.ndarray, dirs: np.ndarray, on_failure: str = "raise"
    ) -> ProjectionResult:
        """Cast rays from one interior *anchor* along unit *dirs* (..., 3).

        Semantics and per-ray status follow ``Projector.project``; this is the
        projection the star-shaped templates (cubed-sphere, O-grid) use.
        """
        ...

    def project(
        self,
        points: np.ndarray,
        entity: EntityRef | None = None,
        *,
        exact: bool = False,
    ) -> Projection:
        """Closest-point projection of *points* (..., 3) onto *entity*.

        ``entity=None`` projects onto the whole body. ``exact=False`` allows a
        backend to answer from a tessellated proxy (ADR 0002 tier 1);
        ``exact=True`` requires the exact geometry (tier 2).
        """
        ...

    def evaluate(self, entity: EntityRef, params: np.ndarray) -> np.ndarray:
        """Positions of *entity* at *params*: (..., 2) on a surface, (...) on a curve."""
        ...

    def fingerprint(self, entity: EntityRef) -> EntityKey:
        """Persistent geometric fingerprint of *entity* (ADR 0003)."""
        ...


class AnchoredCaster:
    """``Projector``-compatible anchored ray caster over any ``Geometry``.

    Exposes ``anchor`` and ``project(dirs, on_failure)`` (the interface the
    k=0 builders call: ``surface.cap``, ``surface.equator``,
    ``surface.remesh``), backed by ``geometry.ray_cast``. With ``snap=True``
    every HIT location is then moved to the exact closest point on the body
    (``geometry.project(..., exact=True)``): for a CAD body the ray cast runs
    on the tessellation, and snapping puts k=0 on the exact surface. For an
    STL body snapping would be a no-op, so it is left off.

    Parameters
    ----------
    geometry :
        Body geometry.
    anchor :
        Ray origin inside the body, shape (3,).
    snap :
        Snap hits to the exact surface.
    """

    def __init__(self, geometry: Geometry, anchor: np.ndarray, snap: bool = False) -> None:
        self.geometry = geometry
        self.anchor = np.asarray(anchor, dtype=np.float64).copy()
        self.snap = snap

    def project(self, dirs: np.ndarray, on_failure: str = "raise") -> ProjectionResult:
        """Cast from the anchor along *dirs* (..., 3); see ``Projector.project``."""
        res = self.geometry.ray_cast(self.anchor, dirs, on_failure)
        if not self.snap:
            return res
        valid = np.all(np.isfinite(res.locations), axis=-1)
        if not valid.any():
            return res
        locations = res.locations.copy()
        distances = res.distances.copy()
        locations[valid] = self.geometry.project(locations[valid], exact=True).points
        distances[valid] = np.linalg.norm(locations[valid] - self.anchor, axis=-1)
        return ProjectionResult(locations=locations, status=res.status, distances=distances)
