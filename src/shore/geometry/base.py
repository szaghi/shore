"""Geometry protocol: the only geometry API the mesh core sees (ADR 0002).

A ``Geometry`` answers batched numpy queries about a body: anchored ray casts
(the existing templates), closest-point projection onto the body or onto one
of its entities, and parametric evaluation where the backend has parameters.
Backends: ``StlGeometry`` (triangulated STL, via trimesh) and, later,
``OcctGeometry`` (STEP/IGES/BREP, via build123d/OCP).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol, runtime_checkable

import numpy as np

from shore.surface.projector import ProjectionResult


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
