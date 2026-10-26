"""Block-to-geometry association (ADR 0003).

Every vertex, edge and face of a block may be *associated* with geometry: the
whole body, a set of CAD surfaces, a chain of CAD curves, or a CAD point. The
nodes of an associated entity are projected onto its geometry after any
operation that moves them, with priority vertex -> edge -> face (ICEM's
projection order): face nodes first, then edge nodes, then the corners, so
the lower-dimensional constraint wins where they meet.

Associations name geometry persistently: ``geometry`` is a key into the
mesh's geometries (``"body"`` for the meshed body) and ``keys`` are
``EntityKey`` fingerprints, re-found in the current geometry by
``match_keys``. An association of kind ``BODY`` has no keys.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from shore.exceptions import MeshInputError, ParameterError
from shore.geometry.base import BODY, EntityKind, EntityRef
from shore.geometry.keys import EntityKey, MatchReport, match_keys

if TYPE_CHECKING:
    from shore.geometry.base import Geometry
    from shore.mesh.block import HexBlock

BODY_GEOMETRY = "body"


@dataclass(frozen=True)
class Association:
    """Persistent link from a block entity to geometry.

    Parameters
    ----------
    geometry :
        Key of the geometry in the mesh's geometry map (``"body"`` for the
        meshed body).
    kind :
        ``BODY`` (whole body), ``SURFACE`` (surface set), ``CURVE`` (curve
        chain) or ``POINT``.
    keys :
        Fingerprints of the associated entities; empty for ``BODY``.
    """

    geometry: str
    kind: EntityKind
    keys: tuple[EntityKey, ...] = ()

    def __post_init__(self) -> None:
        if self.kind is EntityKind.BODY and self.keys:
            raise ParameterError("A BODY association carries no entity keys.")
        if self.kind is not EntityKind.BODY and not self.keys:
            raise ParameterError(f"A {self.kind.value} association needs entity keys.")
        if any(k.kind is not self.kind for k in self.keys):
            raise ParameterError(f"All keys of a {self.kind.value} association must match it.")

    @classmethod
    def body(cls, geometry: str = BODY_GEOMETRY) -> Association:
        """Association with the whole body."""
        return cls(geometry, EntityKind.BODY)

    @classmethod
    def of(cls, geometry_name: str, geometry: Geometry, refs: list[EntityRef]) -> Association:
        """Association with entities *refs* of *geometry* (all of one kind)."""
        if not refs:
            raise ParameterError("Association.of needs at least one entity.")
        kind = refs[0].kind
        if kind is EntityKind.BODY:
            return cls.body(geometry_name)
        return cls(geometry_name, kind, tuple(geometry.fingerprint(r) for r in refs))

    def to_json(self) -> dict[str, Any]:
        return {
            "geometry": self.geometry,
            "kind": self.kind.value,
            "keys": [k.to_json() for k in self.keys],
        }

    @classmethod
    def from_json(cls, data: dict[str, Any], *, source: str = "<association>") -> Association:
        try:
            keys = tuple(EntityKey.from_json(k, source=source) for k in data.get("keys", []))
            return cls(str(data["geometry"]), EntityKind(data["kind"]), keys)
        except (KeyError, ValueError, ParameterError) as err:
            raise MeshInputError(f"{source}: malformed association {data!r} ({err})") from err


class AssociationResolver:
    """Resolves associations to current ``EntityRef``\\ s and projects nodes.

    Parameters
    ----------
    geometries :
        Geometry map, keyed like ``Association.geometry``.
    exact :
        Project with ``exact=True`` (tier 2) instead of the tessellation.
    """

    def __init__(self, geometries: dict[str, Geometry], *, exact: bool = False) -> None:
        self.geometries = geometries
        self.exact = exact
        self._refs: dict[Association, list[EntityRef]] = {}
        self.reports: dict[Association, MatchReport] = {}

    def refs(self, assoc: Association) -> list[EntityRef]:
        """Current entities of *assoc*; raises if any key is lost."""
        if assoc not in self._refs:
            geometry = self._geometry(assoc)
            if assoc.kind is EntityKind.BODY:
                self._refs[assoc] = [BODY]
            else:
                report = match_keys(list(assoc.keys), geometry)
                self.reports[assoc] = report
                if not report.all_matched:
                    raise MeshInputError(
                        f"association with '{assoc.geometry}' no longer resolves"
                        f" ({report.summary()}); re-associate the entity."
                    )
                self._refs[assoc] = [m.ref for m in report.matches]  # type: ignore[misc]
        return self._refs[assoc]

    def project(self, assoc: Association, points: np.ndarray) -> np.ndarray:
        """Project *points* (..., 3) onto *assoc*; nearest entity of a set wins."""
        geometry = self._geometry(assoc)
        flat = np.asarray(points, dtype=np.float64).reshape(-1, 3)
        best_p = flat.copy()
        best_d = np.full(len(flat), np.inf)
        for ref in self.refs(assoc):
            proj = geometry.project(flat, ref, exact=self.exact)
            closer = proj.distances < best_d
            best_p[closer] = proj.points[closer]
            best_d[closer] = proj.distances[closer]
        return best_p.reshape(np.shape(points))

    def _geometry(self, assoc: Association) -> Geometry:
        try:
            return self.geometries[assoc.geometry]
        except KeyError as err:
            raise MeshInputError(
                f"association refers to geometry '{assoc.geometry}', which this mesh lacks"
            ) from err


def snap_block(block: HexBlock, resolver: AssociationResolver) -> None:
    """Project every associated node of *block* onto its geometry, in place.

    Faces first, then edges, then vertices, so that the lower-dimensional
    association wins at shared nodes (vertex -> edge -> face priority).
    """
    from shore.mesh.edge import _edge_slice
    from shore.mesh.face import _face_slice

    for face in block.faces():
        if face.association is not None:
            view = _face_slice(block, face.name)
            view[...] = resolver.project(face.association, view)
    for edge in block.edges():
        if edge.association is not None:
            view = _edge_slice(block, edge.name)
            view[...] = resolver.project(edge.association, view)
    for vertex in block.vertices():
        if vertex.association is not None:
            i, j, k = vertex.ijk
            block.nodes[i, j, k] = resolver.project(vertex.association, block.nodes[i, j, k])


# k_lo face, its 4 edges and 4 corners (k = 0)
_WALL_EDGES = ("j_lo_k_lo", "j_hi_k_lo", "i_lo_k_lo", "i_hi_k_lo")


def associate_wall(block: HexBlock, assoc: Association) -> None:
    """Associate *block*'s k_lo face, its 4 edges and 4 corners with *assoc*."""
    block.face("k_lo").association = assoc
    for edge in block.edges():
        if edge.name in _WALL_EDGES:
            edge.association = assoc
    for vertex in block.vertices():
        if vertex.ijk[2] == 0:
            vertex.association = assoc
