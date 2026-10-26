"""Persistent entity identity: fingerprints and matching (ADR 0003).

OCCT shapes carry no identity that survives re-reading an edited CAD file,
and ``EntityRef`` indices are only valid for one load. An ``EntityKey`` is
a geometric fingerprint of an entity, normalised by the model's bounding
box so that a uniform rescale leaves it unchanged:

- ``kind`` and ``geom_type`` (``Plane``, ``Cylinder``, ``BSplineSurface``,
  ``Line``, ``Circle``, ...) must match exactly;
- ``centroid``, ``measure`` (area, length, 0 for a point) and ``extent``
  (bounding-box size) are compared numerically.

``match_keys`` re-finds stored keys in a geometry by a global one-to-one
assignment per (kind, type) group and reports each key as ``same``
(unchanged within ``tol_same``), ``changed`` (matched, but moved or
resized, as after a parametric edit) or ``lost`` (no candidate within
``tol_match``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from scipy.optimize import linear_sum_assignment

from shore.exceptions import MeshInputError
from shore.geometry.base import EntityKind, EntityRef

if TYPE_CHECKING:
    from shore.geometry.base import Geometry

MatchStatus = Literal["same", "changed", "lost"]


@dataclass(frozen=True)
class EntityKey:
    """Geometric fingerprint of one entity, normalised by the model bbox.

    Coordinates are ``(x - bbox_min) / diag``; ``measure`` is area / diag**2
    (surface), length / diag (curve), 0 (point, body uses area).
    """

    kind: EntityKind
    geom_type: str
    centroid: tuple[float, float, float]
    measure: float
    extent: tuple[float, float, float]
    name: str | None = None

    def to_json(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "kind": self.kind.value,
            "type": self.geom_type,
            "centroid": list(self.centroid),
            "measure": self.measure,
            "extent": list(self.extent),
        }
        if self.name is not None:
            out["name"] = self.name
        return out

    @classmethod
    def from_json(cls, data: dict[str, Any], *, source: str = "<key>") -> EntityKey:
        try:
            return cls(
                kind=EntityKind(data["kind"]),
                geom_type=str(data["type"]),
                centroid=_triple(data["centroid"]),
                measure=float(data["measure"]),
                extent=_triple(data["extent"]),
                name=data.get("name"),
            )
        except (KeyError, TypeError, ValueError) as err:
            raise MeshInputError(f"{source}: malformed entity key {data!r} ({err})") from err

    def distance(self, other: EntityKey) -> float:
        """Dissimilarity in normalised units; ``inf`` across kinds or types."""
        if self.kind is not other.kind or self.geom_type != other.geom_type:
            return float("inf")
        dc = np.linalg.norm(np.subtract(self.centroid, other.centroid))
        de = np.linalg.norm(np.subtract(self.extent, other.extent))
        big = max(abs(self.measure), abs(other.measure))
        dm = abs(self.measure - other.measure) / big if big > 0.0 else 0.0
        return float(dc + de + dm)


def _triple(v: Any) -> tuple[float, float, float]:
    a, b, c = (float(x) for x in v)
    return a, b, c


def make_key(
    kind: EntityKind,
    geom_type: str,
    centroid: np.ndarray,
    measure: float,
    bbox: np.ndarray,
    model_bounds: np.ndarray,
    name: str | None = None,
) -> EntityKey:
    """Normalise raw entity properties into an ``EntityKey``.

    *bbox* and *model_bounds* are ``(2, 3)`` ``[min, max]`` boxes; *measure*
    is raw area / length (0 for points).
    """
    lo = np.asarray(model_bounds[0], dtype=np.float64)
    diag = float(np.linalg.norm(np.asarray(model_bounds[1]) - lo)) or 1.0
    power = {EntityKind.SURFACE: 2, EntityKind.BODY: 2, EntityKind.CURVE: 1, EntityKind.POINT: 0}
    c = (np.asarray(centroid, dtype=np.float64) - lo) / diag
    e = (np.asarray(bbox[1], dtype=np.float64) - np.asarray(bbox[0])) / diag
    return EntityKey(
        kind=kind,
        geom_type=geom_type,
        centroid=_triple(np.round(c, 12)),
        measure=round(float(measure) / diag ** power[kind], 12),
        extent=_triple(np.round(e, 12)),
        name=name,
    )


@dataclass
class KeyMatch:
    """Where a stored key landed in a geometry."""

    key: EntityKey
    ref: EntityRef | None
    status: MatchStatus
    cost: float


@dataclass
class MatchReport:
    """Outcome of ``match_keys``: one ``KeyMatch`` per stored key, in order."""

    matches: list[KeyMatch] = field(default_factory=list)

    def count(self, status: MatchStatus) -> int:
        return sum(m.status == status for m in self.matches)

    @property
    def all_matched(self) -> bool:
        return self.count("lost") == 0

    def summary(self) -> str:
        return (
            f"{len(self.matches)} keys: {self.count('same')} same,"
            f" {self.count('changed')} changed, {self.count('lost')} lost"
        )


def match_keys(
    keys: list[EntityKey],
    geometry: Geometry,
    *,
    tol_same: float = 1e-9,
    tol_match: float = 0.25,
) -> MatchReport:
    """Re-find *keys* in *geometry* by one-to-one assignment.

    Keys are grouped by ``(kind, geom_type)``; within a group the
    assignment minimises the total ``EntityKey.distance``. A pairing with
    cost above *tol_match* is rejected (``lost``); one at or below
    *tol_same* is ``same``; anything in between is ``changed``. Names, when
    both sides have one, must agree.
    """
    report = MatchReport([KeyMatch(k, None, "lost", float("inf")) for k in keys])
    current: dict[tuple[EntityKind, str], list[tuple[EntityRef, EntityKey]]] = {}
    for kind in EntityKind:
        for ref in geometry.entities(kind):
            fp = geometry.fingerprint(ref)
            current.setdefault((fp.kind, fp.geom_type), []).append((ref, fp))
    groups: dict[tuple[EntityKind, str], list[int]] = {}
    for i, k in enumerate(keys):
        groups.setdefault((k.kind, k.geom_type), []).append(i)
    for group, idx in groups.items():
        cands = current.get(group, [])
        if not cands:
            continue
        cost = np.array([[keys[i].distance(fp) for _, fp in cands] for i in idx])
        for r, c in zip(*_assign(cost), strict=True):
            k = keys[idx[r]]
            ref, fp = cands[c]
            named_clash = k.name is not None and fp.name is not None and k.name != fp.name
            if cost[r, c] > tol_match or named_clash:
                continue
            status: MatchStatus = "same" if cost[r, c] <= tol_same else "changed"
            report.matches[idx[r]] = KeyMatch(k, ref, status, float(cost[r, c]))
    return report


def _assign(cost: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Optimal one-to-one assignment; infinite costs are never chosen."""
    finite = np.where(np.isfinite(cost), cost, 1e12)
    rows, cols = linear_sum_assignment(finite)
    keep = np.isfinite(cost[rows, cols])
    return rows[keep], cols[keep]
