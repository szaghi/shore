"""Anchor-point resolution strategies for ray-cast projections.

Phase 1: resolve_anchor with centroid / bbox / incenter / user / auto strategies.
"""

from __future__ import annotations

import logging

import numpy as np
import trimesh

from shore.exceptions import MeshInputError, ParameterError

log = logging.getLogger(__name__)

_STRATEGIES = frozenset({"auto", "centroid", "bbox", "incenter", "user"})


def resolve_anchor(
    mesh: trimesh.Trimesh,
    strategy: str = "auto",
    user_point: np.ndarray | None = None,
) -> np.ndarray:
    """Resolve a validated ray-cast anchor point for *mesh*.

    Parameters
    ----------
    mesh :
        Triangulated body surface (already hygiene-processed is recommended).
    strategy :
        ``'centroid'`` — volumetric centroid (legacy; may lie outside on
        non-convex bodies).
        ``'bbox'`` — bounding-box centre (always inside a convex body).
        ``'incenter'`` — centre of the minimum enclosing sphere
        (``trimesh.nsphere.minimum_nsphere``); robust for smooth convex bodies.
        ``'user'`` — caller-supplied *user_point*; validated before use.
        ``'auto'`` — try incenter, validate, fall back to bbox, validate, raise.
    user_point :
        Required when *strategy* is ``'user'``.

    Returns
    -------
    np.ndarray, shape (3,)
        Validated anchor point confirmed inside the mesh.

    Raises
    ------
    ParameterError
        If *strategy* is unknown or *user_point* is absent for ``'user'``.
    MeshInputError
        If the resolved anchor fails the containment check, or if ``'auto'``
        exhausts all fallbacks.
    """
    if strategy not in _STRATEGIES:
        raise ParameterError(
            f"Unknown anchor strategy {strategy!r}.  "
            f"Choose from {sorted(_STRATEGIES)}."
        )

    if strategy == "user":
        if user_point is None:
            raise ParameterError(
                "strategy='user' requires user_point to be supplied."
            )
        anchor = np.asarray(user_point, dtype=np.float64).copy()
        _validate_anchor(mesh, anchor, label="user-supplied")
        return anchor

    if strategy == "centroid":
        anchor = mesh.centroid.copy()
        _validate_anchor(mesh, anchor, label="centroid")
        return anchor

    if strategy == "bbox":
        anchor = _bbox_center(mesh)
        _validate_anchor(mesh, anchor, label="bbox")
        return anchor

    if strategy == "incenter":
        anchor = _incenter(mesh)
        _validate_anchor(mesh, anchor, label="incenter")
        return anchor

    # strategy == "auto": try incenter → bbox → raise
    for label, fn in (("incenter", _incenter), ("bbox", _bbox_center)):
        try:
            anchor = fn(mesh)
            _validate_anchor(mesh, anchor, label=label)
            log.debug("auto anchor resolved via '%s': %s", label, anchor)
            return anchor
        except MeshInputError:
            log.debug("anchor strategy '%s' failed containment; trying next.", label)
            continue

    raise MeshInputError(
        "No valid anchor found via 'auto' strategy (tried incenter and bbox).  "
        "The mesh may not be star-shaped.  Supply anchor_strategy='user' and "
        "provide a known interior point."
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _bbox_center(mesh: trimesh.Trimesh) -> np.ndarray:
    return np.mean(mesh.bounds, axis=0).copy()


def _incenter(mesh: trimesh.Trimesh) -> np.ndarray:
    """Centre of the minimum enclosing sphere of the mesh vertices."""
    center, _ = trimesh.nsphere.minimum_nsphere(mesh)
    return np.asarray(center, dtype=np.float64).copy()


def _validate_anchor(
    mesh: trimesh.Trimesh,
    anchor: np.ndarray,
    label: str,
) -> None:
    """Raise MeshInputError if *anchor* is not inside *mesh*."""
    inside = mesh.contains(anchor.reshape(1, 3))
    if not bool(inside[0]):
        raise MeshInputError(
            f"Anchor ({label}) at {anchor} is outside the mesh.  "
            "The mesh may not be watertight, or the body is non-convex."
        )
