"""Ray-mesh projector with cached BVH and per-ray status tracking.

Phase 0: Projector class + ProjectionResult decoupling.
Phase 2: Full status array (BACKFACE, TOO_NEAR, TOO_FAR, DEGENERATE) + on_failure.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from enum import Enum

import numpy as np
import trimesh

from shore.exceptions import MeshInputError

log = logging.getLogger(__name__)


class ProjectionStatus(Enum):
    """Per-ray outcome of a Projector.project call."""

    HIT = "hit"
    MISSED = "missed"
    BACKFACE = "backface"
    DEGENERATE = "degenerate"
    TOO_NEAR = "too_near"
    TOO_FAR = "too_far"


@dataclass
class ProjectionResult:
    """Outcome of a batch ray-cast from a single anchor point.

    Parameters
    ----------
    locations : np.ndarray, shape (n, 3)
        Hit positions for HIT rays; NaN for non-HIT rays.
    status : np.ndarray, shape (n,), dtype=object
        Per-ray ProjectionStatus value.
    distances : np.ndarray, shape (n,)
        Anchor-to-hit distance; NaN for non-HIT rays.
    """

    locations: np.ndarray
    status: np.ndarray
    distances: np.ndarray

    @property
    def hit_mask(self) -> np.ndarray:
        """Boolean mask where True means HIT."""
        return self.status == ProjectionStatus.HIT

    @property
    def all_hit(self) -> bool:
        """True iff every ray has status HIT."""
        return bool(np.all(self.hit_mask))


class Projector:
    """Ray-mesh intersector with a single cached BVH.

    Holds a ``RayMeshIntersector`` built once at construction time so that
    multiple ``project`` calls share the same acceleration structure.

    Parameters
    ----------
    mesh :
        Triangulated body surface.
    anchor :
        Ray origin, shape (3,).
    """

    def __init__(self, mesh: trimesh.Trimesh, anchor: np.ndarray) -> None:
        self.mesh = mesh
        self.anchor = np.asarray(anchor, dtype=np.float64).copy()
        self._intersector = trimesh.ray.ray_triangle.RayMeshIntersector(mesh)
        bounds = mesh.bounds  # shape (2, 3)
        self._bbox_diag = float(np.linalg.norm(bounds[1] - bounds[0]))

    def project(
        self,
        dirs: np.ndarray,
        on_failure: str = "raise",
    ) -> ProjectionResult:
        """Cast rays from anchor in directions *dirs*.

        Parameters
        ----------
        dirs :
            Unit direction vectors, shape (..., 3).  Flattened internally.
        on_failure :
            ``'raise'`` — raise MeshInputError if any ray is not HIT.
            ``'nearest'`` — fall back to trimesh closest-point for non-HIT rays
            and flag them in the status array.
            ``'skip'`` — leave non-HIT locations as NaN.

        Returns
        -------
        ProjectionResult
        """
        original_shape = dirs.shape[:-1]
        n_rays = int(np.prod(original_shape)) if original_shape else 1
        flat_dirs = dirs.reshape(n_rays, 3)

        origins = np.tile(self.anchor, (n_rays, 1))

        raw_locs, raw_ray_idx, raw_face_idx = self._intersector.intersects_location(
            ray_origins=origins,
            ray_directions=flat_dirs,
            multiple_hits=False,
        )

        # --- initialise output arrays (NaN = not yet resolved) ---
        locations = np.full((n_rays, 3), np.nan)
        distances = np.full(n_rays, np.nan)
        status = np.full(n_rays, ProjectionStatus.MISSED, dtype=object)

        if len(raw_ray_idx) > 0:
            # Detect duplicate ray indices (phantom double-hits with multiple_hits=False)
            _, counts = np.unique(raw_ray_idx, return_counts=True)
            if np.any(counts > 1):
                n_dup = int(np.sum(counts > 1))
                log.warning(
                    "%d ray(s) returned duplicate hits from intersects_location; keeping first.",
                    n_dup,
                )
                keep = np.zeros(len(raw_ray_idx), dtype=bool)
                seen: set[int] = set()
                for i, ri in enumerate(raw_ray_idx.tolist()):
                    if ri not in seen:
                        keep[i] = True
                        seen.add(ri)
                raw_locs = raw_locs[keep]
                raw_ray_idx = raw_ray_idx[keep]
                raw_face_idx = raw_face_idx[keep]

            # Back-face detection: for inside-out casting (anchor inside body),
            # a valid hit has dot(outward_normal, outward_ray) > 0.
            # A backface hit (inverted/wrong-side triangle) has dot <= 0.
            face_normals = self.mesh.face_normals[raw_face_idx]
            dots = np.einsum("ij,ij->i", face_normals, flat_dirs[raw_ray_idx])
            is_backface = dots <= 0.0

            dists = np.linalg.norm(raw_locs - self.anchor, axis=-1)
            is_too_near = dists < 0.05 * self._bbox_diag
            is_too_far = dists > 5.0 * self._bbox_diag
            has_nan = np.any(np.isnan(raw_locs), axis=-1)

            # Priority: DEGENERATE > BACKFACE > TOO_NEAR > TOO_FAR > HIT
            ray_status = np.full(len(raw_ray_idx), ProjectionStatus.HIT, dtype=object)
            ray_status[is_too_far] = ProjectionStatus.TOO_FAR
            ray_status[is_too_near] = ProjectionStatus.TOO_NEAR
            ray_status[is_backface] = ProjectionStatus.BACKFACE
            ray_status[has_nan] = ProjectionStatus.DEGENERATE

            for ri, loc, dist, s in zip(
                raw_ray_idx.tolist(), raw_locs, dists, ray_status.tolist(), strict=False
            ):
                locations[ri] = loc
                distances[ri] = dist
                status[ri] = s

        # --- handle failures ---
        not_hit_mask = status != ProjectionStatus.HIT
        n_fail = int(np.sum(not_hit_mask))

        if n_fail > 0:
            if on_failure == "raise":
                counts_map = Counter(status[not_hit_mask].tolist())
                counts_str = ", ".join(
                    f"{k.value}={v}" for k, v in counts_map.items()
                )
                raise MeshInputError(
                    f"{n_fail} of {n_rays} rays failed ({counts_str}).  "
                    "The body may not be star-shaped w.r.t. the anchor, or the "
                    "STL has holes.  Try a different anchor strategy or repair the mesh."
                )
            elif on_failure == "nearest":
                fail_idx = np.where(not_hit_mask)[0]
                query = self.anchor + flat_dirs[fail_idx] * self._bbox_diag
                close_pts, _, _ = trimesh.proximity.closest_point(self.mesh, query)
                locations[fail_idx] = close_pts
                distances[fail_idx] = np.linalg.norm(
                    close_pts - self.anchor, axis=-1
                )
            # on_failure == "skip": leave locations as NaN

        return ProjectionResult(
            locations=locations,
            status=status,
            distances=distances,
        )
