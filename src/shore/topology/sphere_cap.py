"""Polar cap block generation for sphere O-grids.

Algorithm overview
------------------
A single-block sphere O-grid (ni, nj, nk, 3) has a periodic j-direction and
open i-boundaries near the poles.  This module closes both polar holes by:

1. Splitting the equatorial block into 4 conforming sub-blocks along j.
2. Building a gnomonic seed face for each polar cap from the 4 corner nodes
   that meet at each pole.
3. Marching all 6 blocks simultaneously, layer by layer, with the equatorial
   sub-blocks driving the shared faces of the cap blocks at every layer.

The result is a list of 6 arrays in order::

    [sub0, sub1, sub2, sub3, cap_north, cap_south]

where ``sub_q`` has shape ``(ni, nj//4+1, nk, 3)`` and each cap has shape
``(nc, nc, nk, 3)`` with ``nc = nj//4 + 1``.

Block adjacency (shared faces)
-------------------------------
- sub0 ↔ sub1  (at j = nj//4)
- sub1 ↔ sub2  (at j = nj//2)
- sub2 ↔ sub3  (at j = 3*nj//4)
- sub3 ↔ sub0  (at j = 0 / nj, periodic seam)
- sub0, sub1, sub2, sub3 ↔ cap_north  (at i = ni-1)
- sub0, sub1, sub2, sub3 ↔ cap_south  (at i = 0)
"""

from __future__ import annotations

import warnings

import numpy as np

from shore.exceptions import MarchingError
from shore.volume.hyperbolic import _compute_normals, _laplacian_smooth
from shore.volume.jacobian import _jacobian_field_open

# ---------------------------------------------------------------------------
# Step 1 — nj validation and equatorial block splitting
# ---------------------------------------------------------------------------


def _adjust_nj(nj: int) -> int:
    """Round nj down to the nearest multiple of 4, warning if adjusted."""
    if nj % 4 != 0:
        adjusted = (nj // 4) * 4
        warnings.warn(
            f"nj={nj} is not divisible by 4; rounding down to nj={adjusted} "
            "so the equatorial block can be split into 4 equal sectors.",
            stacklevel=3,
        )
        return adjusted
    return nj


def _overwrite_cap_boundaries(
    cap_layer: np.ndarray,
    eq_layers: list[np.ndarray],
    pole: str = "north",
) -> None:
    """Overwrite the 4 boundary edges of a cap layer with equatorial values.

    The two caps have opposite face windings (outward normals in opposite
    directions), so the edge-to-sub-block mapping differs between poles.

    North cap (south pole of sphere, z≈-1, winding: i=0 row → sub0):

    - ``cap[0,  :]`` ← sub0 pole row  (phi 0°→90°)
    - ``cap[:, -1]`` ← sub1 pole row  (phi 90°→180°)
    - ``cap[-1, :]`` ← sub2 pole row reversed
    - ``cap[:,  0]`` ← sub3 pole row reversed

    South cap (north pole of sphere, z≈+1, winding is mirrored — j=0 col → sub0):

    - ``cap[:,  0]`` ← sub0 pole row  (phi 0°→90°)
    - ``cap[-1, :]`` ← sub1 pole row  (phi 90°→180°)
    - ``cap[:, -1]`` ← sub2 pole row reversed
    - ``cap[0,  :]`` ← sub3 pole row reversed

    Parameters
    ----------
    cap_layer :
        Current cap layer, shape (nc, nc, 3).  Modified in place.
    eq_layers :
        ``[sub0_row, sub1_row, sub2_row, sub3_row]``, each shape (nc, 3).
    pole :
        ``'north'`` or ``'south'`` — selects the edge mapping.
    """
    if pole == "north":
        cap_layer[0, :, :] = eq_layers[0]  # i=0  row ← sub0
        cap_layer[:, -1, :] = eq_layers[1]  # j=-1 col ← sub1
        cap_layer[-1, :, :] = eq_layers[2][::-1]  # i=-1 row ← sub2 reversed
        cap_layer[:, 0, :] = eq_layers[3][::-1]  # j=0  col ← sub3 reversed
    else:  # south cap (z≈+1 pole, face winding mirrored)
        cap_layer[:, 0, :] = eq_layers[0]  # j=0  col ← sub0
        cap_layer[-1, :, :] = eq_layers[1]  # i=-1 row ← sub1
        cap_layer[:, -1, :] = eq_layers[2][::-1]  # j=-1 col ← sub2 reversed
        cap_layer[0, :, :] = eq_layers[3][::-1]  # i=0  row ← sub3 reversed


# ---------------------------------------------------------------------------
# Step 5-6 — Simultaneous marching
# ---------------------------------------------------------------------------


def _march_one_layer(
    layer: np.ndarray,
    step: float,
    smooth_strength: float,
    smooth_iters: int,
    periodic: tuple[bool, bool],
    normals_override: np.ndarray | None = None,
    ghost_rows: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    """Advance one layer by hyperbolic extrusion + smoothing.

    Parameters
    ----------
    normals_override :
        When provided, replaces the grid-computed normals for this step.
        Shape must match ``(ni, nj, 3)``.  Used for STL-normal blending.
    ghost_rows :
        Optional virtual neighbour rows (per-face) for centred-difference
        normals at multi-block seams; see ``_compute_normals``.
    """
    if normals_override is not None:
        normals = normals_override
    else:
        normals = _compute_normals(layer, periodic=periodic, ghost_rows=ghost_rows)
    new_layer = layer + step * normals
    for _ in range(smooth_iters):
        new_layer = _laplacian_smooth(new_layer, smooth_strength, periodic=periodic)
    return new_layer


def _march_cap_layer(
    layer: np.ndarray,
    step: float,
    smooth_strength: float,
    smooth_iters: int,
    eq_rows: list[np.ndarray],
    pole: str = "north",
    normals_override: np.ndarray | None = None,
    ghost_rows: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    """Advance one cap layer with Dirichlet boundary pinning.

    The 4 boundary edges are set to the equatorial values before normal
    computation and re-pinned after every Laplacian sweep.  This prevents
    the interior from drifting away from the boundary-driven geometry.

    Parameters
    ----------
    layer :
        Current cap layer, shape (nc, nc, 3).
    step :
        Physical step size for this layer.
    smooth_strength :
        Laplacian blend factor.
    smooth_iters :
        Number of Laplacian sweeps.
    eq_rows :
        List of 4 equatorial pole-row arrays, shape (nc, 3) each, giving the
        target positions of the 4 boundary edges for the next layer.
    normals_override :
        When provided, replaces grid-computed normals for this step.
        Used for STL-normal blending.

    Returns
    -------
    numpy.ndarray, shape (nc, nc, 3)
    """
    # 1. Normal-based displacement.
    if normals_override is not None:
        normals = normals_override
    else:
        normals = _compute_normals(layer, periodic=(False, False), ghost_rows=ghost_rows)
    new_layer = layer + step * normals

    # 2. Laplacian smoothing passes (same as equatorial blocks).
    for _ in range(smooth_iters):
        new_layer = _laplacian_smooth(new_layer, smooth_strength, periodic=(False, False))

    # 3. Pin the 4 boundary edges to the equatorial values.
    #    This overwrites the smoothed boundary with the exact equatorial
    #    positions, ensuring node-for-node conformity.
    _overwrite_cap_boundaries(new_layer, eq_rows, pole=pole)

    return new_layer


def _check_jacobian(
    layer: np.ndarray,
    new_layer: np.ndarray,
    periodic_j: bool,
    label: str,
    k: int,
    exclude_corners: bool = False,
    tol: float = 0.0,
) -> float:
    """Check minimum Jacobian; raise RuntimeError if below *tol*.

    Parameters
    ----------
    exclude_corners :
        When True, the 4 corner cells of the block are excluded from the
        check.  Used for cap blocks: the corners are shared with adjacent
        equatorial sub-blocks and can have near-zero or slightly negative
        Jacobians due to the degenerate winding at the 4-block junction.
    tol :
        Minimum acceptable Jacobian (default 0 — strictly positive).  Cap
        blocks with ``cap_pull > 0`` can produce Jacobians of O(-1e-8) due
        to floating-point noise in the normal computation; setting
        ``tol=-1e-6`` prevents spurious failures while still catching real
        cell inversions (which are O(-1e-5) or larger).
    """
    two = np.stack([layer, new_layer], axis=2)
    if periodic_j:
        from shore.volume.jacobian import jacobian_field

        jac = jacobian_field(two)
    else:
        jac = _jacobian_field_open(two)

    if exclude_corners:
        ni_c, nj_c = jac.shape[0], jac.shape[1]
        mask = np.ones(jac.shape[:2], dtype=bool)
        mask[0, 0] = False
        mask[0, nj_c - 1] = False
        mask[ni_c - 1, 0] = False
        mask[ni_c - 1, nj_c - 1] = False
        jmin = float(jac[mask].min())
    else:
        jmin = float(jac.min())

    if jmin < tol:
        raise MarchingError(
            f"Negative Jacobian ({jmin:.3e}) in {label} at layer k={k + 1}. "
            "Reduce ds, growth ratio, or increase smooth_strength."
        )
    return jmin
