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

from shore.exceptions import MarchingError, ParameterError
from shore.volume.hyperbolic import (
    SMOOTH_MODES,
    _compute_normals,
    _laplacian_smooth,
    _limit_smoothing,
)
from shore.volume.hyperbolic_scheme import hyperbolic_layer
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
    ghost_rows: dict[str, np.ndarray] | None = None,
    smooth_mode: str = "limited",
    march_scheme: str = "algebraic",
) -> np.ndarray:
    """Advance one layer by normal extrusion + smoothing, or one hyperbolic step.

    Parameters
    ----------
    ghost_rows :
        Optional virtual neighbour rows (per-face) for centred-difference
        normals at multi-block seams; see ``_compute_normals``.
    smooth_mode :
        ``'tangential'`` (default): keep only the tangential part of the
        smoothing motion, capped by the step (``_limit_smoothing``).
        ``'legacy'``: the full Laplacian motion (folds thin first layers).
    """
    if march_scheme == "hyperbolic":
        return hyperbolic_layer(layer, step, periodic, ghost_rows)
    normals = _compute_normals(layer, periodic=periodic, ghost_rows=ghost_rows)
    raw = layer + step * normals
    new_layer = raw
    for _ in range(smooth_iters):
        new_layer = _laplacian_smooth(new_layer, smooth_strength, periodic=periodic)
    return _apply_smooth_mode(raw, new_layer, normals, step, smooth_mode)


def _apply_smooth_mode(
    raw: np.ndarray, smoothed: np.ndarray, normals: np.ndarray, step: float, mode: str
) -> np.ndarray:
    if mode == "limited":
        return _limit_smoothing(raw, smoothed, normals, step)
    if mode == "legacy":
        return smoothed
    raise ParameterError(f"smooth_mode must be one of {SMOOTH_MODES}, got {mode!r}.")


def _march_cap_layer(
    layer: np.ndarray,
    step: float,
    smooth_strength: float,
    smooth_iters: int,
    eq_rows: list[np.ndarray],
    pole: str = "north",
    ghost_rows: dict[str, np.ndarray] | None = None,
    smooth_mode: str = "limited",
    march_scheme: str = "algebraic",
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

    Returns
    -------
    numpy.ndarray, shape (nc, nc, 3)
    """
    if march_scheme == "hyperbolic":
        # the pinned equator rows are the Dirichlet boundary of the implicit solve
        target = layer.copy()
        _overwrite_cap_boundaries(target, eq_rows, pole=pole)
        new_layer = hyperbolic_layer(
            layer, step, (False, False), ghost_rows, boundary_target=target
        )
        _overwrite_cap_boundaries(new_layer, eq_rows, pole=pole)
        return new_layer

    # 1. Normal-based displacement.
    normals = _compute_normals(layer, periodic=(False, False), ghost_rows=ghost_rows)
    raw = layer + step * normals

    # 2. Laplacian smoothing passes (same as equatorial blocks), constrained
    #    by smooth_mode before the boundary re-pin below.
    new_layer = raw
    for _ in range(smooth_iters):
        new_layer = _laplacian_smooth(new_layer, smooth_strength, periodic=(False, False))
    new_layer = _apply_smooth_mode(raw, new_layer, normals, step, smooth_mode)

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


def _check_folds(layer: np.ndarray, new_layer: np.ndarray, label: str, k: int) -> None:
    """Raise if a cell between *layer* and *new_layer* is folded (ADR 0052).

    Every cell counts, the block's corner cells at the 3-block junctions
    included, by its minimum corner scaled Jacobian (``shore.quality``):
    scale-free, so thin wall cells are judged like thick ones.
    ``_check_jacobian`` skips those corner cells and tests one corner against
    an absolute tolerance; folds there went through unreported.
    """
    from shore.quality import scaled_jacobian

    sj = scaled_jacobian(np.stack([layer, new_layer], axis=2))[..., 0]
    if sj.min() > 0.0:
        return
    i, j = np.unravel_index(int(np.argmin(sj)), sj.shape)
    raise MarchingError(
        f"Folded cell (scaled Jacobian {float(sj.min()):.3f}) in {label} at layer "
        f"k={k + 1}, cell ({i}, {j}), {int((sj <= 0.0).sum())} folded in the layer. "
        "Cubed-sphere folds start at the cap corners: use a larger cap "
        "(theta_cap 30 deg or more), march_scheme='hyperbolic', or a smaller growth."
    )
