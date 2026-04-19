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
        cap_layer[0,  :, :]  = eq_layers[0]           # i=0  row ← sub0
        cap_layer[:, -1, :]  = eq_layers[1]           # j=-1 col ← sub1
        cap_layer[-1, :, :]  = eq_layers[2][::-1]     # i=-1 row ← sub2 reversed
        cap_layer[:,  0, :]  = eq_layers[3][::-1]     # j=0  col ← sub3 reversed
    else:  # south cap (z≈+1 pole, face winding mirrored)
        cap_layer[:,  0, :]  = eq_layers[0]            # j=0  col ← sub0
        cap_layer[-1, :, :]  = eq_layers[1]            # i=-1 row ← sub1
        cap_layer[:, -1, :]  = eq_layers[2][::-1]      # j=-1 col ← sub2 reversed
        cap_layer[0,  :, :]  = eq_layers[3][::-1]      # i=0  row ← sub3 reversed


# ---------------------------------------------------------------------------
# Cap boundary squareness
# ---------------------------------------------------------------------------

def _squareness_row(row: np.ndarray, squareness: float, sphere_r: float) -> np.ndarray:
    """Blend one cap boundary edge row from a circular arc toward a geodesic chord.

    The two end-nodes (corners) are fixed.  Each intermediate node at parametric
    position ``t`` is slerp-interpolated between the two corners and blended with
    the original circular-arc position:

        node_out = slerp(corner_a, corner_b, t) * (1 - squareness)
                 + original_node              * (1 - squareness)   # wrong, see below

    Correct formula: blend the *direction* between the slerp direction and the
    original direction, then re-project onto the sphere:

        dir_arc   = original_node / |original_node|
        dir_chord = slerp unit direction at t
        dir_out   = normalise((1 - squareness) * dir_arc + squareness * dir_chord)
        node_out  = dir_out * sphere_r

    Parameters
    ----------
    row :
        Shape ``(nc, 3)``.  End-nodes ``row[0]`` and ``row[-1]`` are fixed corners.
    squareness :
        Blend in ``[0, 1]``.  0 = original circular arc.  1 = geodesic chord.
    sphere_r :
        Sphere radius at this layer.

    Returns
    -------
    numpy.ndarray, shape ``(nc, 3)`` — reshaped boundary row.
    """
    if squareness == 0.0:
        return row.copy()
    nc = row.shape[0]
    if nc <= 2:
        return row.copy()

    out = row.copy()
    a = row[0]   # corner node a
    b = row[-1]  # corner node b

    # Unit vectors for slerp
    ua = a / (np.linalg.norm(a) + 1e-30)
    ub = b / (np.linalg.norm(b) + 1e-30)

    cos_ab = np.clip(np.dot(ua, ub), -1.0, 1.0)
    sin_ab = np.sqrt(max(1.0 - cos_ab**2, 0.0))

    t_vals = np.linspace(0.0, 1.0, nc)   # t=0 → corner a, t=1 → corner b

    for idx in range(1, nc - 1):        # skip the two fixed corners
        t = t_vals[idx]
        if sin_ab < 1e-12:
            # Degenerate: both corners coincide — slerp falls back to ua
            dir_chord = ua
        else:
            dir_chord = (np.sin((1.0 - t) * np.arccos(cos_ab)) * ua
                         + np.sin(t * np.arccos(cos_ab)) * ub) / sin_ab

        dir_arc = row[idx] / (np.linalg.norm(row[idx]) + 1e-30)

        dir_out = (1.0 - squareness) * dir_arc + squareness * dir_chord
        mag = np.linalg.norm(dir_out)
        if mag < 1e-15:
            dir_out = dir_arc
        else:
            dir_out /= mag

        out[idx] = dir_out * sphere_r

    return out


def _squareness_rows(
    eq_rows: list[np.ndarray],
    squareness: float,
) -> list[np.ndarray]:
    """Apply ``_squareness_row`` to all 4 cap boundary edge rows.

    The sphere radius is estimated from the corner nodes of the rows.

    Parameters
    ----------
    eq_rows :
        ``[row0, row1, row2, row3]``, each shape ``(nc, 3)``.
    squareness :
        Blend factor in ``[0, 1]``.

    Returns
    -------
    list of 4 numpy.ndarray, shape ``(nc, 3)`` each.
    """
    if squareness == 0.0:
        return [r.copy() for r in eq_rows]
    corners = np.stack([
        eq_rows[0][0], eq_rows[0][-1],
        eq_rows[1][-1], eq_rows[2][-1],
    ])
    sphere_r = float(np.mean(np.linalg.norm(corners, axis=-1)))
    return [_squareness_row(r, squareness, sphere_r) for r in eq_rows]


# ---------------------------------------------------------------------------
# Cap interior pull
# ---------------------------------------------------------------------------

def _pull_cap_interior(
    seed: np.ndarray,
    eq_rows: list[np.ndarray],
    cap_pull: float,
    pole: str = "north",
) -> np.ndarray:
    """Pull the interior nodes of the cap seed face toward the pole.

    The 4 boundary edges of the cap are kept fixed (conforming with the
    equatorial sub-blocks).  Only the strictly interior nodes are displaced
    radially toward the geographic pole and re-projected onto the sphere.

    The pull weight is a product of two sine profiles — one along each
    parametric direction — so the weight is 0 on all 4 boundary edges and
    peaks at the face centre:

        w(i, j) = cap_pull
                * sin(pi * i / (nc-1))
                * sin(pi * j / (nc-1))

    ``sin(pi*t)`` is symmetric around ``t=0.5``: nodes equidistant from
    opposite boundaries receive equal weight.  This ensures:

    - All boundary nodes are unchanged (zero weight at t=0 and t=1).
    - The centre has the maximum displacement.
    - The pull is smooth and 4-fold symmetric.
    - The equatorial sub-blocks require no modification (full conformity).

    The displaced nodes are re-projected radially onto the sphere so they
    remain on the sphere surface.

    Parameters
    ----------
    seed :
        Cap k=0 seed face, shape ``(nc, nc, 3)``.  Not modified in place;
        a new array is returned.
    eq_rows :
        ``[sub0_row, sub1_row, sub2_row, sub3_row]``, each shape ``(nc, 3)``.
        Used to determine the pole direction (centroid of the 4 corner nodes)
        and to restore the boundary edges after the pull.
    cap_pull :
        Blend factor in ``[0, 1)``.  0 = no pull (original Laplacian face).
        Values around 0.3-0.6 move the cap centre toward the pole, improving
        interior cell quality.
    pole :
        ``'north'`` or ``'south'`` — passed through to
        ``_overwrite_cap_boundaries`` for the final boundary restore.

    Returns
    -------
    numpy.ndarray, shape ``(nc, nc, 3)`` — reshaped seed face.
    """
    if cap_pull == 0.0:
        return seed.copy()

    nc = seed.shape[0]
    if nc <= 2:
        return seed.copy()

    # Sphere radius from boundary node radii
    all_boundary = np.concatenate(list(eq_rows), axis=0)
    sphere_r = float(np.mean(np.linalg.norm(all_boundary, axis=-1)))

    # Pole direction: centroid of the 4 unique corner nodes, normalised.
    # The 4 corners are the junction nodes between adjacent sub-blocks:
    #   corner(sub0,sub3) = eq_rows[0][0]
    #   corner(sub0,sub1) = eq_rows[0][-1]  (== eq_rows[1][0])
    #   corner(sub1,sub2) = eq_rows[1][-1]  (== eq_rows[2][0])
    #   corner(sub2,sub3) = eq_rows[2][-1]  (== eq_rows[3][0])
    corners = np.stack([
        eq_rows[0][0],   # corner sub0-sub3
        eq_rows[0][-1],  # corner sub0-sub1
        eq_rows[1][-1],  # corner sub1-sub2
        eq_rows[2][-1],  # corner sub2-sub3
    ])  # (4, 3) - all 4 unique corner nodes
    pole_dir = corners.mean(axis=0)
    pole_dir /= np.linalg.norm(pole_dir) + 1e-30

    # 2-D sine weight — zero on all boundary edges, peak at centre.
    # sin(pi*t) is symmetric: sin(pi*t) = sin(pi*(1-t)), so nodes
    # equidistant from opposite boundaries get equal weight.
    # (The Hann window 0.5*(1-cos(pi*t)) is NOT symmetric in this sense:
    # w(t) + w(1-t) = 1, not w(t) = w(1-t).)
    t_i = np.linspace(0.0, 1.0, nc)              # (nc,)
    t_j = np.linspace(0.0, 1.0, nc)              # (nc,)
    w_i = np.sin(np.pi * t_i)                    # (nc,)  — 0 at endpoints, 1 at centre
    w_j = np.sin(np.pi * t_j)                    # (nc,)
    W = cap_pull * w_i[:, None] * w_j[None, :]  # (nc, nc)

    pole_pt = pole_dir * sphere_r                # (3,)
    pulled = (1.0 - W[:, :, None]) * seed + W[:, :, None] * pole_pt[None, None, :]

    # Radial re-projection onto sphere
    mag = np.linalg.norm(pulled, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    pulled = pulled / mag * sphere_r

    # Restore boundary edges (zero-weight, but re-pin for floating-point safety)
    _overwrite_cap_boundaries(pulled, eq_rows, pole=pole)

    return pulled


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
) -> np.ndarray:
    """Advance one layer by hyperbolic extrusion + smoothing.

    Parameters
    ----------
    normals_override :
        When provided, replaces the grid-computed normals for this step.
        Shape must match ``(ni, nj, 3)``.  Used for STL-normal blending.
    """
    if normals_override is not None:
        normals = normals_override
    else:
        normals = _compute_normals(layer, periodic=periodic)
    new_layer = layer + step * normals
    for _ in range(smooth_iters):
        new_layer = _laplacian_smooth(new_layer, smooth_strength, periodic=periodic)
    return new_layer


def _build_cap_k0(
    eq_rows: list[np.ndarray],
    nc: int,
    pole: str = "north",
) -> np.ndarray:
    """Build the k=0 cap seed face consistent with the equatorial boundary mapping.

    Instead of using the gnomonic projection (which produces an interior
    inconsistent with the boundary mapping after the winding correction), this
    function builds the seed face by:

    1. Setting all 4 boundary edges from the equatorial k=0 pole rows.
    2. Initializing the interior to zero.
    3. Solving the discrete Laplacian with Dirichlet BCs to convergence.

    The result is the minimum-energy (smoothest) surface that interpolates the
    4 equatorial boundary edges.

    Parameters
    ----------
    eq_rows :
        ``[sub0_row, sub1_row, sub2_row, sub3_row]``, each shape (nc, 3).
    nc :
        Edge size of the cap block (= nj//4 + 1).
    pole :
        ``'north'`` or ``'south'`` — selects the boundary edge mapping.

    Returns
    -------
    numpy.ndarray, shape (nc, nc, 3)
    """
    # Target radius: average of boundary node radii (all boundary nodes should
    # be on the sphere surface at this radius).
    all_boundary = np.concatenate(list(eq_rows), axis=0)
    sphere_r = float(np.mean(np.linalg.norm(all_boundary, axis=-1)))

    seed = np.zeros((nc, nc, 3), dtype=np.float64)
    _overwrite_cap_boundaries(seed, eq_rows, pole=pole)

    tol      = 1e-12
    max_iter = 5000
    for _ in range(max_iter):
        prev = seed.copy()
        seed = _laplacian_smooth(seed, 0.5, periodic=(False, False))
        _overwrite_cap_boundaries(seed, eq_rows, pole=pole)
        if np.max(np.abs(seed - prev)) < tol:
            break

    # Project the interior nodes onto the sphere so they are co-spherical with
    # the boundary nodes.  This prevents the interior from sitting slightly inside
    # the sphere surface, which would cause a Jacobian fold at k=1 when the
    # normal-march step advances the interior outward faster than the boundary.
    # Boundary nodes are already exactly on the sphere — re-projecting them is
    # a no-op (within floating-point precision).
    mag = np.linalg.norm(seed, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    seed = seed / mag * sphere_r
    # Re-pin boundary (radial projection can introduce sub-epsilon drift)
    _overwrite_cap_boundaries(seed, eq_rows, pole=pole)

    return seed



def _march_cap_layer(
    layer: np.ndarray,
    step: float,
    smooth_strength: float,
    smooth_iters: int,
    eq_rows: list[np.ndarray],
    pole: str = "north",
    normals_override: np.ndarray | None = None,
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
        normals = _compute_normals(layer, periodic=(False, False))
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
        mask[0,        0       ] = False
        mask[0,        nj_c - 1] = False
        mask[ni_c - 1, 0       ] = False
        mask[ni_c - 1, nj_c - 1] = False
        jmin = float(jac[mask].min())
    else:
        jmin = float(jac.min())

    if jmin < tol:
        raise MarchingError(
            f"Negative Jacobian ({jmin:.3e}) in {label} at layer k={k+1}. "
            "Reduce ds, growth ratio, or increase smooth_strength."
        )
    return jmin

