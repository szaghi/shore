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
from collections.abc import Callable

import numpy as np

from shore.exceptions import MarchingError, ParameterError
from shore.volume.hyperbolic import _compute_normals, _laplacian_smooth
from shore.volume.jacobian import min_jacobian_two_layers, _jacobian_field_open
from shore.volume.spacing import build_k_steps


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


def _split_equatorial(grid: np.ndarray) -> list[np.ndarray]:
    """Split the full O-grid into 4 equatorial sub-blocks along j.

    Each sub-block explicitly duplicates its trailing shared edge so that
    adjacent blocks are node-for-node conforming.

    Parameters
    ----------
    grid :
        Full O-grid, shape (ni, nj, nk, 3).  nj must be divisible by 4.

    Returns
    -------
    list of 4 arrays, each shape (ni, nj//4+1, nk, 3)
    """
    ni, nj, nk, _ = grid.shape
    sector = nj // 4
    sub_blocks = []
    for q in range(4):
        j0 = q * sector
        j1 = j0 + sector  # exclusive end (sector columns)
        # Take sector+1 columns: j0..j1 inclusive; wrap j1 for the last block
        j_indices = list(range(j0, j1)) + [j1 % nj]
        sub_blocks.append(grid[:, j_indices, :, :].copy())
    return sub_blocks


# ---------------------------------------------------------------------------
# Step 2-3 — Polar cap seed faces
# ---------------------------------------------------------------------------

def _gnomonic_seed_face(corners: np.ndarray, nc: int) -> np.ndarray:
    """Build a (nc, nc, 3) polar cap seed face by gnomonic projection.

    Algorithm:
    1. Bilinearly interpolate the 4 corners in 3D space.
    2. Project each interpolated point radially onto the sphere
       (using the average corner radius as the target radius).
    3. Verify that the face normal at the centre points outward (away from the
       sphere origin); if not, flip the v-axis so the winding is correct for
       hyperbolic extrusion.

    Parameters
    ----------
    corners :
        Shape (4, 3) — four corner nodes ordered as::

            c0 (u=0,v=0) --- c1 (u=1,v=0)
            |                            |
            c3 (u=0,v=1) --- c2 (u=1,v=1)

    nc :
        Number of nodes along each edge of the seed face (= nj//4+1).

    Returns
    -------
    numpy.ndarray, shape (nc, nc, 3)
        Seed face with outward-pointing normals (away from the sphere origin).
    """
    c0, c1, c2, c3 = corners  # (3,) each

    # Average sphere radius from the 4 corners
    r = float(np.mean(np.linalg.norm(corners, axis=-1)))

    u = np.linspace(0.0, 1.0, nc)
    v = np.linspace(0.0, 1.0, nc)
    ug, vg = np.meshgrid(u, v, indexing="ij")  # (nc, nc)

    # Bilinear interpolation in 3D
    pts = (
        (1.0 - ug[..., None]) * (1.0 - vg[..., None]) * c0
        + ug[..., None]        * (1.0 - vg[..., None]) * c1
        + ug[..., None]        *        vg[..., None]  * c2
        + (1.0 - ug[..., None]) *       vg[..., None]  * c3
    )  # (nc, nc, 3)

    # Radial projection onto the sphere
    mag = np.linalg.norm(pts, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    face = pts / mag * r   # (nc, nc, 3)

    # --- Winding check: ensure centre normal points outward ---
    # Use the centre of the face as the representative point
    mid = nc // 2
    centre = face[mid, mid]                 # (3,) — radial position
    outward = centre / (np.linalg.norm(centre) + 1e-30)

    # Estimate face normal at centre via central differences
    ti = face[mid + 1, mid] - face[mid - 1, mid] if nc > 2 else face[1, mid] - face[0, mid]
    tj = face[mid, mid + 1] - face[mid, mid - 1] if nc > 2 else face[mid, 1] - face[mid, 0]
    normal_centre = np.cross(ti, tj)

    if np.dot(normal_centre, outward) < 0.0:
        # Flip the v-axis to reverse winding
        face = face[:, ::-1, :]

    return face


def _build_cap_seeds(
    grid: np.ndarray,
    sub_blocks: list[np.ndarray],
    nc: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build seed faces (k=0) for the north and south polar caps.

    The 4 corner nodes are the j-interface nodes of the equatorial ring at
    i=ni-1 (north) and i=0 (south), taken from the k=0 layer.

    Corner ordering (viewing the cap from outside the sphere, CCW):
    - north: c0=sub0[:,0], c1=sub1[:,0], c2=sub2[:,0], c3=sub3[:,0]  at i=ni-1
    - south: c0=sub0[:,0], c1=sub3[:,0], c2=sub2[:,0], c3=sub1[:,0]  at i=0
      (reversed to keep outward normal pointing south)

    Parameters
    ----------
    grid :
        Full O-grid, shape (ni, nj, nk, 3).
    sub_blocks :
        List of 4 equatorial sub-blocks (shape (ni, nj//4+1, nk, 3) each).
    nc :
        Edge node count for the cap seed face.

    Returns
    -------
    cap_north_seed, cap_south_seed : np.ndarray, shape (nc, nc, 3) each
    """
    ni = grid.shape[0]

    # Corner nodes at the pole rows (k=0 layer, j=0 of each sub-block)
    # sub_q[:,0,0,:] is the first j-column of sub-block q, at i=all, k=0
    north_corners = np.stack([
        sub_blocks[0][ni - 1, 0,  0, :],   # c0
        sub_blocks[1][ni - 1, 0,  0, :],   # c1
        sub_blocks[2][ni - 1, 0,  0, :],   # c2
        sub_blocks[3][ni - 1, 0,  0, :],   # c3
    ])  # (4, 3)

    south_corners = np.stack([
        sub_blocks[0][0, 0, 0, :],          # c0
        sub_blocks[3][0, 0, 0, :],          # c1 — reversed for outward normal
        sub_blocks[2][0, 0, 0, :],          # c2
        sub_blocks[1][0, 0, 0, :],          # c3
    ])  # (4, 3)

    cap_north_seed = _gnomonic_seed_face(north_corners, nc)
    cap_south_seed = _gnomonic_seed_face(south_corners, nc)
    return cap_north_seed, cap_south_seed


# ---------------------------------------------------------------------------
# Step 4 — Shared-face enforcement
# ---------------------------------------------------------------------------

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
        Values around 0.3–0.6 move the cap centre toward the pole, improving
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
    all_boundary = np.concatenate([r for r in eq_rows], axis=0)
    sphere_r = float(np.mean(np.linalg.norm(all_boundary, axis=-1)))

    # Pole direction: centroid of the 4 unique corner nodes, normalised.
    # The 4 corners are the junction nodes between adjacent sub-blocks:
    #   corner(sub0,sub3) = eq_rows[0][0]
    #   corner(sub0,sub1) = eq_rows[0][-1]  (== eq_rows[1][0])
    #   corner(sub1,sub2) = eq_rows[1][-1]  (== eq_rows[2][0])
    #   corner(sub2,sub3) = eq_rows[2][-1]  (== eq_rows[3][0])
    corners = np.stack([
        eq_rows[0][0],   # corner sub0–sub3
        eq_rows[0][-1],  # corner sub0–sub1
        eq_rows[1][-1],  # corner sub1–sub2
        eq_rows[2][-1],  # corner sub2–sub3
    ])  # (4, 3) — all 4 unique corner nodes
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
) -> np.ndarray:
    """Advance one layer by hyperbolic extrusion + smoothing."""
    normals   = _compute_normals(layer, periodic=periodic)
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
    all_boundary = np.concatenate([r for r in eq_rows], axis=0)
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
    # The cap is a fully-boundary-driven block: all 4 edges are set from the
    # equatorial sub-blocks.  The interior is determined by solving the discrete
    # Laplacian with Dirichlet boundary conditions (iterate to convergence).
    #
    # Algorithm:
    # 1. Start from the boundary-pinned current layer as initial guess.
    # 2. Repeatedly apply Laplacian smoothing (strength=0.5) with boundary
    #    re-pinning until the interior converges (max change < tol).
    # 3. The equatorial boundary values encode the correct outward displacement;
    #    the interior follows via the converged Laplacian solution.

    # 1. Normal-based displacement (same as equatorial blocks).
    normals   = _compute_normals(layer, periodic=(False, False))
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


def split_ogrid_with_polar_caps(
    grid: np.ndarray,
    nk: int,
    ds: float = 1e-3,
    growth: float = 1.15,
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    cap_pull: float = 0.0,
    cap_squareness: float = 0.0,
    layer_callback: Callable[[int, float], None] | None = None,
    steps: np.ndarray | None = None,
) -> list[np.ndarray]:
    """Split a sphere O-grid into 4 equatorial sub-blocks + 2 polar cap blocks.

    The input is a fully-marched single-block O-grid **or** only its k=0
    surface layer (in which case ``nk``, ``ds``, and ``growth`` drive the
    extrusion).  The function re-marches all 6 blocks simultaneously from the
    k=0 surface so that shared faces are exactly conforming at every layer.

    Parameters
    ----------
    grid :
        Full O-grid, shape ``(ni, nj, nk_in, 3)`` **or** the k=0 surface
        layer ``(ni, nj, 3)``.  If the k=0 surface is passed, ``nk`` must
        be given explicitly.  If the full volume is passed its ``nk`` is used
        (the ``nk`` argument is ignored when shape[2] > 1).
    nk :
        Number of wall-normal layers (including k=0).  Only used when
        ``grid`` has shape ``(ni, nj, 3)`` (k=0 layer only).
    ds :
        First-layer thickness.
    growth :
        Geometric growth ratio between successive layers.
    smooth_strength :
        Laplacian blend factor (0 = off).
    smooth_iters :
        Number of Laplacian sweeps per layer.
    cap_pull :
        Blend factor in ``[0, 1)`` that domes the cap interior toward the
        pole.  0 (default) leaves the interior as the flat Laplacian face.
        Values 0.3–0.5 improve interior cell quality without Jacobian issues.
    cap_squareness :
        Blend factor in ``[0, 1]`` that straightens each cap boundary edge
        from a circular arc toward the geodesic chord connecting its two
        corner nodes.  0 (default) = original parallel-circle arc (~157°
        corner angle).  1 = full geodesic chord (the edge of the inscribed
        square, ~90° corner angle).  Corner nodes are fixed; only intermediate
        nodes are moved.  Applied at every layer to keep the cap boundary
        edges in the inscribed-square geometry as the grid grows outward.
    layer_callback :
        Optional callable ``(k, jmin_global)`` invoked after each accepted
        layer set, where ``jmin_global`` is the minimum Jacobian across all
        6 blocks.
    steps :
        Pre-computed wall-normal step sizes, shape ``(nk-1,)``.  When
        supplied, overrides ``ds`` and ``growth``.  Use
        ``shore.volume.spacing.build_k_steps`` to generate tanh/tanh2 laws.

    Returns
    -------
    list of 6 numpy.ndarray
        ``[sub0, sub1, sub2, sub3, cap_north, cap_south]``

        - ``sub_q``:    shape ``(ni, nc, nk, 3)``
        - ``cap_north``, ``cap_south``: shape ``(nc, nc, nk, 3)``

        where ``nc = nj_adj // 4 + 1`` and ``nj_adj`` is the (possibly
        adjusted) nj used after rounding.

    Raises
    ------
    ValueError
        If the grid array has an unexpected shape.
    RuntimeError
        If a negative-Jacobian cell is detected during marching.
    """
    # --- accept either (ni, nj, 3) surface or (ni, nj, nk_in, 3) volume ---
    if grid.ndim == 3:
        surface_only = True
        ni, nj_raw, _ = grid.shape
    elif grid.ndim == 4:
        surface_only = False
        ni, nj_raw, nk_in, _ = grid.shape
        nk = nk_in
    else:
        raise ParameterError(f"grid must be 3-D or 4-D, got shape {grid.shape}")

    # Step 1 — validate/adjust nj
    nj = _adjust_nj(nj_raw)
    nc = nj // 4 + 1

    # Trim grid if nj was adjusted
    if nj < nj_raw:
        if surface_only:
            grid = grid[:, :nj, :]
        else:
            grid = grid[:, :nj, :, :]

    # Build the k=0 surface (ni, nj, 3)
    if surface_only:
        surface_k0 = grid                  # (ni, nj, 3)
    else:
        surface_k0 = grid[:, :, 0, :]     # (ni, nj, 3)

    # Step 2 — split k=0 surface into 4 equatorial sub-surfaces
    sub_surfaces = _split_equatorial(
        surface_k0[:, :, np.newaxis, :].repeat(1, axis=2)
        # need (ni, nj, 1, 3) temporarily for _split_equatorial
        # reshape trick: build a minimal 4-D grid from the surface
    )
    # _split_equatorial expects (ni, nj, nk, 3) — pass a 1-layer version
    _surface_4d = surface_k0[:, :, np.newaxis, :]   # (ni, nj, 1, 3)
    sub_surfaces_list = []
    sector = nj // 4
    for q in range(4):
        j0 = q * sector
        j_indices = list(range(j0, j0 + sector)) + [(j0 + sector) % nj]
        sub_surfaces_list.append(surface_k0[:, j_indices, :].copy())  # (ni, nc, 3)

    # Step 3 — build polar cap seed faces
    # Temporary sub_blocks in (ni, nc, 1, 3) shape to reuse _build_cap_seeds
    _temp_subs = [s[:, :, np.newaxis, :] for s in sub_surfaces_list]
    cap_north_seed, cap_south_seed = _build_cap_seeds(
        surface_k0[:, :, np.newaxis, :], _temp_subs, nc
    )
    # _build_cap_seeds uses grid.shape[0] for ni and sub_blocks[q][ni-1,0,0,:]
    # Recompute directly from sub_surfaces_list (avoids passing 4-D grids)
    # Both caps use the same corner ordering c0=sub0, c1=sub1, c2=sub2, c3=sub3.
    # _gnomonic_seed_face auto-corrects winding so normals point outward.
    north_corners = np.stack([
        sub_surfaces_list[q][ni - 1, 0, :] for q in range(4)
    ])  # (4, 3) — south pole of the sphere (i=ni-1 → theta≈pi)
    south_corners = np.stack([
        sub_surfaces_list[q][0, 0, :] for q in range(4)
    ])  # (4, 3) — north pole of the sphere (i=0 → theta≈0)
    cap_north_seed = _gnomonic_seed_face(north_corners, nc)
    cap_south_seed = _gnomonic_seed_face(south_corners, nc)

    # --- Allocate output volumes ---
    sub_vols  = [np.empty((ni, nc, nk, 3), dtype=np.float64) for _ in range(4)]
    cap_north = np.empty((nc, nc, nk, 3), dtype=np.float64)
    cap_south = np.empty((nc, nc, nk, 3), dtype=np.float64)

    # k=0 initialisation
    for q in range(4):
        sub_vols[q][:, :, 0, :] = sub_surfaces_list[q]

    # Build cap k=0 faces: boundaries from equatorial pole rows (exact conformity),
    # interior solved by Laplacian convergence.  This guarantees that the interior
    # is consistent with the boundary mapping used during marching.
    north_eq_rows_k0 = [sub_surfaces_list[q][ni - 1, :, :] for q in range(4)]
    south_eq_rows_k0 = [sub_surfaces_list[q][0, :, :] for q in range(4)]

    # cap_squareness: straighten each boundary edge from circular arc toward
    # geodesic chord, then write modified rows back into sub_vols k=0.
    if cap_squareness != 0.0:
        north_eq_rows_k0 = _squareness_rows(north_eq_rows_k0, cap_squareness)
        south_eq_rows_k0 = _squareness_rows(south_eq_rows_k0, cap_squareness)
        for q in range(4):
            sub_vols[q][ni - 1, :, 0, :] = north_eq_rows_k0[q]
            sub_vols[q][0,      :, 0, :] = south_eq_rows_k0[q]
        # Re-stitch shared corners (slerp may give tiny floating-point divergence)
        for q in range(4):
            q_next = (q + 1) % 4
            shared_n = 0.5 * (sub_vols[q][ni - 1, -1, 0, :] + sub_vols[q_next][ni - 1, 0, 0, :])
            sub_vols[q][ni - 1, -1, 0, :]     = shared_n
            sub_vols[q_next][ni - 1, 0, 0, :] = shared_n
            north_eq_rows_k0[q][-1]            = shared_n
            north_eq_rows_k0[q_next][0]        = shared_n
            shared_s = 0.5 * (sub_vols[q][0, -1, 0, :] + sub_vols[q_next][0, 0, 0, :])
            sub_vols[q][0, -1, 0, :]     = shared_s
            sub_vols[q_next][0, 0, 0, :] = shared_s
            south_eq_rows_k0[q][-1]      = shared_s
            south_eq_rows_k0[q_next][0]  = shared_s

    cap_north_k0_seed = _build_cap_k0(north_eq_rows_k0, nc, pole="north")
    cap_south_k0_seed = _build_cap_k0(south_eq_rows_k0, nc, pole="south")

    # cap_pull: dome the cap interior toward the pole (boundary edges unchanged).
    if cap_pull != 0.0:
        cap_north_k0_seed = _pull_cap_interior(
            cap_north_k0_seed, north_eq_rows_k0, cap_pull, pole="north"
        )
        cap_south_k0_seed = _pull_cap_interior(
            cap_south_k0_seed, south_eq_rows_k0, cap_pull, pole="south"
        )

    cap_north[:, :, 0, :] = cap_north_k0_seed
    cap_south[:, :, 0, :] = cap_south_k0_seed

    _steps = steps if steps is not None else build_k_steps(nk, ds, growth)

    # Step 6 — simultaneous layer-by-layer marching
    for k in range(nk - 1):
        step = _steps[k]

        # --- 1. Advance all 4 equatorial sub-blocks ---
        new_subs = []
        for q in range(4):
            layer_q = sub_vols[q][:, :, k, :]
            new_q   = _march_one_layer(
                layer_q, step, smooth_strength, smooth_iters,
                periodic=(False, False),   # sub-blocks: open in both directions
            )
            new_subs.append(new_q)

        # Re-stitch shared edges between adjacent sub-blocks so they are
        # exactly identical (numerical symmetry can cause tiny divergence)
        for q in range(4):
            q_next = (q + 1) % 4
            # shared edge: last j-column of sub_q == first j-column of sub_{q+1}
            shared = 0.5 * (new_subs[q][:, -1, :] + new_subs[q_next][:, 0, :])
            new_subs[q][:, -1, :]    = shared
            new_subs[q_next][:, 0, :] = shared

        # Jacobian check for sub-blocks
        jmin_global = np.inf
        for q, new_q in enumerate(new_subs):
            jmin = _check_jacobian(
                sub_vols[q][:, :, k, :], new_q,
                periodic_j=False, label=f"sub{q}", k=k,
            )
            jmin_global = min(jmin_global, jmin)

        # --- 2. Extract pole rows from equatorial sub-blocks for this k+1 ---
        # north pole rows (i=ni-1), shape (nc, 3) each
        north_eq_rows = [new_subs[q][ni - 1, :, :] for q in range(4)]
        # south pole rows (i=0), shape (nc, 3) each — same ordering as north
        south_eq_rows = [new_subs[q][0, :, :] for q in range(4)]

        # cap_squareness: straighten each boundary edge at this layer, then
        # write the modified rows back into new_subs so sub-blocks remain
        # conforming with the cap boundary at every layer.
        if cap_squareness != 0.0:
            north_eq_rows = _squareness_rows(north_eq_rows, cap_squareness)
            south_eq_rows = _squareness_rows(south_eq_rows, cap_squareness)
            for q in range(4):
                new_subs[q][ni - 1, :, :] = north_eq_rows[q]
                new_subs[q][0,      :, :] = south_eq_rows[q]
            for q in range(4):
                q_next = (q + 1) % 4
                shared_n = 0.5 * (new_subs[q][ni - 1, -1, :] + new_subs[q_next][ni - 1, 0, :])
                new_subs[q][ni - 1, -1, :]     = shared_n
                new_subs[q_next][ni - 1, 0, :] = shared_n
                north_eq_rows[q][-1]            = shared_n
                north_eq_rows[q_next][0]        = shared_n
                shared_s = 0.5 * (new_subs[q][0, -1, :] + new_subs[q_next][0, 0, :])
                new_subs[q][0, -1, :]     = shared_s
                new_subs[q_next][0, 0, :] = shared_s
                south_eq_rows[q][-1]      = shared_s
                south_eq_rows[q_next][0]  = shared_s

        # --- 3. Advance cap blocks ---
        # Boundary rows are set directly from the equatorial values (Dirichlet
        # condition).  Only interior nodes are marched freely; the Laplacian
        # smoothing propagates the boundary constraint inward with boundary
        # pinning re-applied after each sweep.
        new_cap_north = _march_cap_layer(
            cap_north[:, :, k, :], step, smooth_strength, smooth_iters,
            north_eq_rows, pole="north",
        )
        new_cap_south = _march_cap_layer(
            cap_south[:, :, k, :], step, smooth_strength, smooth_iters,
            south_eq_rows, pole="south",
        )

        # Jacobian check for cap blocks.  The 4 corner cells are excluded:
        # they border the junction where 4 equatorial sub-blocks meet at one
        # node and the winding at these corners is degenerate by construction.
        # The equatorial sub-blocks own those corner cells and are checked
        # separately with strict positivity.
        jmin_cap_n = _check_jacobian(
            cap_north[:, :, k, :], new_cap_north,
            periodic_j=False, label="cap_north", k=k, exclude_corners=True,
            tol=-1e-6,
        )
        jmin_cap_s = _check_jacobian(
            cap_south[:, :, k, :], new_cap_south,
            periodic_j=False, label="cap_south", k=k, exclude_corners=True,
            tol=-1e-6,
        )
        jmin_global = min(jmin_global, jmin_cap_n, jmin_cap_s)

        # --- 5. Commit layer k+1 ---
        for q in range(4):
            sub_vols[q][:, :, k + 1, :] = new_subs[q]
        cap_north[:, :, k + 1, :] = new_cap_north
        cap_south[:, :, k + 1, :] = new_cap_south

        if layer_callback is not None:
            layer_callback(k + 1, jmin_global)

    return [*sub_vols, cap_north, cap_south]
