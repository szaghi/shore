"""Wing-tip caps of the C-grid (ADR 0059).

The C-grid wraps a wing's sides and leaves its span ends open. A cap closes
one tip as an overset (Chimera) component: a body-fitted grid on the flat
tip face, extruded outwards.

- **Footprint, in the tip plane:** an O-H layout of five blocks. A
  ``core`` sits inside the section, inset by ``core_frac`` of the
  half-thickness from both sides. A ring of four blocks joins it to the
  section's wall: ``upper`` (the C-grid's i = 0 side), ``le`` round the
  leading edge between the two LE cuts, ``lower``, and ``te`` across the
  trailing-edge base (a blunt or rounded TE). The ring's lines run straight
  from the core to the wall nodes, which are the C-grid's own tip-station
  nodes. The core's corners sit on 45 deg diagonals from the LE cuts and TE
  corners, so the three blocks at each corner share it at about 90 / 135 /
  135 deg.
- **Volume:** k leaves the tip plane along its outward normal, at the
  distances ``k_distances`` (the wrap's wall-normal schedule). Ring lines
  lean outwards by up to ``flare_deg``, growing linearly from the core to
  the section's edge, so the cap overhangs the wing's side.
- **Faces:** every block is right-handed. ``k_lo`` is the wall (the tip
  face). The faces between cap blocks are SHARED, with axis maps found from
  their corner nodes. Every other face is FREE: an overset fringe.

The caps are built, not marched; ``check_cap_blocks`` refuses a cap with a
folded cell or a block without cells (ADR 0058).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from shore.exceptions import ParameterError
from shore.volume.spacing import build_i_coords

if TYPE_CHECKING:
    from shore.mesh.block import HexBlock


def check_cap_blocks(blocks: list) -> None:
    """Refuse a wing-tip cap with a folded cell or a block without cells (ADR 0058).

    The caps are built, not marched, so no marcher guard sees them. Every
    block must have at least one cell along each axis and every cell a
    minimum corner scaled Jacobian above 0 (``shore.quality``).

    Raises
    ------
    ParameterError
        Naming the first invalid block and its worst cell.
    """
    from shore.quality import scaled_jacobian

    for b in blocks:
        shape = tuple(int(n) for n in b.nodes.shape[:3])
        if min(shape) < 2:
            raise ParameterError(
                f"Tip cap block {b.label} has no cells (nodes {shape}): the cap "
                "construction cannot close this section; no cap was built (ADR 0058)."
            )
        sj = scaled_jacobian(b.nodes)
        if sj.min() > 0.0:
            continue
        i, j, k = np.unravel_index(int(np.argmin(sj)), sj.shape)
        raise ParameterError(
            f"Folded cell (scaled Jacobian {float(sj.min()):.3f}) in tip cap block "
            f"{b.label} at cell ({i}, {j}, {k}), {int((sj <= 0.0).sum())} of {sj.size} "
            "folded: the cap construction cannot close this section; no cap was "
            "built (ADR 0058)."
        )


#: Distance of a capped C-grid's tip stations from the span ends, as a
#: fraction of the span: the stations sit on the tip faces, where the caps'
#: walls go, instead of 2 % inside them (ADR 0059).
TIP_STATION_MARGIN = 1e-5

#: Cap blocks per tip, in order: the core and the ring round it.
CAP_BLOCKS: tuple[str, ...] = ("core", "upper", "le", "lower", "te")


def section_normal(cross_section: np.ndarray, outward: np.ndarray) -> np.ndarray:
    """Unit normal of the section's plane, on the side of *outward*."""
    pts = cross_section - cross_section.mean(axis=0)
    n = np.linalg.svd(pts, full_matrices=False)[2][-1]
    n = n / np.linalg.norm(n)
    return n if float(n @ outward) > 0.0 else -n


def build_cap_blocks(
    cross_section: np.ndarray,
    *,
    outward: np.ndarray,
    le_cut_upper_idx: int,
    le_cut_lower_idx: int,
    k_distances: np.ndarray,
    n_cap_j: int = 12,
    j_spacing: str = "tanh",
    j_beta: float = 4.0,
    core_frac: float = 0.5,
    flare_deg: float = 30.0,
    tip_inset: float = 0.0,
) -> dict[str, np.ndarray]:
    """The O-H cap of one wing tip (ADR 0059): a core and four ring blocks.

    Parameters
    ----------
    cross_section :
        ``(ni, 3)`` wall nodes of the tip station, from the upper TE
        (i = 0) round the LE to the lower TE, in one plane. The TE must be
        blunt or rounded (``auto_round_sharp_te``): the ``te`` block spans
        its base. These nodes become the cap's edge, bit for bit when
        ``tip_inset`` is 0.
    outward :
        Any vector pointing away from the wing: it picks the side of the
        section's plane the cap grows on (``section_normal``).
    le_cut_upper_idx, le_cut_lower_idx :
        i-indices of the LE cuts, symmetric about the middle node
        (``find_symmetric_le_cuts``); they bound the ``le`` block, and the
        ``te`` block gets the same node count across the TE base.
    k_distances :
        ``(nk,)`` distances of the layers from the wall, from 0 and
        increasing.
    n_cap_j :
        Ring blocks' node count from the core to the wall.
    j_spacing, j_beta :
        Spacing law of that direction (``"uniform"`` / ``"tanh"`` /
        ``"tanh2"``; ``"tanh"`` clusters at the wall) and its strength.
    core_frac :
        The core's half-thickness as a fraction of the section's, in
        (0, 1).
    flare_deg :
        Outward lean of the ring's k-lines at the section's edge, in
        [0, 90) degrees, growing linearly from 0 at the core.
    tip_inset :
        Distance of the wall from the section, along the outward normal.

    Returns
    -------
    dict[str, ndarray]
        ``CAP_BLOCKS`` -> ``(ni, nj, nk, 3)`` right-handed node arrays, the
        wall at ``k = 0``. In the ring blocks j runs from the core to the
        section's wall.

    Raises
    ------
    ParameterError
        On a sharp TE, cuts that are not symmetric or do not straddle the
        LE, a section that is not single-valued along its chord, or
        parameters out of range.
    """
    cs = np.asarray(cross_section, dtype=np.float64)
    ni = cs.shape[0]
    if np.linalg.norm(cs[0] - cs[-1]) <= 1e-12 * np.ptp(cs, axis=0).max():
        raise ParameterError("build_cap_blocks: the section's TE is sharp; round it first.")
    if not (0 < le_cut_upper_idx < le_cut_lower_idx < ni - 1):
        raise ParameterError("build_cap_blocks: LE cuts out of range")
    if le_cut_lower_idx != ni - 1 - le_cut_upper_idx:
        raise ParameterError("build_cap_blocks: the LE cuts must be symmetric in i")
    if not 0.0 < core_frac < 1.0:
        raise ParameterError(f"build_cap_blocks: core_frac must be in (0, 1), got {core_frac}")
    if j_spacing not in ("uniform", "tanh", "tanh2"):
        raise ParameterError(f"build_cap_blocks: unknown j_spacing {j_spacing!r}")
    if n_cap_j < 2:
        raise ParameterError(f"build_cap_blocks: n_cap_j must be >= 2, got {n_cap_j}")
    if not 0.0 <= flare_deg < 90.0:
        raise ParameterError(f"build_cap_blocks: flare_deg must be in [0, 90), got {flare_deg}")
    dist = np.asarray(k_distances, dtype=np.float64)
    if dist.ndim != 1 or dist.size < 2 or dist[0] != 0.0 or np.any(np.diff(dist) <= 0.0):
        raise ParameterError("build_cap_blocks: k_distances must start at 0 and increase")
    n = section_normal(cs, np.asarray(outward, dtype=np.float64))
    wall = cs + tip_inset * n

    # In-plane frame: u along the chord (LE -> TE), v across it.
    te_mid = 0.5 * (wall[0] + wall[-1])
    le = wall[int(np.argmax(np.linalg.norm(wall - te_mid, axis=1)))]
    c = te_mid - le
    c = c - (c @ n) * n
    c = c / np.linalg.norm(c)
    t = np.cross(n, c)
    u, v = (wall - le) @ c, (wall - le) @ t
    up, lo = le_cut_upper_idx, le_cut_lower_idx
    ile = int(np.argmin(u))
    if not up < ile < lo:
        raise ParameterError("build_cap_blocks: the LE cuts must straddle the leading edge")
    us, vs = u[: ile + 1][::-1], v[: ile + 1][::-1]  # one side, LE -> TE
    ul, vl = u[ile:], v[ile:]
    if np.any(np.diff(us) <= 0.0) or np.any(np.diff(ul) <= 0.0):
        raise ParameterError("build_cap_blocks: the section is not single-valued along its chord")

    def mid(x: np.ndarray) -> np.ndarray:
        return 0.5 * (np.interp(x, us, vs) + np.interp(x, ul, vl))

    def half(x: np.ndarray) -> np.ndarray:
        return 0.5 * np.abs(np.interp(x, us, vs) - np.interp(x, ul, vl))

    def point(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return le + np.multiply.outer(x, c) + np.multiply.outer(y, t)

    def fractions(arc: np.ndarray) -> np.ndarray:
        d = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(arc, axis=0), axis=1))]
        return d / d[-1]

    # Wall arcs of the four ring blocks.
    side_a = wall[: up + 1]  # TE -> LE cut, the i = 0 side
    side_b = wall[lo:][::-1]  # TE -> LE cut, the other side
    nose = wall[up : lo + 1]  # round the LE
    nb = nose.shape[0]
    s = np.linspace(0.0, 1.0, nb)[1:-1, None]
    base = np.vstack([wall[-1], (1.0 - s) * wall[-1] + s * wall[0], wall[0]])  # across the TE

    # The core: between the two sides, inset by core_frac of the half-thickness.
    x_le = 0.5 * (u[up] + u[lo])
    x_te = 0.5 * (u[0] + u[-1])
    x_a = x_le + (1.0 - core_frac) * half(np.array(x_le))
    x_b = x_te - (1.0 - core_frac) * half(np.array(x_te))
    if not x_a < x_b:
        raise ParameterError("build_cap_blocks: the section is too short for the core")
    sign = 1.0 if float(np.mean(v[:ile])) > float(np.mean(v[ile + 1 :])) else -1.0
    xa_side = x_b + (x_a - x_b) * fractions(side_a)
    xb_side = x_b + (x_a - x_b) * fractions(side_b)
    core_a = point(xa_side, mid(xa_side) + sign * core_frac * half(xa_side))  # TE -> LE
    core_b = point(xb_side, mid(xb_side) - sign * core_frac * half(xb_side))
    qa_te, qa_le, qb_le, qb_te = core_a[0], core_a[-1], core_b[-1], core_b[0]
    core_nose = _segment(qa_le, qb_le, fractions(nose))  # side a -> side b at the LE
    core_base = _segment(qb_te, qa_te, fractions(base))  # side b -> side a at the TE

    # Core: i from the TE to the LE, j from side b to side a (Coons patch).
    core = _coons(core_b, core_a, core_base, core_nose[::-1])

    raw = build_i_coords(n_cap_j, 0.0, 1.0, spacing=j_spacing, beta=j_beta)
    tj = 1.0 - raw[::-1] if j_spacing == "tanh" else raw
    footprints = {
        "core": (core, np.zeros(core.shape[:2]), np.zeros(core.shape)),
        "upper": _ring(core_a, side_a, tj),
        "le": _ring(core_nose, nose, tj),
        "lower": _ring(core_b, side_b, tj),
        "te": _ring(core_base, base, tj),
    }
    tan_f = np.tan(np.radians(flare_deg))
    blocks: dict[str, np.ndarray] = {}
    for name in CAP_BLOCKS:
        fp, r, m = footprints[name]
        d = n + (tan_f * r)[..., None] * m
        d = d / np.linalg.norm(d, axis=-1, keepdims=True)
        x = fp[:, :, None, :] + dist[None, None, :, None] * d[:, :, None, :]
        e_i, e_j = x[1, 0, 0] - x[0, 0, 0], x[0, 1, 0] - x[0, 0, 0]
        blocks[name] = x if float(np.cross(e_i, e_j) @ n) > 0.0 else x[::-1].copy()
    return blocks


def _segment(a: np.ndarray, b: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Points from *a* to *b* at fractions *f*, the ends exact."""
    out = (1.0 - f)[:, None] * a + f[:, None] * b
    out[0], out[-1] = a, b
    return out


def _coons(bottom: np.ndarray, top: np.ndarray, left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Coons patch: rows bottom/top along i, columns left/right along j."""
    ni, nj = bottom.shape[0], left.shape[0]
    xi = np.linspace(0.0, 1.0, ni)[:, None, None]
    et = np.linspace(0.0, 1.0, nj)[None, :, None]
    p = (
        (1 - et) * bottom[:, None]
        + et * top[:, None]
        + (1 - xi) * left[None]
        + xi * right[None]
        - (1 - xi) * (1 - et) * bottom[0]
        - (1 - xi) * et * top[0]
        - xi * (1 - et) * bottom[-1]
        - xi * et * top[-1]
    )
    p[:, 0], p[:, -1], p[0], p[-1] = bottom, top, left, right
    return p


def _ring(
    inner: np.ndarray, wall: np.ndarray, tj: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A ring block's footprint, straight lines from the core to the wall.

    Returns the footprint, each node's fraction of the way to the wall and
    its line's in-plane direction: the flare's weight and direction.
    """
    fp = (1.0 - tj)[None, :, None] * inner[:, None] + tj[None, :, None] * wall[:, None]
    fp[:, 0], fp[:, -1] = inner, wall
    m = wall - inner
    m = m / np.linalg.norm(m, axis=-1, keepdims=True)
    return fp, np.broadcast_to(tj, fp.shape[:2]), np.broadcast_to(m[:, None], fp.shape)


def wire_cap_blocks(blocks: list[HexBlock]) -> None:
    """Make the faces two cap blocks share SHARED, with their axis maps.

    The i and j faces only (``shore.mesh.seams.wire_by_corners``): the k
    faces are the wall and the outer fringe.
    """
    from shore.mesh.seams import wire_by_corners

    wire_by_corners(blocks, ("i_lo", "i_hi", "j_lo", "j_hi"))


def find_symmetric_le_cuts(
    cross_section: np.ndarray,
    le_frac: float,
    span_axis_world: int,
) -> tuple[int, int]:
    """Find the upper / lower LE-cut indices in a body cross-section.

    The cross-section sweeps upper-TE-cut → suction → LE → pressure →
    lower-TE-cut.  This locates the i-indices on each side of the LE
    point that correspond to chord position ``le_frac``, snapped to be
    **symmetric** around the LE point in i-station count (so
    ``upper_idx + lower_idx == 2 * le_idx``).  Symmetry is required by
    the cap-fill builder for the chord-cut SHARED face.

    Parameters
    ----------
    cross_section :
        ``(ni, 3)`` polyline.
    le_frac :
        Chord-fraction of the LE-cut, in ``(0, 1)``.
    span_axis_world :
        Which world axis (0 / 1 / 2) is the span direction; the
        chord axis is the in-plane axis with the largest extent.

    Returns
    -------
    (le_cut_upper_idx, le_cut_lower_idx) : tuple of int

    Raises
    ------
    ParameterError
        On invalid ``le_frac`` or when the i-axis is too coarse to
        resolve the requested cut.
    """
    if not (0.0 < le_frac < 1.0):
        raise ParameterError(f"find_symmetric_le_cuts: le_frac must be in (0, 1), got {le_frac}")
    if span_axis_world not in (0, 1, 2):
        raise ParameterError(
            f"find_symmetric_le_cuts: span_axis_world must be 0/1/2, got {span_axis_world}"
        )

    ni = cross_section.shape[0]
    extents = cross_section.max(axis=0) - cross_section.min(axis=0)
    extents[span_axis_world] = -np.inf
    chord_axis = int(np.argmax(extents))

    x = cross_section[:, chord_axis]
    x_le = float(x.min())
    x_te = float(x.max())
    chord_len = x_te - x_le
    le_idx = int(np.argmin(x))
    x_target = x_le + le_frac * chord_len

    upper_idx = le_idx
    while upper_idx > 0 and x[upper_idx] < x_target:
        upper_idx -= 1
    # Symmetry around the body's i-axis midpoint (which corresponds to
    # the LE for a typical airfoil i-distribution).  This guarantees
    # equal i-station counts for cap_upper and cap_lower regardless of
    # whether ``ni_body`` is even or odd.
    lower_idx = ni - 1 - upper_idx

    if lower_idx <= upper_idx:
        raise ParameterError(
            f"find_symmetric_le_cuts: chord fraction {le_frac} resolves to "
            f"a single i-station ({upper_idx}); pick a larger fraction or "
            f"use a finer i-axis."
        )
    if upper_idx == le_idx or lower_idx == le_idx:
        raise ParameterError(
            f"find_symmetric_le_cuts: chord fraction {le_frac} resolves to "
            f"the LE point itself; pick a larger fraction or use a finer i-axis."
        )
    return upper_idx, lower_idx


def detect_te_thickness(
    cross_section: np.ndarray,
    span_axis_world: int,
) -> float:
    """Return the TE thickness of a cross-section polyline.

    Sharp TE → ``0.0`` (i=0 and i=ni-1 coincide).  Blunt TE → distance
    between i=0 and i=ni-1 in the airfoil's thickness direction.
    """
    delta = cross_section[-1, :] - cross_section[0, :]
    in_plane = delta.copy()
    in_plane[span_axis_world] = 0.0
    return float(np.linalg.norm(in_plane))


def auto_round_sharp_te(
    cross_section: np.ndarray,
    te_thickness: float,
    span_axis_world: int,
) -> np.ndarray:
    """Synthesise a blunt cross-section by splitting the sharp TE node.

    Returns a copy of ``cross_section`` with the i=0 and i=ni-1 nodes
    moved apart by ``te_thickness`` along the airfoil-thickness
    direction (perpendicular to chord, in the airfoil's 2D plane), so
    the TE-base segment has finite extent.  Used by ``CGridTipCap`` to
    auto-round sharp-TE inputs *only inside the cap*; the user's
    body wrap is left untouched.

    The thickness direction is the in-plane axis with the smaller
    spread (the same heuristic used by ``_detect_thickness_axis``).
    The split is anti-symmetric: i=0 moves by +te_thickness/2 and
    i=ni-1 by -te_thickness/2 so the TE-base midpoint stays at the
    original sharp-TE location.

    Parameters
    ----------
    cross_section :
        ``(ni, 3)`` polyline.  Must have ``cross_section[0]`` and
        ``cross_section[-1]`` coincident (sharp TE).
    te_thickness :
        Desired TE thickness.  Must be > 0.
    span_axis_world :
        World axis (0 / 1 / 2) that is the spanwise direction; used
        to identify the airfoil-thickness in-plane axis.

    Returns
    -------
    ndarray, shape ``(ni, 3)``, same dtype as input.

    Raises
    ------
    ParameterError
        On invalid inputs or when the cross-section is not actually
        sharp.
    """
    if te_thickness <= 0:
        raise ParameterError(f"auto_round_sharp_te: te_thickness must be > 0, got {te_thickness}")
    if span_axis_world not in (0, 1, 2):
        raise ParameterError(
            f"auto_round_sharp_te: span_axis_world must be 0/1/2, got {span_axis_world}"
        )
    delta = cross_section[-1, :] - cross_section[0, :]
    if float(np.linalg.norm(delta)) > 1e-9:
        raise ParameterError(
            f"auto_round_sharp_te: cross-section is not sharp (TE endpoints "
            f"differ by {float(np.linalg.norm(delta)):.6g}); auto-round "
            f"only operates on sharp TE inputs."
        )

    extents = cross_section.max(axis=0) - cross_section.min(axis=0)
    extents[span_axis_world] = -np.inf
    chord_axis = int(np.argmax(extents))
    candidates = [a for a in (0, 1, 2) if a != span_axis_world and a != chord_axis]
    thickness_axis = candidates[0]

    out = cross_section.copy()
    half = 0.5 * te_thickness
    out[0, thickness_axis] += +half
    out[-1, thickness_axis] += -half
    return out
