"""The collar block of two intersecting bodies (ADR 0010).

``build_collar`` wraps one Chimera volume block around the junction of two
bodies:

1. **Wall.** ``intersection_loop`` (ADR 0009) gives the junction curve and
   ``march_surface`` (ADR 0008) grows a surface grid on each body away from
   the other; the two halves meet at the junction row ``i = len(steps_b)``.
   j runs around the loop (periodic, seam node not repeated), i from B's far
   edge to A's far edge, k leaves the walls.
2. **Near-wall layers.** The hyperbolic scheme (ADR 0005) while the layers
   are thinner than half the junction spacing *h*. Beyond that, orthogonal
   offsetting collapses the cells beside the concave crease (at a 90 deg
   junction the cell next to it has width ``h - d`` at thickness *d*).
3. **Outer boundary.** Per loop station, the near-wall i-line is offset by
   the remaining thickness in its best-fit plane. The swallowtail that the
   concave junction makes is trimmed at its self-intersection, and the
   nodes are placed by arc length with the trimmed crease mapped to the
   junction row.
4. **Interior.** Per loop station, TFI between the near-wall line and the
   outer boundary with the k schedule, then Winslow smoothing with
   Thomas-Middlecoff control functions (the clustering of the boundaries is
   kept, as in ADR 0007), solved by Picard iteration with SuperLU.

Steps 2-4 are ``extrude_collar``. With ``plane`` (ADR 0013) the collar is
the half-collar of a half model: the open curve of ``intersection_curve``
(ADR 0011), both halves marched with their ends on the plane (ADR 0012),
and the volume extruded from the wall mirrored across the plane into a
closed loop, of which the first ``n_loop`` stations are kept; the j-faces
are put on the plane and marked ``FaceBC.SYMMETRY``. ``loop_spacing`` and
``loop_curvature`` distribute the nodes along the junction curve (ADR 0020).

``split_collar`` (ADR 0024) is the collar of a loop with sharp corners
(turning more than ``CollarParams.sharp_deg``: a sharp trailing edge), where
the surface march around the corner folds. The loop is cut at each sharp
corner into pieces. Each piece is extended at both ends by a wake line on A
along the corner's bisector, and the wall over the wake is a virtual sheet
B': the corner's B-line translated along the wake. The two sides of every
sheet are joined by a SHARED seam (a C-grid wake cut). Each piece is
extruded as one open block and split into six (wake, body, wake by B, A) so
that every face has one boundary condition. Loops without sharp corners
give ``build_collar``'s block.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from shore.exceptions import ParameterError
from shore.geometry.base import BODY, Geometry
from shore.mesh.block import HexBlock
from shore.mesh.face import AxisMap, FaceBC
from shore.mesh.spacing import Spacing1D
from shore.split import BlockSplit, split_blocks
from shore.surface.intersect import DEFAULT as TRACE
from shore.surface.intersect import _loop, intersection_curve, intersection_loop
from shore.surface.march import _Ends, _march_between, _Mirror, march_surface
from shore.surface.tfi import _control
from shore.volume.hyperbolic_scheme import hyperbolic_layer


@dataclass(frozen=True)
class CollarParams:
    """Collar settings (not per-case tuning knobs)."""

    hyperbolic_fraction: float = 0.5  # march hyperbolically while thinner than this x h
    max_iterations: int = 400
    tolerance: float = 1e-12  # max node update, relative to the collar thickness
    control_limit: float = 1.5  # |P|, |Q| bound: keeps the 5-point operator an M-matrix
    oblique_deg: float = 1.0  # warn when an open curve's end is this far off the plane normal
    sharp_deg: float = 105.0  # split_collar: corners turning more than this are cut (ADR 0024)
    wake_ratio: float = 1.2  # split_collar: largest ratio of adjacent wake cells


DEFAULT = CollarParams()


def build_collar(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    *,
    n_loop: int,
    steps_a: np.ndarray,
    steps_b: np.ndarray,
    steps_k: np.ndarray,
    label: str = "collar",
    params: CollarParams = DEFAULT,
    plane: tuple[np.ndarray, np.ndarray] | None = None,
    loop_spacing: Spacing1D | None = None,
    loop_curvature: float = 0.0,
) -> HexBlock:
    """The collar block around the junction of *body_a* and *body_b*.

    Parameters
    ----------
    body_a, body_b :
        The two bodies (``Geometry.project`` with normals).
    seed : np.ndarray, shape (3,)
        A point near the junction; selects the loop when there are several.
    n_loop :
        Nodes around the loop (j), >= 3.
    steps_a, steps_b :
        Surface marching steps away from the junction on A and on B.
    steps_k :
        Wall-normal steps; the collar is ``sum(steps_k)`` thick.
    label :
        Block label.
    params :
        Collar settings.
    plane :
        ``(point, normal)`` of a symmetry plane: build the half-collar of a
        half model, on the seed's side. ``None``: the closed collar.
    loop_spacing, loop_curvature :
        Node distribution along the junction curve, passed to
        ``intersection_loop`` / ``intersection_curve`` as *spacing* and
        *curvature* (ADR 0020; default: equally spaced in arc length).

    Returns
    -------
    HexBlock
        Nodes ``(len(steps_b) + len(steps_a) + 1, n_loop, len(steps_k) + 1, 3)``;
        face ``k_lo`` is WALL. Closed collar: ``periodic_j``, the other faces
        FREE. Half-collar: open in j, faces ``j_lo`` / ``j_hi`` on the plane
        and SYMMETRY, the others FREE.

    Raises
    ------
    ParameterError
        Malformed steps or plane, ``n_loop < 3``, or no intersection near
        *seed*.

    Warns
    -----
    UserWarning
        The half-collar's curve meets the plane obliquely (its end tangent
        more than ``params.oblique_deg`` off the plane normal): the bodies
        are not mirror-symmetric there, and the block's quality is not
        guaranteed.
    """
    if isinstance(n_loop, bool) or not isinstance(n_loop, (int, np.integer)) or n_loop < 3:
        raise ParameterError(f"n_loop must be an integer >= 3, got {n_loop!r}.")
    steps_a, steps_b, steps_k = (
        _steps(name, s)
        for name, s in (("steps_a", steps_a), ("steps_b", steps_b), ("steps_k", steps_k))
    )

    if plane is None:
        curve = intersection_loop(
            body_a, body_b, seed, int(n_loop), spacing=loop_spacing, curvature=loop_curvature
        )
        return _closed_collar(body_a, body_b, curve, steps_a, steps_b, steps_k, label, params)
    mirror = _Mirror.of(plane)
    curve = intersection_curve(
        body_a,
        body_b,
        seed,
        int(n_loop),
        plane=plane,
        spacing=loop_spacing,
        curvature=loop_curvature,
    )
    if _end_angle_deg(curve, mirror.normal) > params.oblique_deg:
        warnings.warn(
            "the junction curve meets the symmetry plane obliquely: the bodies are not "
            "mirror-symmetric there, and the half-collar's quality is not guaranteed.",
            UserWarning,
            stacklevel=2,
        )
    on_a = march_surface(body_a, curve, steps_a, plane=plane)
    on_b = march_surface(body_b, curve[::-1], steps_b, plane=plane)[::-1]
    wall = _wall(on_a, on_b)
    closed = np.concatenate([wall, mirror.reflect(wall[:, -2:0:-1])], axis=1)
    nodes = extrude_collar(closed, steps_k, len(steps_b), params=params)[:, : int(n_loop)].copy()
    nodes[:, [0, -1]] = mirror.onto(nodes[:, [0, -1]])
    block = HexBlock(label, nodes, periodic_i=False, periodic_j=False)
    block._set_face_bc("j_lo", FaceBC.SYMMETRY)
    block._set_face_bc("j_hi", FaceBC.SYMMETRY)
    return block


def split_collar(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    *,
    n_loop: int,
    steps_a: np.ndarray,
    steps_b: np.ndarray,
    steps_k: np.ndarray,
    label: str = "collar",
    params: CollarParams = DEFAULT,
    loop_spacing: Spacing1D | None = None,
    loop_curvature: float = 0.0,
    wake_length: float | None = None,
) -> list[HexBlock]:
    """The collar of a junction whose loop has sharp corners, cut there (ADR 0024).

    Parameters
    ----------
    body_a, body_b, seed, n_loop, steps_a, steps_b, steps_k, label, params, loop_spacing, \
loop_curvature :
        As for ``build_collar`` (closed loops only).
    wake_length :
        Length of the wake lines on A behind each sharp corner (default
        ``sum(steps_a)``, the width of the collar on A).

    Returns
    -------
    list[HexBlock]
        Without sharp corners, ``[build_collar(...)]``. Otherwise six blocks
        per piece *p* of the loop between consecutive sharp corners,
        ``f"{label}_{p}_{part}_{side}"`` with *part* ``wake_lo`` / ``body``
        / ``wake_hi`` (j: the wake behind the piece's last corner, the
        piece, the wake behind its first corner) and *side* ``b`` / ``a``
        (i below / above the junction row). ``k_lo`` is WALL, except on the
        ``wake_*_b`` blocks: the sheet behind a corner, SHARED with the
        other side's sheet block (along i same, along j reversed). The
        seams inside a piece are SHARED; the other faces FREE.

    Raises
    ------
    ParameterError
        As ``build_collar``; a bad *wake_length*; a sharp corner that is not
        an edge of B (on A's edge the wake has no surface to run on).
    """
    if isinstance(n_loop, bool) or not isinstance(n_loop, (int, np.integer)) or n_loop < 3:
        raise ParameterError(f"n_loop must be an integer >= 3, got {n_loop!r}.")
    steps_a, steps_b, steps_k = (
        _steps(name, s)
        for name, s in (("steps_a", steps_a), ("steps_b", steps_b), ("steps_k", steps_k))
    )
    length = float(np.sum(steps_a)) if wake_length is None else wake_length
    if not (isinstance(length, (int, float)) and np.isfinite(length) and length > 0.0):
        raise ParameterError(f"wake_length must be a positive finite length, got {wake_length!r}.")

    loop, corners = _loop(
        body_a, body_b, seed, int(n_loop), loop_spacing, loop_curvature, BODY, BODY, TRACE
    )
    n = len(loop)
    ahead = np.roll(loop, -1, axis=0) - loop
    t = ahead / np.linalg.norm(ahead, axis=1, keepdims=True)
    t_in = np.roll(t, 1, axis=0)
    turn = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", t_in, t), -1.0, 1.0)))
    sharp = [c for c in corners if turn[c] > params.sharp_deg]
    if not sharp:
        return [_closed_collar(body_a, body_b, loop, steps_a, steps_b, steps_k, label, params)]

    wakes, planes = {}, {}
    for c in sharp:
        d = t_in[c] - t[c]
        d /= np.linalg.norm(d)
        first_cell = 0.5 * float(np.linalg.norm(ahead[c - 1]) + np.linalg.norm(ahead[c]))
        cells = _wake_cells(first_cell, length, params.wake_ratio)
        wakes[c] = np.asarray(body_a.project(loop[c] + np.cumsum(cells)[:, None] * d).points)
        planes[c] = _sheet_plane(body_b, loop, c, d)

    ends = list(zip(sharp, [*sharp[1:], sharp[0]], strict=True))
    marched = []
    for first, last in ends:
        piece = loop[np.arange(first, last + (n if last <= first else 0) + 1) % n]
        between = _Ends(_Mirror.of(planes[last]), _Mirror.of(planes[first]))
        marched.append((piece, _march_between(body_b, piece[::-1], steps_b, between)[::-1]))
    # one B-line per corner, shared by the pieces that end there
    lines = {c: [] for c in sharp}
    for (first, last), (_, on_b) in zip(ends, marched, strict=True):
        lines[first].append(on_b[0])
        lines[last].append(on_b[-1])
    line = {c: np.mean(v, axis=0) for c, v in lines.items()}
    for c in sharp:
        line[c][0] = loop[c]

    junction = len(steps_b)
    blocks: list[HexBlock] = []
    sheets: dict[int, list[HexBlock]] = {c: [] for c in sharp}
    for p, ((first, last), (piece, on_b)) in enumerate(zip(ends, marched, strict=True)):
        on_b = on_b.copy()
        on_b[0], on_b[-1] = line[first], line[last]
        front = np.vstack([wakes[first][::-1], piece, wakes[last]])
        on_a = march_surface(body_a, front, steps_a)
        on_b = np.concatenate([
            line[first][None] + (wakes[first][::-1] - loop[first])[:, None],
            on_b,
            line[last][None] + (wakes[last] - loop[last])[:, None],
        ])  # fmt: skip
        wall = _wall(on_a, on_b)  # j reversed: the wake behind `last` comes first
        h = _junction_spacing(wall, junction)
        nodes = _extrude(wall, steps_k, h, junction, params, periodic=False)
        lo, hi = len(wakes[last]), wall.shape[1] - 1 - len(wakes[first])
        whole = HexBlock(f"{label}_{p}", nodes, periodic_i=False, periodic_j=False)
        cuts = [BlockSplit(whole.label, "i", [junction]), BlockSplit(whole.label, "j", [lo, hi])]
        cut = split_blocks([whole], cuts)
        named = {}
        for chunk in cut:
            side = "b" if "__i0" in chunk.label else "a"
            part = ("wake_lo", "body", "wake_hi")[int(chunk.label.split("__j")[1][0])]
            chunk.label = f"{label}_{p}_{part}_{side}"
            named[(part, side)] = chunk
        sheets[last].append(named[("wake_lo", "b")])
        sheets[first].append(named[("wake_hi", "b")])
        blocks += [named[(part, side)] for side in "ba" for part in ("wake_lo", "body", "wake_hi")]
    flip_j = AxisMap.from_pairs(("i", False), ("j", True))
    for c in sharp:
        lo_side, hi_side = sheets[c]
        lo_side._set_face_bc("k_lo", FaceBC.SHARED, hi_side.face("k_lo"), flip_j)
        hi_side._set_face_bc("k_lo", FaceBC.SHARED, lo_side.face("k_lo"), flip_j)
    return blocks


def _closed_collar(
    body_a: Geometry,
    body_b: Geometry,
    curve: np.ndarray,
    steps_a: np.ndarray,
    steps_b: np.ndarray,
    steps_k: np.ndarray,
    label: str,
    params: CollarParams,
) -> HexBlock:
    """``build_collar``'s closed collar over the loop *curve*."""
    on_a = march_surface(body_a, curve, steps_a, periodic=True)
    on_b = march_surface(body_b, curve[::-1], steps_b, periodic=True)[::-1]
    nodes = extrude_collar(_wall(on_a, on_b), steps_k, len(steps_b), params=params)
    return HexBlock(label, nodes, periodic_i=False, periodic_j=True)


def _wall(on_a: np.ndarray, on_b: np.ndarray) -> np.ndarray:
    """The collar wall (ni, nj, 3) from the two surface grids: i from B's far edge to
    A's, the junction row at ``i = m_b``, j reversed (the loop runs along -j)."""
    wall = np.concatenate([on_b[:, ::-1], on_a[:, 1:]], axis=1)
    return np.ascontiguousarray(np.transpose(wall, (1, 0, 2))[:, ::-1])


def _wake_cells(first: float, length: float, ratio: float) -> np.ndarray:
    """Cells from *first*, growing by at most *ratio*, that add up to *length*."""
    if first >= length:
        return np.array([length])
    count = int(np.ceil(np.log1p(length * (ratio - 1.0) / first) / np.log(ratio)))
    lo, hi = 1.0, ratio
    for _ in range(200):  # the growth that fits exactly: total is increasing in it
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if first * np.sum(mid ** np.arange(count)) < length else (lo, mid)
    return first * (0.5 * (lo + hi)) ** np.arange(count)


def _sheet_plane(
    body_b: Geometry, loop: np.ndarray, c: int, d: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The plane of the sheet behind sharp corner *c*: through B's edge and along *d*."""
    mids = 0.5 * (loop[[c - 1, (c + 1) % len(loop)]] + loop[c])
    n1, n2 = np.asarray(body_b.project(mids).normals)
    edge = np.cross(n1, n2)
    if np.linalg.norm(edge) < 0.1:
        raise ParameterError(
            f"the sharp corner at {loop[c].tolist()} is not an edge of body B: "
            "split_collar cuts sharp edges of B only."
        )
    normal = np.cross(d, edge)
    return loop[c], normal / np.linalg.norm(normal)


def _junction_spacing(wall: np.ndarray, junction: int) -> float:
    """The smallest i-spacing at the junction row (the collar's *h*)."""
    return float(
        min(
            np.linalg.norm(wall[junction + 1] - wall[junction], axis=-1).min(),
            np.linalg.norm(wall[junction - 1] - wall[junction], axis=-1).min(),
        )
    )


def extrude_collar(
    wall: np.ndarray, steps_k: np.ndarray, junction: int, *, params: CollarParams = DEFAULT
) -> np.ndarray:
    """The collar volume over a j-periodic *wall* (steps 2-4 of the method).

    Parameters
    ----------
    wall : np.ndarray, shape (ni, nj, 3)
        The collar wall: i from B's far edge through the junction row to A's
        far edge, j around the loop (periodic, seam node not repeated).
    steps_k :
        Wall-normal steps.
    junction :
        The junction row (``len(steps_b)``), 0 < junction < ni - 1.
    params :
        Collar settings.

    Returns
    -------
    np.ndarray, shape (ni, nj, len(steps_k) + 1, 3)
        ``[:, :, 0]`` is *wall*.
    """
    w = np.ascontiguousarray(np.asarray(wall, dtype=np.float64))
    if w.ndim != 3 or w.shape[2] != 3 or w.shape[1] < 3:
        raise ParameterError(f"wall must have shape (ni, nj >= 3, 3), got {w.shape}.")
    if not 0 < junction < w.shape[0] - 1:
        raise ParameterError(f"junction must be an inner i row, got {junction}.")
    steps_k = _steps("steps_k", steps_k)
    return _extrude(w, steps_k, _junction_spacing(w, junction), junction, params)


def _end_angle_deg(curve: np.ndarray, normal: np.ndarray) -> float:
    """Largest angle between *normal* and the curve's second-order end tangents."""
    worst = 0.0
    for r0, r1, r2 in ((curve[0], curve[1], curve[2]), (curve[-1], curve[-2], curve[-3])):
        t = -3.0 * r0 + 4.0 * r1 - r2
        cos = abs(float(t @ normal)) / float(np.linalg.norm(t))
        worst = max(worst, float(np.degrees(np.arccos(min(1.0, cos)))))
    return worst


def _steps(name: str, steps: np.ndarray) -> np.ndarray:
    s = np.asarray(steps, dtype=np.float64)
    if s.ndim != 1 or s.size == 0 or not np.all(np.isfinite(s)) or np.any(s <= 0.0):
        raise ParameterError(f"{name} must be a non-empty 1-D array of positive finite steps.")
    return s


def _extrude(
    wall: np.ndarray,
    steps_k: np.ndarray,
    h: float,
    junction: int,
    params: CollarParams,
    *,
    periodic: bool = True,
) -> np.ndarray:
    ni, nj, _ = wall.shape
    nk = len(steps_k)
    cum = np.r_[0.0, np.cumsum(steps_k)]
    k0 = max(1, int(np.searchsorted(cum, params.hyperbolic_fraction * h)) - 1)
    k0 = min(k0, nk)

    nodes = np.empty((ni, nj, nk + 1, 3))
    nodes[:, :, 0] = wall
    for k in range(k0):
        nodes[:, :, k + 1] = hyperbolic_layer(nodes[:, :, k], float(steps_k[k]), (False, periodic))
    if k0 == nk:
        return nodes

    inner = nodes[:, :, k0]
    normals = _normals(inner, periodic)
    rest = cum[-1] - cum[k0]
    outer = np.stack(
        [_outer_line(inner[:, j], normals[:, j], rest, junction) for j in range(nj)], axis=1
    )
    w = (cum[k0:] - cum[k0]) / rest
    sheets = _tfi(inner, outer, w)  # (ni, nj, nk - k0 + 1, 3)
    if sheets.shape[2] > 2:
        sheets = _winslow(np.transpose(sheets, (1, 0, 2, 3)), cum[-1], params)
        sheets = np.transpose(sheets, (1, 0, 2, 3))
    nodes[:, :, k0:] = sheets
    return nodes


def _normals(layer: np.ndarray, periodic: bool = True) -> np.ndarray:
    """Unit normals of a layer, j-periodic or not (one-sided at the open ends)."""
    ti = np.gradient(layer, axis=0)
    if periodic:
        tj = 0.5 * (np.roll(layer, -1, axis=1) - np.roll(layer, 1, axis=1))
    else:
        tj = np.gradient(layer, axis=1)
    n = np.cross(ti, tj)
    return n / np.linalg.norm(n, axis=-1, keepdims=True)


def _arc(line: np.ndarray) -> np.ndarray:
    """Cumulative arc length of a polyline (not normalised)."""
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))]


def _outer_line(inner: np.ndarray, normals: np.ndarray, d: float, junction: int) -> np.ndarray:
    """The outer boundary of one loop station: offset, trimmed, crease on the junction row."""
    offset = inner + d * normals
    centre = offset.mean(axis=0)
    _, _, frame = np.linalg.svd(np.vstack([inner, offset]) - centre, full_matrices=False)
    local = (offset - centre) @ frame.T  # in-plane (x, y), small out-of-plane z

    poly, crease = _trim(local[:, :2], junction)
    s = _arc(poly)
    s_inner = _arc(inner)
    frac_b = s_inner[: junction + 1] / s_inner[junction]
    frac_a = (s_inner[junction:] - s_inner[junction]) / (s_inner[-1] - s_inner[junction])
    target = np.r_[s[crease] * frac_b, s[crease] + (s[-1] - s[crease]) * frac_a[1:]]
    xy = np.stack([np.interp(target, s, poly[:, c]) for c in range(2)], axis=-1)

    # the out-of-plane part follows the same B / A split of the original indices
    n = len(inner)
    index = np.r_[junction * frac_b, junction + (n - 1 - junction) * frac_a[1:]]
    z = np.interp(index, np.arange(n), local[:, 2])
    return centre + np.c_[xy, z] @ frame


def _trim(poly: np.ndarray, junction: int) -> tuple[np.ndarray, int]:
    """Cut the swallowtails out of a 2-D polyline; also the index of the junction's image.

    The earliest segment that crosses a later one is joined to the latest
    such segment at their crossing, until no two segments cross. A junction
    vertex inside a cut loop maps to the crossing point that replaced it.
    """
    pts = poly.copy()
    origin = np.arange(len(pts))  # original vertex index, -1 for a crossing point
    while True:
        hit = _first_crossing(pts)
        if hit is None:
            break
        a, b, t = hit
        x = pts[a] + t * (pts[a + 1] - pts[a])
        pts = np.vstack([pts[: a + 1], x[None], pts[b + 1 :]])
        origin = np.r_[origin[: a + 1], -1, origin[b + 1 :]]
    hits = np.flatnonzero(origin == junction)
    if hits.size:
        return pts, int(hits[0])
    before = int(np.flatnonzero((origin >= 0) & (origin < junction)).max())
    return pts, before + int(np.flatnonzero(origin[before:] == -1)[0])


def _first_crossing(pts: np.ndarray) -> tuple[int, int, float] | None:
    """(a, b, t) of the crossing with the smallest a, then the largest b > a + 1."""
    p, r = pts[:-1], pts[1:] - pts[:-1]
    m = len(p)
    qp = p[None, :, :] - p[:, None, :]  # q_b - p_a
    den = r[:, None, 0] * r[None, :, 1] - r[:, None, 1] * r[None, :, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (qp[..., 0] * r[None, :, 1] - qp[..., 1] * r[None, :, 0]) / den
        u = (qp[..., 0] * r[:, None, 1] - qp[..., 1] * r[:, None, 0]) / den
    a_idx, b_idx = np.indices((m, m))
    cross = (
        (b_idx >= a_idx + 2)
        & (np.abs(den) >= 1e-300)
        & (t > 0.0)
        & (t < 1.0)
        & (u > 0.0)
        & (u < 1.0)
    )
    if not cross.any():
        return None
    a = int(np.flatnonzero(cross.any(axis=1))[0])
    b = int(np.flatnonzero(cross[a])[-1])
    return a, b, float(t[a, b])


def _tfi(inner: np.ndarray, outer: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Per-station TFI between *inner* and *outer* (ni, nj, 3), blending weights *w* in k.

    The far i-edges are straight lines with the same k blending, so the TFI
    edge corrections vanish and only the k interpolation is left.
    """
    return inner[:, :, None] + w[None, None, :, None] * (outer - inner)[:, :, None]


def _winslow(x: np.ndarray, scale: float, params: CollarParams) -> np.ndarray:
    """Winslow + Thomas-Middlecoff on a stack of 2-D sheets (nj, ni, nk, 3), boundaries fixed.

    ``alpha (x_ii + P x_i) - 2 beta x_ik + gamma (x_kk + Q x_k) = 0`` with
    ``alpha = |x_k|^2``, ``beta = x_i . x_k``, ``gamma = |x_i|^2``; Picard
    iteration, one SuperLU factorisation of the block-diagonal system per
    step for all sheets and all three coordinates.
    """
    x = x.copy()
    ns, ni, nk, _ = x.shape
    lim = params.control_limit
    t = np.linspace(0.0, 1.0, nk)[None, None, :]
    s = np.linspace(0.0, 1.0, ni)[None, :, None]
    p_fn = (1 - t) * np.stack([_control(e, lim) for e in x[:, :, 0]])[:, :, None]
    p_fn = p_fn + t * np.stack([_control(e, lim) for e in x[:, :, -1]])[:, :, None]
    q_fn = (1 - s) * np.stack([_control(e, lim) for e in x[:, 0]])[:, None, :]
    q_fn = q_fn + s * np.stack([_control(e, lim) for e in x[:, -1]])[:, None, :]
    p_in, q_in = p_fn[:, 1:-1, 1:-1], q_fn[:, 1:-1, 1:-1]

    mi, mk = ni - 2, nk - 2
    shape = (ns, mi, mk)
    idx = np.arange(ns * mi * mk).reshape(shape)
    sheet = np.arange(ns)[:, None, None]
    for _ in range(params.max_iterations):
        x_i = 0.5 * (x[:, 2:, 1:-1] - x[:, :-2, 1:-1])
        x_k = 0.5 * (x[:, 1:-1, 2:] - x[:, 1:-1, :-2])
        alpha = np.einsum("...j,...j->...", x_k, x_k)
        beta = np.einsum("...j,...j->...", x_i, x_k)
        gamma = np.einsum("...j,...j->...", x_i, x_i)

        # alpha (E + W - 2 c) + alpha P (E - W) / 2 + gamma (N + S - 2 c) + gamma Q (N - S) / 2
        #   - beta (NE - NW - SE + SW) / 2 = 0
        half_b = 0.5 * beta
        stencil = (
            (alpha * (1.0 + 0.5 * p_in), 1, 0),
            (alpha * (1.0 - 0.5 * p_in), -1, 0),
            (gamma * (1.0 + 0.5 * q_in), 0, 1),
            (gamma * (1.0 - 0.5 * q_in), 0, -1),
            (-half_b, 1, 1),
            (half_b, 1, -1),
            (half_b, -1, 1),
            (-half_b, -1, -1),
        )
        rhs = np.zeros((*shape, 3))
        rows_, cols_, vals_ = [idx.ravel()], [idx.ravel()], [(-2.0 * (alpha + gamma)).ravel()]
        ii = np.arange(mi)[None, :, None]
        kk = np.arange(mk)[None, None, :]
        for coef, di, dk in stencil:
            ti, tk = ii + di, kk + dk
            inside = np.broadcast_to((ti >= 0) & (ti < mi) & (tk >= 0) & (tk < mk), shape)
            dst = idx[sheet, np.clip(ti, 0, mi - 1), np.clip(tk, 0, mk - 1)]
            rows_.append(idx[inside])
            cols_.append(dst[inside])
            vals_.append(coef[inside])
            edge = ~inside
            if edge.any():
                bs = np.broadcast_to(sheet, shape)[edge]
                bi = np.broadcast_to(ti + 1, shape)[edge]
                bk = np.broadcast_to(tk + 1, shape)[edge]
                rhs[edge] -= coef[edge][:, None] * x[bs, bi, bk]
        n = ns * mi * mk
        matrix = sp.csc_matrix(
            (np.concatenate(vals_), (np.concatenate(rows_), np.concatenate(cols_))), shape=(n, n)
        )
        new = spla.splu(matrix).solve(rhs.reshape(-1, 3)).reshape((*shape, 3))
        change = float(np.max(np.abs(new - x[:, 1:-1, 1:-1])))
        x[:, 1:-1, 1:-1] = new
        if change <= params.tolerance * scale:
            break
    return x
