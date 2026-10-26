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
are put on the plane and marked ``FaceBC.SYMMETRY``.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from shore.exceptions import ParameterError
from shore.geometry.base import Geometry
from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC
from shore.surface.intersect import intersection_curve, intersection_loop
from shore.surface.march import _Mirror, march_surface
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

    mirror = None if plane is None else _Mirror.of(plane)
    if mirror is None:
        curve = intersection_loop(body_a, body_b, seed, int(n_loop))
        march = {"periodic": True}
    else:
        curve = intersection_curve(body_a, body_b, seed, int(n_loop), plane=plane)
        if _end_angle_deg(curve, mirror.normal) > params.oblique_deg:
            warnings.warn(
                "the junction curve meets the symmetry plane obliquely: the bodies are not "
                "mirror-symmetric there, and the half-collar's quality is not guaranteed.",
                UserWarning,
                stacklevel=2,
            )
        march = {"plane": plane}
    on_a = march_surface(body_a, curve, steps_a, **march)
    on_b = march_surface(body_b, curve[::-1], steps_b, **march)[::-1]
    wall = np.concatenate([on_b[:, ::-1], on_a[:, 1:]], axis=1)
    wall = np.ascontiguousarray(np.transpose(wall, (1, 0, 2))[:, ::-1])
    junction = len(steps_b)

    if mirror is None:
        nodes = extrude_collar(wall, steps_k, junction, params=params)
        return HexBlock(label, nodes, periodic_i=False, periodic_j=True)
    closed = np.concatenate([wall, mirror.reflect(wall[:, -2:0:-1])], axis=1)
    nodes = extrude_collar(closed, steps_k, junction, params=params)[:, : int(n_loop)].copy()
    nodes[:, [0, -1]] = mirror.onto(nodes[:, [0, -1]])
    block = HexBlock(label, nodes, periodic_i=False, periodic_j=False)
    block._set_face_bc("j_lo", FaceBC.SYMMETRY)
    block._set_face_bc("j_hi", FaceBC.SYMMETRY)
    return block


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
    h = float(
        min(
            np.linalg.norm(w[junction + 1] - w[junction], axis=-1).min(),
            np.linalg.norm(w[junction - 1] - w[junction], axis=-1).min(),
        )
    )
    return _extrude(w, steps_k, h, junction, params)


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
    wall: np.ndarray, steps_k: np.ndarray, h: float, junction: int, params: CollarParams
) -> np.ndarray:
    ni, nj, _ = wall.shape
    nk = len(steps_k)
    cum = np.r_[0.0, np.cumsum(steps_k)]
    k0 = max(1, int(np.searchsorted(cum, params.hyperbolic_fraction * h)) - 1)
    k0 = min(k0, nk)

    nodes = np.empty((ni, nj, nk + 1, 3))
    nodes[:, :, 0] = wall
    for k in range(k0):
        nodes[:, :, k + 1] = hyperbolic_layer(nodes[:, :, k], float(steps_k[k]), (False, True))
    if k0 == nk:
        return nodes

    inner = nodes[:, :, k0]
    normals = _normals(inner)
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


def _normals(layer: np.ndarray) -> np.ndarray:
    """Unit normals of a j-periodic layer (one-sided at the open i-ends)."""
    ti = np.gradient(layer, axis=0)
    tj = 0.5 * (np.roll(layer, -1, axis=1) - np.roll(layer, 1, axis=1))
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
