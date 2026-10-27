"""Intersection loop of two bodies, the start curve of a collar grid (ADR 0009).

``intersection_loop`` finds the closed curve where two bodies meet, near a
seed point, with ``Geometry.project`` only (so STL, CAD and analytic bodies
work alike):

1. alternating projections onto A and B take the seed onto the curve;
2. the curve is traced with predictor steps along ``t = n_b x n_a`` (the
   bodies' outward normals), each corrected back onto both bodies by
   alternating projections, until the trace returns to its start;
3. the traced polyline is resampled to *n* nodes equally spaced in arc
   length, and every node is polished back onto both bodies.

The ``n_b x n_a`` order makes ``march_surface(A, loop)`` march on A away
from B, and ``march_surface(B, loop[::-1])`` on B away from A: the two halves
of a collar (ADR 0008). Alternating projections converge fast where the
bodies cross at a clear angle; grazing intersections are out of scope.

``intersection_curve`` is the open variant for half models (ADR 0011): the
curve on the seed's side of a plane, traced both ways from the seed until
it reaches the plane (within the on-curve tolerance: half-model bodies stop
at the plane, so the trace lands on it without crossing), each end put
exactly on the plane by Newton steps along the curve tangent.

Both take a node distribution along the curve (ADR 0020): a ``Spacing1D``
law by arc length and/or ``curvature`` clustering, which equidistributes
``1 + alpha * kappa * L / (2 pi)`` on the traced polyline (refined first
where it turns, so that the curvature and the arc lengths are resolved)
before the law is applied, with adjacent cells limited to
``TraceParams.max_cell_ratio``. The defaults keep the plain arc-length
resampling above, bit for bit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from shore.exceptions import ParameterError
from shore.geometry.base import BODY, EntityRef, Geometry

if TYPE_CHECKING:
    from shore.mesh.spacing import Spacing1D


@dataclass(frozen=True)
class TraceParams:
    """Tracer settings (not per-case tuning knobs); lengths relative to the model diagonal."""

    step: float = 4e-3  # predictor step along the curve
    tolerance: float = 1e-13  # alternating-projection convergence of the output nodes
    trace_tolerance: float = 1e-10  # looser convergence of the traced points
    on_curve: float = 1e-9  # a point is on the curve when this close to both bodies
    polish_iterations: int = 200
    max_steps: int = 200_000
    # node distribution (ADR 0020)
    refine_turn_deg: float = 2.0  # polyline refined until it turns less per segment
    refine_min: float = 1e-5  # ... or its segments are this short (a corner)
    curvature_window: float = 0.5  # curvature averaged over this fraction of a uniform cell
    max_cell_ratio: float = 1.2  # adjacent cells of the curvature distribution


DEFAULT = TraceParams()


def intersection_loop(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    n: int,
    *,
    spacing: Spacing1D | None = None,
    curvature: float = 0.0,
    entity_a: EntityRef = BODY,
    entity_b: EntityRef = BODY,
    params: TraceParams = DEFAULT,
) -> np.ndarray:
    """The closed intersection loop of *body_a* and *body_b* nearest *seed*.

    Parameters
    ----------
    body_a, body_b :
        The two bodies; only ``project`` (with normals) is used.
    seed : np.ndarray, shape (3,)
        A point near the junction; selects the loop when there are several.
    n :
        Number of nodes, >= 3.
    spacing :
        Node distribution by arc length from node 0 round to node 0 again:
        ``uniform`` or ``tanh2`` (default: uniform).
    curvature :
        Clustering weight *alpha* >= 0 where the loop bends (ADR 0020).
    entity_a, entity_b :
        The entities intersected (default: the whole bodies).
    params :
        Tracer settings.

    Returns
    -------
    np.ndarray, shape (n, 3)
        Nodes distributed by *spacing* and *curvature* (by default equally
        spaced in arc length), on both bodies, ordered along
        ``n_b x n_a``; ``[0]`` is where the seed lands; the last node is not
        repeated.

    Raises
    ------
    ParameterError
        n < 3; *seed* not a 3-vector; a law other than ``uniform`` or
        ``tanh2``, or a bad *curvature*; no intersection near *seed*; a
        trace that does not close.
    """
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 3:
        raise ParameterError(f"n must be an integer >= 3, got {n!r}.")
    default = _check_distribution(spacing, curvature, ("uniform", "tanh2"))
    x = np.asarray(seed, dtype=np.float64)
    if x.shape != (3,) or not np.all(np.isfinite(x)):
        raise ParameterError(f"seed must be a finite 3-vector, got shape {x.shape}.")

    bounds = np.asarray(body_a.bounds, dtype=np.float64)
    diag = float(np.linalg.norm(bounds[1] - bounds[0]))
    tracer = _Tracer(body_a, body_b, entity_a, entity_b, diag, params)

    start = tracer.polish(x[None])[0]
    if not tracer.on_curve(start[None]).all():
        raise ParameterError(f"no intersection of the two bodies near the seed {x.tolist()}.")
    polyline = tracer.trace(start)

    # resample by arc length (closing segment included) and polish onto both bodies
    if not default:
        polyline = tracer.refine(polyline, closed=True)
    closed = np.vstack([polyline, polyline[:1]])
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
    if default:
        target = np.linspace(0.0, s[-1], n + 1)[:-1]
    else:
        k = _windowed_curvature(polyline, s, n, params) if curvature > 0.0 else None
        target = _distribution(s, n, spacing, curvature, k, params, closed=True)[:-1]
    nodes = np.stack([np.interp(target, s, closed[:, k]) for k in range(3)], axis=-1)
    nodes = tracer.polish(nodes)
    nodes[0] = start
    return nodes


def intersection_curve(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    n: int,
    *,
    plane: tuple[np.ndarray, np.ndarray],
    spacing: Spacing1D | None = None,
    curvature: float = 0.0,
    entity_a: EntityRef = BODY,
    entity_b: EntityRef = BODY,
    params: TraceParams = DEFAULT,
) -> np.ndarray:
    """The open intersection curve of *body_a* and *body_b* ending on *plane*.

    Parameters
    ----------
    body_a, body_b :
        The two bodies; only ``project`` (with normals) is used. They may end
        at the plane (a half model).
    seed : np.ndarray, shape (3,)
        A point near the junction, off the plane; selects the curve and the
        side of the plane that is kept.
    n :
        Number of nodes, >= 3.
    plane :
        ``(point, normal)``, two 3-vectors; the normal's sign does not matter.
    spacing :
        Node distribution by arc length from the first end to the last:
        ``uniform``, ``tanh`` or ``tanh2`` (default: uniform).
    curvature :
        Clustering weight *alpha* >= 0 where the curve, or its mirror image
        through the plane, bends (ADR 0020).
    entity_a, entity_b :
        The entities intersected (default: the whole bodies).
    params :
        Tracer settings.

    Returns
    -------
    np.ndarray, shape (n, 3)
        Nodes distributed by *spacing* and *curvature* (by default equally
        spaced in arc length), on both bodies and on the seed's side, ordered
        along ``n_b x n_a``; ``[0]`` and ``[-1]`` on the plane.

    Raises
    ------
    ParameterError
        n < 3; *seed* not a 3-vector; *plane* not a point and a non-zero
        normal; a law other than ``uniform``, ``tanh`` or ``tanh2``, or a bad
        *curvature*; *seed* on the plane; no intersection near *seed*; a curve
        that closes, or meets the plane tangentially, or does not reach it.
    """
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 3:
        raise ParameterError(f"n must be an integer >= 3, got {n!r}.")
    default = _check_distribution(spacing, curvature, ("uniform", "tanh", "tanh2"))
    x = np.asarray(seed, dtype=np.float64)
    if x.shape != (3,) or not np.all(np.isfinite(x)):
        raise ParameterError(f"seed must be a finite 3-vector, got shape {x.shape}.")
    try:
        point, normal = (np.asarray(v, dtype=np.float64) for v in plane)
    except (TypeError, ValueError):
        raise ParameterError("plane must be (point, normal), two 3-vectors.") from None
    if point.shape != (3,) or normal.shape != (3,):
        raise ParameterError("plane must be (point, normal), two 3-vectors.")
    length = float(np.linalg.norm(normal))
    if not (np.all(np.isfinite(point)) and np.isfinite(length) and length > 0.0):
        raise ParameterError("plane must have a finite point and a non-zero normal.")
    normal = normal / length

    bounds = np.asarray(body_a.bounds, dtype=np.float64)
    diag = float(np.linalg.norm(bounds[1] - bounds[0]))
    tracer = _Tracer(body_a, body_b, entity_a, entity_b, diag, params)
    side = float((x - point) @ normal)
    if abs(side) <= params.on_curve * diag:
        raise ParameterError("the seed lies on the plane: it must select a side.")
    oriented = normal * np.sign(side)  # points to the kept side

    start = tracer.polish(x[None])[0]
    if not tracer.on_curve(start[None]).all():
        raise ParameterError(f"no intersection of the two bodies near the seed {x.tolist()}.")
    if float((start - point) @ oriented) <= params.on_curve * diag:
        raise ParameterError("the seed's junction point is not on the seed's side of the plane.")
    ahead = tracer.trace_to_plane(start, 1.0, point, oriented)
    behind = tracer.trace_to_plane(start, -1.0, point, oriented)
    polyline = np.vstack([behind[::-1], ahead[1:]])  # along n_b x n_a

    if not default:
        polyline = tracer.refine(polyline, closed=False)
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(polyline, axis=0), axis=1))]
    if default:
        target = np.linspace(0.0, s[-1], n)[1:-1]
    else:
        k = None
        if curvature > 0.0:
            # measured on the mirrored closed loop, so the plane ends see their bend
            mirrored = (
                polyline[-2:0:-1] - 2.0 * ((polyline[-2:0:-1] - point) @ normal)[:, None] * normal
            )
            loop = np.vstack([polyline, mirrored])
            sl = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(np.vstack([loop, loop[:1]]),
                                                             axis=0), axis=1))]  # fmt: skip
            k = _windowed_curvature(loop, sl, 2 * (n - 1), params)[: len(polyline)]
        target = _distribution(s, n, spacing, curvature, k, params, closed=False)[1:-1]
    inner = np.stack([np.interp(target, s, polyline[:, k]) for k in range(3)], axis=-1)
    return np.vstack([polyline[:1], tracer.polish(inner), polyline[-1:]])


class _Tracer:
    def __init__(
        self,
        body_a: Geometry,
        body_b: Geometry,
        entity_a: EntityRef,
        entity_b: EntityRef,
        diag: float,
        params: TraceParams,
    ) -> None:
        self.a, self.b, self.ea, self.eb = body_a, body_b, entity_a, entity_b
        self.diag, self.params = diag, params

    def polish(self, x: np.ndarray, tolerance: float | None = None) -> np.ndarray:
        """Alternating projections onto A and B until the points stop moving."""
        tol = (self.params.tolerance if tolerance is None else tolerance) * self.diag
        for _ in range(self.params.polish_iterations):
            y = self.b.project(self.a.project(x, self.ea).points, self.eb).points
            moved = float(np.max(np.linalg.norm(y - x, axis=-1)))
            x = y
            if moved <= tol:
                break
        return x

    def refine(self, polyline: np.ndarray, *, closed: bool) -> np.ndarray:
        """Bisect the segments next to vertices that turn more than
        ``refine_turn_deg``, new points polished onto both bodies, until none
        does or the segments reach ``refine_min`` (a corner)."""
        max_turn = np.radians(self.params.refine_turn_deg)
        min_len = self.params.refine_min * self.diag
        for _ in range(30):
            pts = np.vstack([polyline, polyline[:1]]) if closed else polyline
            seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
            turn = _turning(polyline, closed)
            turn = np.r_[turn, turn[:1]] if closed else turn
            split = (np.maximum(turn[:-1], turn[1:]) > max_turn) & (seg > min_len)
            if not split.any():
                break
            mids = self.polish(
                0.5 * (pts[:-1][split] + pts[1:][split]), self.params.trace_tolerance
            )
            where = np.flatnonzero(split) + 1
            polyline = np.insert(pts, where, mids, axis=0)
            if closed:
                polyline = polyline[:-1]
        return polyline

    def on_curve(self, x: np.ndarray) -> np.ndarray:
        tol = self.params.on_curve * self.diag
        da = self.a.project(x, self.ea).distances
        db = self.b.project(x, self.eb).distances
        return (da <= tol) & (db <= tol)

    def direction(self, x: np.ndarray) -> np.ndarray:
        na = self.a.project(x[None], self.ea).normals[0]
        nb = self.b.project(x[None], self.eb).normals[0]
        t = np.cross(nb, na)
        norm = float(np.linalg.norm(t))
        if norm < 1e-8:
            raise ParameterError(
                "the bodies meet tangentially here (parallel normals): the curve has no "
                "tracing direction (grazing intersections are out of scope, ADR 0009)."
            )
        return t / norm

    def trace(self, start: np.ndarray) -> np.ndarray:
        """Predictor-corrector along n_b x n_a until the trace closes on *start*."""
        h = self.params.step * self.diag
        points = [start]
        travelled = 0.0
        for _ in range(self.params.max_steps):
            x = points[-1]
            # a traced point only steers the next step: the output nodes are polished later
            y = self.polish((x + h * self.direction(x))[None], self.params.trace_tolerance)[0]
            travelled += float(np.linalg.norm(y - x))
            if travelled > 4.0 * h and np.linalg.norm(y - start) < 1.5 * h:
                return np.array(points)
            points.append(y)
        raise ParameterError(
            f"the intersection trace did not close within {self.params.max_steps} steps."
        )

    def trace_to_plane(
        self, start: np.ndarray, sign: float, point: np.ndarray, normal: np.ndarray
    ) -> np.ndarray:
        """Trace along ``sign * n_b x n_a`` until the plane (*normal* to the kept side).

        A step that reaches the plane within the on-curve tolerance, or
        crosses it, ends the trace: on a half model the bodies stop at the
        plane, so the trace lands on it and would stall there.
        """
        h = self.params.step * self.diag
        reach = self.params.on_curve * self.diag
        points = [start]
        travelled = 0.0
        for _ in range(self.params.max_steps):
            x = points[-1]
            y = self.polish((x + sign * h * self.direction(x))[None], self.params.trace_tolerance)
            y = y[0]
            if float((y - point) @ normal) <= reach:
                points.append(self.end_on_plane(y, point, normal))
                return np.array(points)
            travelled += float(np.linalg.norm(y - x))
            if travelled > 4.0 * h and np.linalg.norm(y - start) < 1.5 * h:
                raise ParameterError(
                    "the intersection curve closes without meeting the plane "
                    "(use intersection_loop)."
                )
            points.append(y)
        raise ParameterError(
            f"the intersection trace did not reach the plane within {self.params.max_steps} steps."
        )

    def end_on_plane(self, x: np.ndarray, point: np.ndarray, normal: np.ndarray) -> np.ndarray:
        """Newton steps along the curve tangent onto the plane, polished onto both bodies."""
        tol = self.params.tolerance * self.diag
        x = self.polish(x[None])[0]  # a traced point is only converged to trace_tolerance
        for _ in range(self.params.polish_iterations):
            t = self.direction(x)
            slope = float(t @ normal)
            if abs(slope) < 1e-8:
                raise ParameterError(
                    "the intersection curve meets the plane tangentially (out of scope, ADR 0011)."
                )
            d = float((x - point) @ normal)
            if abs(d) <= tol:
                break
            x = self.polish((x - d / slope * t)[None])[0]
        return x


# -- node distribution along the curve (ADR 0020) --------------------------------


def _check_distribution(spacing: Spacing1D | None, curvature: float, laws: tuple[str, ...]) -> bool:
    """Validate *spacing* and *curvature*; True for the default resampling."""
    try:
        alpha = float(curvature)
    except (TypeError, ValueError):
        raise ParameterError(f"curvature must be a number >= 0, got {curvature!r}.") from None
    if not (np.isfinite(alpha) and alpha >= 0.0):
        raise ParameterError(f"curvature must be a finite number >= 0, got {curvature!r}.")
    if spacing is None:
        return alpha == 0.0
    law = getattr(spacing, "law", None)
    if law not in laws:
        raise ParameterError(f"spacing law along this curve must be one of {laws}, got {law!r}.")
    at = getattr(spacing, "cluster_at", None)
    if at is not None and not (np.isfinite(at) and 0.0 <= at <= 1.0):
        raise ParameterError(f"spacing cluster_at must lie in [0, 1], got {at!r}.")
    return law == "uniform" and alpha == 0.0


def _turning(polyline: np.ndarray, closed: bool) -> np.ndarray:
    """Turning angle (radians) at each vertex; 0 at the ends of an open polyline."""
    if closed:
        prev, here, nxt = np.roll(polyline, 1, axis=0), polyline, np.roll(polyline, -1, axis=0)
    else:
        prev, here, nxt = polyline[:-2], polyline[1:-1], polyline[2:]
    t1, t2 = here - prev, nxt - here
    c = np.einsum("ij,ij->i", t1, t2) / (np.linalg.norm(t1, axis=1) * np.linalg.norm(t2, axis=1))
    angle = np.arccos(np.clip(c, -1.0, 1.0))
    return angle if closed else np.r_[0.0, angle, 0.0]


def _windowed_curvature(loop: np.ndarray, s: np.ndarray, n: int, params: TraceParams) -> np.ndarray:
    """Curvature at each vertex of a closed polyline, averaged over a window
    of ``curvature_window`` uniform cells (``s[-1] / n``).

    Each vertex's turning angle is spread evenly over its two half-segments;
    the window integrates that density exactly, so the estimate is smooth
    whatever the window's size relative to the segments. *s* is the arc
    length at the vertices, closing vertex included.
    """
    length = s[-1]
    w = params.curvature_window * length / n
    turn = _turning(loop, closed=True)
    mid = 0.5 * (s[:-1] + s[1:])  # vertex i lies between mid[i - 1] and mid[i]
    # cumulative turning, linear between midpoints, over three periods
    knots = np.r_[mid - length, mid, mid + length]
    cum = np.cumsum(np.tile(turn, 3))
    sv = s[:-1]
    return (np.interp(sv + 0.5 * w, knots, cum) - np.interp(sv - 0.5 * w, knots, cum)) / w


def _limit_ratio(h: np.ndarray, ratio: float, closed: bool) -> np.ndarray:
    """Shrink cells until adjacent ones differ by at most *ratio*."""
    log_r = np.log(ratio)
    for _ in range(4):
        x = np.log(np.r_[h, h, h] if closed else h)
        i = np.arange(len(x)) * log_r
        x = np.minimum(x, np.minimum.accumulate(x - i) + i)  # h_k <= r h_(k-1)
        x = np.minimum(x, (np.minimum.accumulate((x + i)[::-1]) - i[::-1])[::-1])
        new = np.minimum(h, np.exp(x[len(h) : 2 * len(h)] if closed else x))
        if np.array_equal(new, h):
            break
        h = new
    return h


def _distribution(
    s: np.ndarray,
    n: int,
    spacing: Spacing1D | None,
    curvature: float,
    kappa: np.ndarray | None,
    params: TraceParams,
    *,
    closed: bool,
) -> np.ndarray:
    """Arc lengths of the nodes (closed: the closing node included)."""
    cells = n if closed else n - 1
    length = s[-1]
    xi = np.linspace(0.0, 1.0, cells + 1)
    if kappa is None:
        x = xi * length
    else:
        k = np.r_[kappa, kappa[:1]] if closed else kappa
        m = 1.0 + curvature * k * length / (2.0 * np.pi)
        sigma = np.r_[0.0, np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(s))]
        x = np.interp(xi, sigma / sigma[-1], s)
        h = _limit_ratio(np.diff(x), params.max_cell_ratio, closed)
        x = np.r_[0.0, np.cumsum(h * (length / h.sum()))]
    law = xi if spacing is None else spacing.build(cells + 1, 0.0, 1.0)
    return np.interp(law, xi, x)
