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

import itertools
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
    # corners (ADR 0021)
    corner_deg: float = 30.0  # a step whose tangent turns more is checked for a corner
    corner_scale: float = 1e-6  # ... down to this step: a failure there is a corner
    ring_iterations: int = 200  # alternating projections onto the curve and a sphere


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
    return _loop(body_a, body_b, seed, n, spacing, curvature, entity_a, entity_b, params)[0]


def _loop(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    n: int,
    spacing: Spacing1D | None,
    curvature: float,
    entity_a: EntityRef,
    entity_b: EntityRef,
    params: TraceParams,
) -> tuple[np.ndarray, list[int]]:
    """``intersection_loop``'s nodes and the indices of its corner nodes (ADR 0024)."""
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
    polyline, corners = tracer.trace(start)
    if corners:  # a node on every corner (ADR 0021)
        nodes, fixed = _resample(
            tracer, polyline, [0, *corners], n, spacing, curvature, closed=True, with_fixed=True
        )
        return nodes, fixed[1:]

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
    return nodes, []


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
    ahead, corners_ahead = tracer.trace_to_plane(start, 1.0, point, oriented)
    behind, corners_behind = tracer.trace_to_plane(start, -1.0, point, oriented)
    polyline = np.vstack([behind[::-1], ahead[1:]])  # along n_b x n_a
    if corners_ahead or corners_behind:  # a node on every corner (ADR 0021)
        last = len(behind) - 1
        corners = sorted([last - i for i in corners_behind] + [last + i for i in corners_ahead])
        fixed = [0, *corners, len(polyline) - 1]
        return _resample(
            tracer, polyline, fixed, n, spacing, curvature, closed=False, plane=(point, normal)
        )

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
        self._sign = 1.0
        self._plane: tuple[np.ndarray, np.ndarray] | None = None

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

    def refine(
        self, polyline: np.ndarray, *, closed: bool, fixed: list[int] | None = None
    ) -> np.ndarray:
        """Bisect the segments next to vertices that turn more than
        ``refine_turn_deg``, new points polished onto both bodies, until none
        does or the segments reach ``refine_min`` (a corner). With *fixed*
        (vertex indices), returns ``(polyline, fixed)`` with them renumbered."""
        index = None if fixed is None else np.asarray(fixed, dtype=np.int64)
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
            if index is not None:
                index = index + np.searchsorted(where, index, side="right")
            polyline = np.insert(pts, where, mids, axis=0)
            if closed:
                polyline = polyline[:-1]
        if index is None:
            return polyline
        return polyline, [int(i) for i in index]  # type: ignore[return-value]

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

    def trace(self, start: np.ndarray) -> tuple[np.ndarray, list[int]]:
        """Predictor-corrector along n_b x n_a until the trace closes on *start*.

        Returns the polyline and the indices of its corners (ADR 0021).
        """
        self._sign = 1.0
        h = self.params.step * self.diag
        points, corners = [start], []
        travelled = 0.0
        t_next = None
        for _ in range(self.params.max_steps):
            x = points[-1]
            t = self.direction(x) if t_next is None else t_next
            t_next = None
            # a traced point only steers the next step: the output nodes are polished later
            y = self.polish((x + h * t)[None], self.params.trace_tolerance)[0]
            if not self._good(x, t, y, h):
                corner = self._fine(x, t, h / 16.0, 1.5 * h)
                if corner is not None:
                    c, t_next = corner
                    travelled += float(np.linalg.norm(c - x))
                    _add_corner(points, corners, c, self.params.on_curve * self.diag)
                    continue
                _check_progress(x, y, h)  # else a tight smooth bend: the plain step
            travelled += float(np.linalg.norm(y - x))
            if travelled > 4.0 * h and np.linalg.norm(y - start) < 1.5 * h:
                return np.array(points), corners
            points.append(y)
        raise ParameterError(
            f"the intersection trace did not close within {self.params.max_steps} steps."
        )

    def trace_to_plane(
        self, start: np.ndarray, sign: float, point: np.ndarray, normal: np.ndarray
    ) -> tuple[np.ndarray, list[int]]:
        """Trace along ``sign * n_b x n_a`` until the plane (*normal* to the kept side).

        A step that reaches the plane within the on-curve tolerance, or
        crosses it, ends the trace: on a half model the bodies stop at the
        plane, so the trace lands on it and would stall there. Returns the
        polyline and the indices of its corners (ADR 0021).
        """
        self._sign = sign
        self._plane = (point, normal)
        h = self.params.step * self.diag
        reach = self.params.on_curve * self.diag
        points, corners = [start], []
        travelled = 0.0
        t_next = None
        try:
            for _ in range(self.params.max_steps):
                x = points[-1]
                t = sign * self.direction(x) if t_next is None else t_next
                t_next = None
                y = self.polish((x + h * t)[None], self.params.trace_tolerance)[0]
                if float((y - point) @ normal) <= reach:
                    points.append(self.end_on_plane(y, point, normal))
                    return np.array(points), corners
                if not self._good(x, t, y, h):
                    corner = self._fine(x, t, h / 16.0, 1.5 * h)
                    if corner is not None:
                        c, t_next = corner
                        travelled += float(np.linalg.norm(c - x))
                        _add_corner(points, corners, c, self.params.on_curve * self.diag)
                        continue
                    stalled = np.linalg.norm(y - x) <= 0.5 * h
                    if stalled and float((y - point) @ normal) <= h:  # on a half body's cut
                        points.append(self.end_on_plane(y, point, normal))
                        return np.array(points), corners
                    _check_progress(x, y, h)
                travelled += float(np.linalg.norm(y - x))
                if travelled > 4.0 * h and np.linalg.norm(y - start) < 1.5 * h:
                    raise ParameterError(
                        "the intersection curve closes without meeting the plane "
                        "(use intersection_loop)."
                    )
                points.append(y)
        finally:
            self._plane = None
        raise ParameterError(
            f"the intersection trace did not reach the plane within {self.params.max_steps} steps."
        )

    # -- corners (ADR 0021) ----------------------------------------------------

    def _tangent(self, x: np.ndarray) -> np.ndarray:
        return self._sign * self.direction(x)

    def _good(self, x: np.ndarray, t: np.ndarray, y: np.ndarray, h: float) -> bool:
        """A plain step: it advanced, along *t*, and the tangent turned less
        than ``corner_deg`` (else a corner, or a bend tighter than the step)."""
        step = y - x
        if np.linalg.norm(step) <= 0.5 * h or _angle_deg(step, t) >= self.params.corner_deg:
            return False
        try:
            return _angle_deg(t, self._tangent(y)) < self.params.corner_deg
        except ParameterError:
            return False

    def _fine(
        self, x: np.ndarray, t: np.ndarray, h: float, length: float
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Walk about *length* from *x* with step *h*, refining recursively where
        a step fails; the first corner found (point, outgoing tangent), or None."""
        travelled = 0.0
        while travelled < length:
            y = self.polish((x + h * t)[None], self.params.trace_tolerance)[0]
            if self._plane is not None and float((y - self._plane[0]) @ self._plane[1]) <= h:
                return None  # the plane comes first (a half body may end there)
            if not self._good(x, t, y, h):
                if h <= self.params.corner_scale * self.diag:
                    return self._turn(x, t, h)
                corner = self._fine(x, t, h / 16.0, 1.5 * h)
                if corner is not None:
                    return corner
                _check_progress(x, y, h)
            travelled += float(np.linalg.norm(y - x))
            x, t = y, self._tangent(y)
        return None

    def _turn(self, x: np.ndarray, t1: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
        """The corner just ahead of *x* (incoming tangent *t1*), located, and the
        outgoing tangent."""
        # hints: back along t1, and fanned around it in A's tangent plane, where the
        # outgoing branch of a thin wedge lies (the 26 lattice directions miss it)
        axis = self.a.project(x[None], self.ea).normals[0]
        hints = [t1, -t1]
        for deg in (5.0, 10.0, 15.0, 20.0, 30.0, 45.0):
            for sign in (1.0, -1.0):
                hints.append(-_rotate(t1, axis, np.radians(sign * deg)))
        ring = self._ring(x, 2.0 * h, tuple(hints))
        # the way out turns most from the incoming tangent; a thin wedge puts it close to
        # the incoming branch, so distance alone cannot tell them apart
        turned = [(q, _angle_deg(self._tangent(q), t1)) for q in ring]
        out = [q for q, angle in turned if angle >= self.params.corner_deg]
        if out:
            p = max(out, key=lambda q: _angle_deg(self._tangent(q), t1))
        else:
            back = x - 2.0 * h * t1
            ahead = [q for q in ring if np.linalg.norm(q - back) > 0.5 * h]
            if not ahead:
                raise ParameterError(
                    f"the intersection curve has a corner near {x.tolist()} with no outgoing "
                    "branch."
                )
            p = max(ahead, key=lambda q: float(np.linalg.norm(q - back)))
        t2 = self._tangent(p)
        c = self._locate(x, t1, p, t2)
        if c is None:
            raise ParameterError(
                f"the intersection trace stalls near {x.tolist()} without a corner there."
            )
        return c, t2

    def _locate(
        self, x: np.ndarray, t1: np.ndarray, p: np.ndarray, t2: np.ndarray
    ) -> np.ndarray | None:
        """The corner between *x* (tangent *t1*) and *p* (tangent *t2*): the
        crossing of the two branches' tangent lines, from points ever closer to
        it; None when the branches turn less than ``corner_deg`` there."""
        c = _crossing(x, t1, p, t2)
        d = max(float(np.linalg.norm(c - x)), float(np.linalg.norm(c - p)))
        while d > 1e-9 * self.diag:
            d /= 8.0
            ring = self._ring(c, d, (t1, -t1, t2, -t2))
            before = [q for q in ring if (q - c) @ t1 < 0.0]
            after = [q for q in ring if (q - c) @ t2 > 0.0]
            if not before or not after:
                break
            p1 = min(before, key=lambda q: _angle_deg(self._tangent(q), t1))
            p2 = min(after, key=lambda q: _angle_deg(self._tangent(q), t2))
            if np.linalg.norm(p1 - p2) <= 0.05 * d:
                break  # below the on-curve tolerance the two branches merge: c stands
            t1, t2 = self._tangent(p1), self._tangent(p2)
            if _angle_deg(t1, t2) < self.params.corner_deg:
                return None
            c = _crossing(p1, t1, p2, t2)
        return self.polish(c[None])[0]

    def _ring(self, c: np.ndarray, r: float, hints: tuple[np.ndarray, ...]) -> list[np.ndarray]:
        """The curve's points at distance *r* from *c*: alternating projections
        onto A, B and the sphere, from 26 directions and the *hints*."""
        dirs = np.vstack([_DIRECTIONS, *(np.asarray(h)[None] / np.linalg.norm(h) for h in hints)])
        x = c + r * dirs
        tol = self.params.tolerance * self.diag
        done = np.zeros(len(x), dtype=bool)
        last = np.full(len(x), np.inf)
        for _ in range(self.params.ring_iterations):
            active = np.flatnonzero(~done)
            if not len(active):
                break
            y = self.b.project(self.a.project(x[active], self.ea).points, self.eb).points
            v = y - c
            norm = np.linalg.norm(v, axis=1, keepdims=True)
            y = np.where(norm > 0.0, c + r * v / np.where(norm > 0.0, norm, 1.0), y)
            moved = np.linalg.norm(y - x[active], axis=1)
            x[active] = y
            last[active] = moved
            done[active[moved <= tol]] = True
        keep = done & self.on_curve(x) & (np.abs(np.linalg.norm(x - c, axis=1) - r) <= 1e-6 * r)
        # a projection's own noise (sliver STL triangles) can sit above the tolerances:
        # a point settled to a small fraction of r is good enough to steer the corner search
        loose = 1e-2 * r
        if not keep.all():
            da = self.a.project(x, self.ea).distances
            db = self.b.project(x, self.eb).distances
            settled = (last <= loose) & (da <= loose) & (db <= loose)
            keep |= settled & (np.abs(np.linalg.norm(x - c, axis=1) - r) <= loose)
        out: list[np.ndarray] = []
        for q in x[keep]:
            if all(np.linalg.norm(q - o) > 0.05 * r for o in out):
                out.append(q)
        return out

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


def _resample(
    tracer: _Tracer,
    polyline: np.ndarray,
    fixed: list[int],
    n: int,
    spacing: Spacing1D | None,
    curvature: float,
    *,
    closed: bool,
    plane: tuple[np.ndarray, np.ndarray] | None = None,
    with_fixed: bool = False,
) -> np.ndarray | tuple[np.ndarray, list[int]]:
    """Nodes of a curve with fixed vertices (the start or ends, and corners):
    each segment between them gets cells in proportion to its length (its
    monitor integral with *curvature*) and the law of ADR 0020 by itself.
    *with_fixed*: also the node index of each fixed vertex."""
    params = tracer.params
    polyline, fixed = tracer.refine(polyline, closed=closed, fixed=fixed)  # type: ignore[misc]
    pts = np.vstack([polyline, polyline[:1]]) if closed else polyline
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    ends = [*fixed, len(pts) - 1] if closed else list(fixed)
    total = n if closed else n - 1
    if total < len(ends) - 1:
        raise ParameterError(
            f"n = {n} is too small for the {len(ends) - 1} segments between the curve's "
            "corners: give each at least one cell."
        )
    kappa = None
    if curvature > 0.0:
        if closed:
            kappa = _windowed_curvature(polyline, s, n, params)
            kappa = np.r_[kappa, kappa[:1]]
        else:
            point, normal = plane  # type: ignore[misc]
            inner = polyline[-2:0:-1]
            loop = np.vstack([polyline, inner - 2.0 * ((inner - point) @ normal)[:, None] * normal])
            closed_loop = np.vstack([loop, loop[:1]])
            sl = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(closed_loop, axis=0), axis=1))]
            kappa = _windowed_curvature(loop, sl, 2 * (n - 1), params)[: len(polyline)]
        m = 1.0 + curvature * kappa * s[-1] / (2.0 * np.pi)
        measure = np.r_[0.0, np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(s))]
    else:
        measure = s
    if kappa is None:
        cells = _allocate(np.diff(measure[ends]), total)
        parts = [np.full(int(c), (s[b] - s[a]) / int(c)) for (a, b), c in
                 zip(itertools.pairwise(ends), cells, strict=True)]  # fmt: skip
    else:
        # a graded spacing function over the whole curve, corners included (ADR 0022)
        m = 1.0 + curvature * kappa * s[-1] / (2.0 * np.pi)
        g = _graded_spacing(s, 1.0 / m, total, float(np.log(params.max_cell_ratio)), closed)
        # integrate 1 / g finer than the cells: the polyline's vertices can be a trace
        # step apart where the cells are much smaller
        fine = np.unique(np.r_[s, np.arange(0.0, s[-1], 0.25 * g.min())])
        gf = np.interp(fine, s, g)
        sigma_f = np.r_[0.0, np.cumsum(0.5 * (1.0 / gf[1:] + 1.0 / gf[:-1]) * np.diff(fine))]
        sigma = np.interp(s, fine, sigma_f)
        cells = _allocate(np.diff(sigma[ends]), total)
        parts = []
        for (a, b), count in zip(itertools.pairwise(ends), cells, strict=True):
            sel = (fine >= s[a]) & (fine <= s[b])
            local = (sigma_f[sel] - sigma[a]) / (sigma[b] - sigma[a])
            x = np.interp(np.linspace(0.0, 1.0, int(count) + 1), local, fine[sel] - s[a])
            parts.append(np.diff(x))
    target, node_of_fixed = [], []
    for (a, b), h in zip(itertools.pairwise(ends), parts, strict=True):
        node_of_fixed.append(sum(len(t) for t in target))
        x = np.r_[0.0, np.cumsum(h)]
        x[-1] = s[b] - s[a]
        if spacing is not None:
            xi = np.linspace(0.0, 1.0, len(x))
            x = np.interp(spacing.build(len(x), 0.0, 1.0), xi, x)
        target.append(s[a] + x[:-1])
    if not closed:
        node_of_fixed.append(sum(len(t) for t in target))
        target.append(s[-1:])
    target = np.concatenate(target)
    nodes = np.stack([np.interp(target, s, pts[:, k]) for k in range(3)], axis=-1)
    free = np.ones(len(nodes), dtype=bool)
    free[node_of_fixed] = False
    nodes[free] = tracer.polish(nodes[free])
    nodes[node_of_fixed] = pts[ends[: len(node_of_fixed)]]
    if with_fixed:
        return nodes, node_of_fixed
    return nodes


def _graded_spacing(
    s: np.ndarray, shape: np.ndarray, cells: int, slope: float, closed: bool
) -> np.ndarray:
    """Spacing function at the vertices *s*: ``C * shape`` capped by cones of
    *slope* (so adjacent cells differ by at most ``1 + slope``), with ``C``
    such that the curve holds *cells* cells (ADR 0022)."""

    def capped(c: float) -> np.ndarray:
        g = c * shape
        if closed:
            length = s[-1]
            ext_s = np.r_[s[:-1] - length, s[:-1], s[:-1] + length]
            ext = np.r_[g[:-1], g[:-1], g[:-1]]
        else:
            ext_s, ext = s, g.copy()
        for i in range(1, len(ext)):
            ext[i] = min(ext[i], ext[i - 1] + slope * (ext_s[i] - ext_s[i - 1]))
        for i in range(len(ext) - 2, -1, -1):
            ext[i] = min(ext[i], ext[i + 1] + slope * (ext_s[i + 1] - ext_s[i]))
        if closed:
            mid = ext[len(s) - 1 : 2 * (len(s) - 1)]
            return np.r_[mid, mid[:1]]
        return ext

    def count(c: float) -> float:
        g = capped(c)
        return float(np.sum(0.5 * (1.0 / g[1:] + 1.0 / g[:-1]) * np.diff(s)))

    lo, hi = 1e-12, 1.0
    while count(hi) > cells:
        hi *= 2.0
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if count(mid) > cells else (lo, mid)
        if hi / lo < 1.0 + 1e-12:
            break
    return capped(hi)


def _allocate(measure: np.ndarray, total: int) -> np.ndarray:
    """Cells per segment in proportion to *measure*, each at least one
    (largest remainders)."""
    share = measure / measure.sum() * total
    cells = np.maximum(1, np.floor(share)).astype(np.int64)
    while cells.sum() > total:
        over = np.where(cells > 1, cells - share, -np.inf)
        cells[int(np.argmax(over))] -= 1
    while cells.sum() < total:
        cells[int(np.argmax(share - cells))] += 1
    return cells


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


# -- corners (ADR 0021) ----------------------------------------------------------

_DIRECTIONS = np.array(
    [
        [i, j, k]
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)
        if (i, j, k) != (0, 0, 0)
    ],
    dtype=np.float64,
)
_DIRECTIONS /= np.linalg.norm(_DIRECTIONS, axis=1, keepdims=True)


def _angle_deg(u: np.ndarray, v: np.ndarray) -> float:
    c = float(u @ v) / (float(np.linalg.norm(u)) * float(np.linalg.norm(v)))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def _rotate(v: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """*v* rotated by *angle* about the unit *axis* (Rodrigues)."""
    k = axis / np.linalg.norm(axis)
    return v * np.cos(angle) + np.cross(k, v) * np.sin(angle) + k * (k @ v) * (1.0 - np.cos(angle))


def _crossing(p1: np.ndarray, t1: np.ndarray, p2: np.ndarray, t2: np.ndarray) -> np.ndarray:
    """Midpoint of the closest points of two lines (of parallel ones: of *p1*, *p2*)."""
    w = p1 - p2
    a, b, c = t1 @ t1, t1 @ t2, t2 @ t2
    d, e = t1 @ w, t2 @ w
    den = a * c - b * b
    if den <= 1e-12 * a * c:
        return 0.5 * (p1 + p2)
    return 0.5 * (p1 + (b * e - c * d) / den * t1 + p2 + (a * e - b * d) / den * t2)


def _check_progress(x: np.ndarray, y: np.ndarray, h: float) -> None:
    if np.linalg.norm(y - x) <= 0.5 * h:
        raise ParameterError(f"the intersection trace stalls near {x.tolist()} without a corner.")


def _add_corner(points: list[np.ndarray], corners: list[int], c: np.ndarray, tol: float) -> None:
    """Append corner *c*; a traced point already on it (a step that landed on
    the corner) becomes the corner instead."""
    if len(points) > 1 and np.linalg.norm(points[-1] - c) <= tol:
        points[-1] = c
    else:
        points.append(c)
    corners.append(len(points) - 1)
