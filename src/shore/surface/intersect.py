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
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from shore.exceptions import ParameterError
from shore.geometry.base import BODY, EntityRef, Geometry


@dataclass(frozen=True)
class TraceParams:
    """Tracer settings (not per-case tuning knobs); lengths relative to the model diagonal."""

    step: float = 4e-3  # predictor step along the curve
    tolerance: float = 1e-13  # alternating-projection convergence of the output nodes
    trace_tolerance: float = 1e-10  # looser convergence of the traced points
    on_curve: float = 1e-9  # a point is on the curve when this close to both bodies
    polish_iterations: int = 200
    max_steps: int = 200_000


DEFAULT = TraceParams()


def intersection_loop(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    n: int,
    *,
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
    entity_a, entity_b :
        The entities intersected (default: the whole bodies).
    params :
        Tracer settings.

    Returns
    -------
    np.ndarray, shape (n, 3)
        Nodes equally spaced in arc length, on both bodies, ordered along
        ``n_b x n_a``; ``[0]`` is where the seed lands; the last node is not
        repeated.

    Raises
    ------
    ParameterError
        n < 3; *seed* not a 3-vector; no intersection near *seed*; a trace
        that does not close.
    """
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 3:
        raise ParameterError(f"n must be an integer >= 3, got {n!r}.")
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
    closed = np.vstack([polyline, polyline[:1]])
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
    target = np.linspace(0.0, s[-1], n + 1)[:-1]
    nodes = np.stack([np.interp(target, s, closed[:, k]) for k in range(3)], axis=-1)
    nodes = tracer.polish(nodes)
    nodes[0] = start
    return nodes


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
