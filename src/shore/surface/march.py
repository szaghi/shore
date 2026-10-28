"""Hyperbolic surface marching: the collar-grid kernel (ADR 0008).

``march_surface`` grows a structured grid over a surface from a curve, one
layer per step (Chan & Steger 1992, SURGRD type). Layer ``j + 1 = r + dr``
is marched from layer ``r`` to the ``n x t`` side (``n`` the surface normal,
``t = r_xi``), with ``dr`` from the linearised system

    t . dr          + dr0 . D(dr)        = 0          (k-lines orthogonal to the layer)
    (n x t) . dr    + (dr0 x n) . D(dr)  = dA         (prescribed cell area)
    n . dr                               = 0          (stay in the tangent plane)

about the local solution ``dr0 = dist * (n x t) / |n x t|`` (``D`` the
central difference along the layer), plus implicit dissipation
``-eps_i Lap(dr)`` and explicit tangential dissipation ``eps_e T(Lap r)``.
One sparse solve per layer (SuperLU, deterministic); the new layer is then
projected onto the surface, which also gives the next layer's normals, so
the geometry is queried once per node and layer.

Stabilisation, ramped by ``rho = min(1, step / h)`` (``h`` the node spacing)
so thin layers march along the local normal:

- area smoothing: ``dA`` is blended toward its neighbours (fronts that
  collide at a concave corner share their area instead of crossing);
- a convergence sensor raises the dissipation where the pure-normal trial
  layer contracts.

The same stabilisation lengthens the steps of converging fronts: marching
north from a sphere's equator overshoots the schedule by ~3.5 % (ADR 0008).
Open curves march their end nodes along their own ``n x t``, or, with
``plane``, hold them on a symmetry plane (ADR 0012): each end is solved as
an interior node whose missing neighbour is the mirror image of its inner
neighbour across the plane, then put on the surface and the plane by
alternating projections. On mirror-symmetric bodies this is the symmetric
half of the periodic march of the mirrored closed curve.

At a crease the node's normal is made orthogonal to the front tangent
(ADR 0023). At a sharp crease, where the front folds back on itself (the
sharp trailing edge of a wing), the central difference is much shorter than
the node spacing; its tangent is then lengthened to the mean arc spacing of
the node, so the prescribed area does not turn into an outsized step
(ADR 0024). The split collar (ADR 0024) also holds the two ends of an open
curve on two different planes (``_Ends``).

At corners where the loop is clustered (ADR 0027) a layer whose step exceeds
``max_step_ratio`` x the front spacing (counting only nodes not closer than
on the starting curve, so converging fronts do not trigger it) is marched in
sub-steps; and a layer that folds against the previous one (a cell
reversed, or a front segment pointing backwards: converging fronts at a
concave corner) is smoothed along the front round the fold and re-projected
until it does not. Marches where neither happens are unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from shore.exceptions import ParameterError
from shore.geometry.base import BODY, EntityRef, Geometry


@dataclass(frozen=True)
class MarchParams:
    """Defaults of the surface marcher (not per-case tuning knobs)."""

    eps_implicit: float = 1.0  # implicit dissipation at rho = 1
    eps_explicit: float = 0.5  # explicit tangential dissipation at rho = 1
    area_smoothing: float = 0.5  # dA blend toward the neighbour mean at rho = 1
    area_iters: int = 2
    convergence_gain: float = 16.0  # extra dissipation per unit of front contraction
    distance_clip: float = 2.0  # marching distance kept in [step / clip, step * clip]
    crease_tolerance: float = 0.1  # |n . t| above this: n made orthogonal to t (ADR 0023)
    sharp_crease: float = 0.5  # at a crease, |t| below this x the arc spacing: t lengthened
    max_step_ratio: float = 24.0  # a step over this x the spacing is sub-stepped (ADR 0027)
    substep_ratio: float = 2.0  # ... into sub-steps of at most this x the spacing
    untangle_iterations: int = 200  # smoothing passes where a layer folds (ADR 0027)
    untangle_reach: int = 2  # nodes smoothed on each side of a folded cell
    end_tolerance: float = 1e-14  # ends on surface and plane, relative to the model diagonal
    end_iterations: int = 100


DEFAULT = MarchParams()


def march_surface(
    geometry: Geometry,
    curve: np.ndarray,
    steps: np.ndarray,
    *,
    periodic: bool = False,
    entity: EntityRef = BODY,
    params: MarchParams = DEFAULT,
    plane: tuple[np.ndarray, np.ndarray] | None = None,
) -> np.ndarray:
    """March a surface grid from *curve* over *entity* of *geometry*.

    Parameters
    ----------
    geometry :
        Any backend; only ``project`` (with normals) is used.
    curve : np.ndarray, shape (n, 3)
        The starting curve, n >= 3, on the surface. With *periodic* a closed
        loop whose last node is not repeated.
    steps : np.ndarray, shape (m,)
        Positive marching distance of each layer.
    periodic :
        Whether *curve* is closed.
    entity :
        The entity marched over (default: the whole body).
    params :
        Marcher settings.
    plane :
        ``(point, normal)``: hold both ends of an open *curve* on this plane
        with a symmetry condition (ADR 0012). The curve's ends must lie on
        it. ``None``: the ends march freely.

    Returns
    -------
    np.ndarray, shape (n, m + 1, 3)
        ``[:, 0]`` is *curve*; layer j + 1 lies on the ``n x t`` side of
        layer j, every node on the surface.

    Raises
    ------
    ParameterError
        *curve* not (n, 3) with n >= 3; *steps* not a non-empty 1-D array of
        positive finite values; *curve* off the surface by more than 1e-6
        of the model diagonal; *plane* not a point and a non-zero normal,
        given with *periodic*, or with the curve's ends off it by more than
        that tolerance.
    """
    r, s, proj, tol = _checked(geometry, curve, steps, entity)
    mirror = None
    if plane is not None:
        if periodic:
            raise ParameterError("plane holds the ends of an open curve: not with periodic=True.")
        mirror = _Mirror.of(plane)
        off = float(np.max(np.abs(mirror.offset(r[[0, -1]]))))
        if off > tol:
            raise ParameterError(
                f"curve ends are off the plane by up to {off:.3e} (tolerance {tol:.3e})."
            )
    return _march(geometry, r, s, proj, tol, periodic, entity, params, mirror)


def _march_between(
    geometry: Geometry,
    curve: np.ndarray,
    steps: np.ndarray,
    ends: _Ends,
    *,
    entity: EntityRef = BODY,
    params: MarchParams = DEFAULT,
) -> np.ndarray:
    """``march_surface`` of an open *curve* whose two ends are held on two planes."""
    r, s, proj, tol = _checked(geometry, curve, steps, entity)
    for end in (0, -1):
        off = abs(float(ends.end(end).offset(r[end])))
        if off > tol:
            raise ParameterError(f"curve end is off its plane by {off:.3e} (tolerance {tol:.3e}).")
    return _march(geometry, r, s, proj, tol, False, entity, params, ends)


def _checked(
    geometry: Geometry, curve: np.ndarray, steps: np.ndarray, entity: EntityRef
) -> tuple[np.ndarray, np.ndarray, object, float]:
    """The curve and steps as arrays, the curve's projection and the surface tolerance."""
    r = np.asarray(curve, dtype=np.float64)
    if r.ndim != 2 or r.shape[1] != 3 or r.shape[0] < 3:
        raise ParameterError(f"curve must have shape (n, 3) with n >= 3, got {r.shape}.")
    s = np.asarray(steps, dtype=np.float64)
    if s.ndim != 1 or s.size == 0 or not np.all(np.isfinite(s)) or np.any(s <= 0.0):
        raise ParameterError("steps must be a non-empty 1-D array of positive finite values.")

    proj = geometry.project(r, entity)
    bounds = np.asarray(geometry.bounds, dtype=np.float64)
    tol = 1e-6 * float(np.linalg.norm(bounds[1] - bounds[0]))
    worst = float(np.max(proj.distances))
    if worst > tol:
        raise ParameterError(
            f"curve nodes are off the surface by up to {worst:.3e} (tolerance {tol:.3e})."
        )
    return r, s, proj, tol


def _march(
    geometry: Geometry,
    r: np.ndarray,
    s: np.ndarray,
    proj: object,
    tol: float,
    periodic: bool,
    entity: EntityRef,
    params: MarchParams,
    mirror: _Mirror | _Ends | None,
) -> np.ndarray:
    """The march of a checked curve *r* (its projection *proj*) by steps *s*."""
    end_tol = params.end_tolerance * tol / 1e-6
    out = np.empty((r.shape[0], s.size + 1, 3))
    out[:, 0] = r
    normals = np.asarray(proj.normals, dtype=np.float64)  # type: ignore[attr-defined]
    for j, step in enumerate(s):
        points = out[:, j]
        count = _substeps(points, out[:, 0], float(step), periodic, mirror, params)
        for _ in range(count):  # sub-steps where the step dwarfs the spacing (ADR 0027)
            new = _layer(points, normals, float(step) / count, periodic, params, mirror)
            landed = geometry.project(new, entity)
            new_points = np.array(landed.points, dtype=np.float64)
            new_normals = np.array(landed.normals, dtype=np.float64)
            new_points, new_normals = _untangle(
                geometry, entity, points, new_points, new_normals, periodic, params
            )
            points, normals = new_points, new_normals
            if mirror is not None:
                for end in (0, -1):
                    plane, its = mirror.end(end), params.end_iterations
                    points[end], normals[end] = _on_surface_and_plane(
                        geometry, entity, points[end], plane, end_tol, its
                    )
        out[:, j + 1] = points
    return out


def _inverted(old: np.ndarray, new: np.ndarray, normals: np.ndarray, periodic: bool) -> np.ndarray:
    """Folded cells (j, j + 1) between two layers (ADR 0027): oriented against the
    majority, or with the new front segment pointing backwards along the old one."""
    if periodic:
        ahead = np.roll(old, -1, axis=0) - old
        new_ahead = np.roll(new, -1, axis=0) - new
        move = new - old
        move = 0.5 * (move + np.roll(move, -1, axis=0))
        nrm = 0.5 * (normals + np.roll(normals, -1, axis=0))
    else:
        ahead, new_ahead = np.diff(old, axis=0), np.diff(new, axis=0)
        move = 0.5 * ((new - old)[:-1] + (new - old)[1:])
        nrm = 0.5 * (normals[:-1] + normals[1:])
    orient = _dot(np.cross(ahead, move), nrm)
    return (orient * np.sign(np.median(orient)) <= 0.0) | (_dot(ahead, new_ahead) <= 0.0)


def _untangle(
    geometry: Geometry,
    entity: EntityRef,
    old: np.ndarray,
    new: np.ndarray,
    normals: np.ndarray,
    periodic: bool,
    params: MarchParams,
) -> tuple[np.ndarray, np.ndarray]:
    """The new layer and its normals, smoothed along the front round any cell that folds
    against *old* (converging fronts, ADR 0027) and re-projected, until none folds.

    A layer without a folded cell is returned as it is."""
    n = len(new)
    reach = np.arange(-params.untangle_reach, params.untangle_reach + 2)
    for _ in range(params.untangle_iterations):
        bad = _inverted(old, new, normals, periodic)
        if not bad.any():
            break
        nodes = np.unique((np.flatnonzero(bad)[:, None] + reach[None, :]).ravel())
        nodes = nodes % n if periodic else nodes[(nodes > 0) & (nodes < n - 1)]
        moved = new.copy()
        moved[nodes] = 0.5 * new[nodes] + 0.25 * (new[(nodes - 1) % n] + new[(nodes + 1) % n])
        landed = geometry.project(moved, entity)
        new = np.array(landed.points, dtype=np.float64)
        normals = np.array(landed.normals, dtype=np.float64)
    return new, normals


def _substeps(
    r: np.ndarray,
    start: np.ndarray,
    step: float,
    periodic: bool,
    mirror: _Mirror | _Ends | None,
    params: MarchParams,
) -> int:
    """Sub-steps of the layer from *r* (ADR 0027): one while *step* is at most
    ``max_step_ratio`` x the smallest spacing of the nodes not closer together than on
    the starting curve *start* (converging fronts are left alone); else sub-steps of at
    most ``substep_ratio`` x that spacing."""
    h = _arc_spacing(r, periodic, mirror)
    h = h[h >= 0.9 * _arc_spacing(start, periodic, mirror)]
    if h.size == 0:
        return 1
    ratio = step / float(np.min(h))
    if ratio <= params.max_step_ratio:
        return 1
    return int(np.ceil(ratio / params.substep_ratio))


@dataclass(frozen=True)
class _Mirror:
    """Reflection across the plane through *point* with unit *normal*."""

    point: np.ndarray
    normal: np.ndarray

    @classmethod
    def of(cls, plane: tuple[np.ndarray, np.ndarray]) -> _Mirror:
        try:
            point, normal = (np.asarray(v, dtype=np.float64) for v in plane)
        except (TypeError, ValueError):
            raise ParameterError("plane must be (point, normal), two 3-vectors.") from None
        if point.shape != (3,) or normal.shape != (3,):
            raise ParameterError("plane must be (point, normal), two 3-vectors.")
        length = float(np.linalg.norm(normal))
        if not (np.all(np.isfinite(point)) and np.isfinite(length) and length > 0.0):
            raise ParameterError("plane must have a finite point and a non-zero normal.")
        return cls(point, normal / length)

    @property
    def matrix(self) -> np.ndarray:
        return np.eye(3) - 2.0 * np.outer(self.normal, self.normal)

    def offset(self, x: np.ndarray) -> np.ndarray:
        return (x - self.point) @ self.normal

    def reflect(self, x: np.ndarray) -> np.ndarray:
        return x - 2.0 * self.offset(x)[..., None] * self.normal

    def onto(self, x: np.ndarray) -> np.ndarray:
        return x - self.offset(x)[..., None] * self.normal

    def extend(self, f: np.ndarray) -> np.ndarray:
        """*f* with a ghost node at each end: the mirror image of the inner neighbour."""
        return np.concatenate([self.reflect(f[1:2]), f, self.reflect(f[-2:-1])])

    def end(self, which: int) -> _Mirror:
        """The plane of end *which* (0 or -1): this one at both ends."""
        return self

    def matrices(self, first: np.ndarray) -> np.ndarray:
        """Reflection matrices of the ghosts past the first (*first*) or last end."""
        return np.broadcast_to(self.matrix, (first.size, 3, 3))


@dataclass(frozen=True)
class _Ends:
    """Two planes, one for each end of an open curve (ADR 0024)."""

    first: _Mirror
    last: _Mirror

    def end(self, which: int) -> _Mirror:
        return self.first if which == 0 else self.last

    def extend(self, f: np.ndarray) -> np.ndarray:
        return np.concatenate([self.first.reflect(f[1:2]), f, self.last.reflect(f[-2:-1])])

    def matrices(self, first: np.ndarray) -> np.ndarray:
        return np.where(first[:, None, None], self.first.matrix, self.last.matrix)


def _on_surface_and_plane(
    geometry: Geometry,
    entity: EntityRef,
    x: np.ndarray,
    mirror: _Mirror,
    tol: float,
    iterations: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Alternating projections onto the plane and the surface; the point and its normal."""
    for _ in range(iterations):
        landed = geometry.project(mirror.onto(x)[None], entity)
        y = landed.points[0]
        moved = float(np.linalg.norm(y - x))
        x = y
        if moved <= tol:
            break
    return x, np.asarray(landed.normals[0], dtype=np.float64)


def _central(f: np.ndarray, periodic: bool, mirror: _Mirror | _Ends | None = None) -> np.ndarray:
    """Central difference along the layer; one-sided at free ends, mirrored at a plane."""
    if mirror is not None:
        ghosted = mirror.extend(f)
        return 0.5 * (ghosted[2:] - ghosted[:-2])
    if periodic:
        return 0.5 * (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0))
    d = np.empty_like(f)
    d[1:-1] = 0.5 * (f[2:] - f[:-2])
    d[0] = f[1] - f[0]
    d[-1] = f[-1] - f[-2]
    return d


def _arc_spacing(r: np.ndarray, periodic: bool, mirror: _Mirror | _Ends | None) -> np.ndarray:
    """Mean length of the two segments at each node (one segment at free ends)."""
    if mirror is not None:
        g = mirror.extend(r)
        ahead, behind = g[2:] - r, r - g[:-2]
    elif periodic:
        ahead, behind = np.roll(r, -1, axis=0) - r, r - np.roll(r, 1, axis=0)
    else:
        seg = np.diff(r, axis=0)
        ahead, behind = np.vstack([seg, seg[-1:]]), np.vstack([seg[:1], seg])
    return 0.5 * (np.linalg.norm(ahead, axis=-1) + np.linalg.norm(behind, axis=-1))


def _dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.einsum("...j,...j->...", a, b)


def _layer(
    r: np.ndarray,
    normal: np.ndarray,
    step: float,
    periodic: bool,
    params: MarchParams,
    mirror: _Mirror | _Ends | None = None,
) -> np.ndarray:
    """One marching step: the new layer before its projection onto the surface."""
    n_nodes = r.shape[0]
    t = _central(r, periodic, mirror)
    h = np.linalg.norm(t, axis=-1)
    # at a crease (a corner of the front on an edge of the surface) the node's normal is
    # one face's and leans along the front tangent; the scheme assumes n orthogonal to t
    along = _dot(normal, t) / h
    crease = np.abs(along) > params.crease_tolerance
    if crease.any():
        normal = normal.copy()
        flat = normal[crease] - (along[crease] / h[crease])[:, None] * t[crease]
        normal[crease] = flat / np.linalg.norm(flat, axis=-1, keepdims=True)
        # a sharp crease folds the front back: |t| << the spacing (ADR 0024)
        arc = _arc_spacing(r, periodic, mirror)
        sharp = crease & (h < params.sharp_crease * arc)
        if sharp.any():
            t = t.copy()
            t[sharp] *= (arc[sharp] / h[sharp])[:, None]
            h = np.linalg.norm(t, axis=-1)
    b = np.cross(normal, t)
    b /= np.linalg.norm(b, axis=-1, keepdims=True)
    rho = np.minimum(1.0, step / h)

    # prescribed area, smoothed toward the neighbours
    area = step * h
    for _ in range(params.area_iters):
        neighbours = 0.5 * (np.roll(area, 1) + np.roll(area, -1))
        if not periodic:
            neighbours[0], neighbours[-1] = area[1], area[-2]
        nu = params.area_smoothing * rho
        area = (1.0 - nu) * area + nu * neighbours
    clip = params.distance_clip
    dist = np.clip(area / h, step / clip, step * clip)
    dr0 = dist[:, None] * b

    # front-contraction sensor on the pure-normal trial layer
    shrink = np.linalg.norm(_central(r + dr0, periodic, mirror), axis=-1) / h
    contraction = np.clip((1.0 - shrink) / np.maximum(rho, 1e-12), 0.0, 1.0)
    boost = 1.0 + params.convergence_gain * contraction
    eps_i = params.eps_implicit * rho * boost
    eps_e = params.eps_explicit * rho * rho * boost

    if mirror is not None:
        ghosted = mirror.extend(r)
        lap = ghosted[2:] + ghosted[:-2] - 2.0 * r
    else:
        lap = np.roll(r, -1, axis=0) + np.roll(r, 1, axis=0) - 2.0 * r
        if not periodic:
            lap[0] = lap[-1] = 0.0
    lap_t = lap - _dot(lap, normal)[:, None] * normal

    # C dr + A D(dr) - eps_i C Lap(dr) = [0, dist h, 0] + C (eps_e lap_t)
    c_mat = np.stack([t, np.cross(normal, t), normal], axis=1)  # (n, 3, 3)
    a_mat = np.stack([dr0, np.cross(dr0, normal), np.zeros_like(dr0)], axis=1)
    rhs = np.zeros((n_nodes, 3))
    rhs[:, 1] = dist * h
    rhs += np.einsum("nab,nb->na", c_mat, eps_e[:, None] * lap_t)

    fixed = np.zeros(n_nodes, dtype=bool)
    if not periodic and mirror is None:
        fixed[[0, -1]] = True  # free ends march along their own n x t
    unknown = np.flatnonzero(~fixed)
    index = -np.ones(n_nodes, dtype=np.int64)
    index[unknown] = np.arange(unknown.size)

    rows, cols, vals = [], [], []

    def block(i: np.ndarray, j: np.ndarray, mats: np.ndarray) -> None:
        r3 = (3 * i[:, None, None] + np.arange(3)[None, :, None]).repeat(3, axis=2)
        c3 = (3 * j[:, None, None] + np.arange(3)[None, None, :]).repeat(3, axis=1)
        rows.append(r3.ravel())
        cols.append(c3.ravel())
        vals.append(mats.ravel())

    u = np.arange(unknown.size)
    block(u, u, c_mat[unknown] * (1.0 + 2.0 * eps_i[unknown])[:, None, None])
    sys_rhs = rhs[unknown].copy()
    for shift, sign in ((1, 0.5), (-1, -0.5)):
        nb = (unknown + shift) % n_nodes
        mats = sign * a_mat[unknown] - eps_i[unknown, None, None] * c_mat[unknown]
        if mirror is not None:  # past an end: the ghost moves as the mirrored inner neighbour
            past = (unknown + shift < 0) | (unknown + shift >= n_nodes)
            nb = np.where(past, unknown - shift, nb)
            mats = np.where(past[:, None, None], mats @ mirror.matrices(unknown + shift < 0), mats)
        free = ~fixed[nb]
        block(u[free], index[nb[free]], mats[free])
        if (~free).any():  # neighbour is a fixed end: known increment
            sys_rhs[~free] -= np.einsum("nab,nb->na", mats[~free], dr0[nb[~free]])

    matrix = sp.csc_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(3 * unknown.size, 3 * unknown.size),
    )
    dr = dr0.copy()
    dr[unknown] = spla.splu(matrix).solve(sys_rhs.ravel()).reshape(-1, 3)
    return r + dr
