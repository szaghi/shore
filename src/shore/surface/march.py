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
Open curves march their end nodes along their own ``n x t``.
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


DEFAULT = MarchParams()


def march_surface(
    geometry: Geometry,
    curve: np.ndarray,
    steps: np.ndarray,
    *,
    periodic: bool = False,
    entity: EntityRef = BODY,
    params: MarchParams = DEFAULT,
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
        of the model diagonal.
    """
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

    out = np.empty((r.shape[0], s.size + 1, 3))
    out[:, 0] = r
    normals = np.asarray(proj.normals, dtype=np.float64)
    for j, step in enumerate(s):
        new = _layer(out[:, j], normals, float(step), periodic, params)
        landed = geometry.project(new, entity)
        out[:, j + 1] = landed.points
        normals = np.asarray(landed.normals, dtype=np.float64)
    return out


def _central(f: np.ndarray, periodic: bool) -> np.ndarray:
    """Central difference along the layer; one-sided at open ends."""
    if periodic:
        return 0.5 * (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0))
    d = np.empty_like(f)
    d[1:-1] = 0.5 * (f[2:] - f[:-2])
    d[0] = f[1] - f[0]
    d[-1] = f[-1] - f[-2]
    return d


def _dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.einsum("...j,...j->...", a, b)


def _layer(
    r: np.ndarray, normal: np.ndarray, step: float, periodic: bool, params: MarchParams
) -> np.ndarray:
    """One marching step: the new layer before its projection onto the surface."""
    n_nodes = r.shape[0]
    t = _central(r, periodic)
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
    shrink = np.linalg.norm(_central(r + dr0, periodic), axis=-1) / h
    contraction = np.clip((1.0 - shrink) / np.maximum(rho, 1e-12), 0.0, 1.0)
    boost = 1.0 + params.convergence_gain * contraction
    eps_i = params.eps_implicit * rho * boost
    eps_e = params.eps_explicit * rho * rho * boost

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
    if not periodic:
        fixed[[0, -1]] = True  # open ends march along their own n x t
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
