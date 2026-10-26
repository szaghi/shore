"""Hyperbolic volume marching: one layer of the Steger-Chaussee scheme (ADR 0005).

The new layer ``r_{k+1} = r_k + dr`` satisfies, at every node,

    r_xi . r_zeta = 0,   r_eta . r_zeta = 0,   det[r_xi, r_eta, r_zeta] = dV,

(the k-lines are orthogonal to the layer, and the cell volume is prescribed).
Linearised about the local solution ``r_zeta0 = dV * n / |n|^2`` with
``n = r_xi x r_eta`` (Steger & Chaussee 1980; Chan & Steger 1992), the system
for ``dr`` becomes

    (I + Ci A d_xi + Ci B d_eta - eps_i Lap) dr = r_zeta0 + eps_e T(Lap r_k)

with ``A = [r_zeta0; 0; r_eta x r_zeta0]``, ``B = [0; r_zeta0; r_zeta0 x r_xi]``,
``Ci = C^-1``, ``C = [r_xi; r_eta; n]``, central differences ``d``, the
5-point Laplacian ``Lap`` and ``T`` the projection onto the layer's tangent
plane. The coupling ``Ci A``, ``Ci B`` scales with step / cell size: thin
near-wall layers march almost exactly along the local normal, thick layers
feel their neighbours, which spreads converging fronts at concave corners.

Stabilisation, all scaled by ``rho = min(1, step / h)`` (``h`` the local
cell size), so the near-wall layers stay orthogonal and of exact thickness:

- volume smoothing: ``dV`` is blended toward its neighbour average;
- a grid-convergence sensor raises dissipation where neighbouring nodes of
  the local-solution layer approach each other (concave regions);
- implicit (``eps_i``) and tangential explicit (``eps_e``) second-difference
  dissipation.

Boundary rows (non-periodic edges) are Dirichlet: they take the local
solution (or the caller's *boundary_target*, for Dirichlet-pinned cap
edges). Seam continuity is the topology marcher's job, as for the algebraic
scheme. The linear system is solved per block layer with SuperLU
(deterministic).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from shore.volume.hyperbolic import _compute_normals


@dataclass(frozen=True)
class HyperbolicParams:
    """Defaults of the hyperbolic scheme (not per-case tuning knobs)."""

    eps_implicit: float = 0.5  # implicit dissipation at rho = 1
    eps_explicit: float = 0.25  # explicit tangential dissipation at rho = 1
    volume_smoothing: float = 0.2  # dV blend toward neighbour mean at rho = 1
    volume_iters: int = 2
    convergence_gain: float = 16.0  # extra dissipation per unit of front convergence
    distance_clip: float = 2.0  # marching distance kept in [step / clip, step * clip]
    excess_bound: float = 0.04  # outer-face volume correction, at most this fraction of dV


DEFAULT = HyperbolicParams()


def _tangents(
    layer: np.ndarray, periodic: tuple[bool, bool], ghost: dict[str, np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    """Derivatives r_xi, r_eta (unit index spacing), central where possible."""
    pi, pj = periodic
    if pi:
        ti = 0.5 * (np.roll(layer, -1, axis=0) - np.roll(layer, 1, axis=0))
    else:
        ti = np.empty_like(layer)
        ti[1:-1] = 0.5 * (layer[2:] - layer[:-2])
        ti[0] = 0.5 * (layer[1] - ghost["i_lo"]) if "i_lo" in ghost else layer[1] - layer[0]
        ti[-1] = 0.5 * (ghost["i_hi"] - layer[-2]) if "i_hi" in ghost else layer[-1] - layer[-2]
    if pj:
        tj = 0.5 * (np.roll(layer, -1, axis=1) - np.roll(layer, 1, axis=1))
    else:
        tj = np.empty_like(layer)
        tj[:, 1:-1] = 0.5 * (layer[:, 2:] - layer[:, :-2])
        tj[:, 0] = (
            0.5 * (layer[:, 1] - ghost["j_lo"]) if "j_lo" in ghost else layer[:, 1] - layer[:, 0]
        )
        tj[:, -1] = (
            0.5 * (ghost["j_hi"] - layer[:, -2]) if "j_hi" in ghost else layer[:, -1] - layer[:, -2]
        )
    return ti, tj


def _neighbour_mean(field: np.ndarray, periodic: tuple[bool, bool]) -> np.ndarray:
    """Mean of the available 4-neighbours (edges use what exists)."""
    total = np.zeros_like(field)
    count = np.zeros(field.shape[:2] + (1,) * (field.ndim - 2))
    for axis, per in ((0, periodic[0]), (1, periodic[1])):
        for shift in (1, -1):
            nb = np.roll(field, shift, axis=axis)
            valid = np.ones(field.shape[:2], dtype=bool)
            if not per:
                edge = [slice(None)] * 2
                edge[axis] = 0 if shift == 1 else -1
                valid[tuple(edge)] = False
            v = valid.reshape(valid.shape + (1,) * (field.ndim - 2))
            total += np.where(v, nb, 0.0)
            count = count + v
    return total / np.maximum(count, 1)


def _central(field: np.ndarray, axis: int, periodic: bool) -> np.ndarray:
    """Half the centred difference along *axis*; past a non-periodic edge the
    neighbour is the node itself (the system's seam extrapolation)."""
    up = np.roll(field, -1, axis)
    down = np.roll(field, 1, axis)
    if not periodic:
        first = [slice(None)] * 2
        last = [slice(None)] * 2
        first[axis], last[axis] = 0, -1
        down[tuple(first)] = field[tuple(first)]
        up[tuple(last)] = field[tuple(last)]
    return 0.5 * (up - down)


def _laplacian(
    field: np.ndarray, periodic: tuple[bool, bool], ghost: dict[str, np.ndarray]
) -> np.ndarray:
    """5-point Laplacian (unit spacing). A non-periodic edge uses its ghost row
    when one is given; otherwise that direction's term is zero there."""
    lap = np.zeros_like(field)
    for axis, (lo, hi) in ((0, ("i_lo", "i_hi")), (1, ("j_lo", "j_hi"))):
        term = np.roll(field, 1, axis) + np.roll(field, -1, axis) - 2.0 * field
        if not periodic[axis]:
            first = [slice(None)] * 2
            last = [slice(None)] * 2
            first[axis], last[axis] = 0, -1
            second = [slice(None)] * 2
            before = [slice(None)] * 2
            second[axis], before[axis] = 1, -2
            f0, fl = tuple(first), tuple(last)
            term[f0] = field[tuple(second)] + ghost[lo] - 2.0 * field[f0] if lo in ghost else 0.0
            term[fl] = ghost[hi] + field[tuple(before)] - 2.0 * field[fl] if hi in ghost else 0.0
        lap += term
    return lap


def hyperbolic_layer(
    layer: np.ndarray,
    step: float,
    periodic: tuple[bool, bool],
    ghost_rows: dict[str, np.ndarray] | None = None,
    normals_override: np.ndarray | None = None,
    boundary_target: np.ndarray | None = None,
    params: HyperbolicParams = DEFAULT,
) -> np.ndarray:
    """Advance *layer* (ni, nj, 3) by one hyperbolic step of size *step*.

    Parameters
    ----------
    layer :
        Current layer.
    step :
        Marching distance of this layer.
    periodic :
        ``(periodic_i, periodic_j)``.
    ghost_rows :
        Neighbour rows across block seams (see ``_compute_normals``).
    normals_override :
        Replaces the local marching direction (normal blending).
    boundary_target :
        Full (ni, nj, 3) array whose non-periodic edge rows are the required
        new boundary positions (Dirichlet-pinned cap edges); interior values
        are ignored.
    """
    ghost = ghost_rows or {}
    ni, nj, _ = layer.shape
    ti, tj = _tangents(layer, periodic, ghost)
    n = np.cross(ti, tj)
    area = np.linalg.norm(n, axis=-1)
    nhat = (
        normals_override
        if normals_override is not None
        else _compute_normals(layer, periodic=periodic, ghost_rows=ghost_rows)
    )

    # local cell size and the stabilisation ramp rho = min(1, step / h)
    h = np.sqrt(np.maximum(area, 1e-300))
    rho = np.minimum(1.0, abs(step) / h)

    # prescribed volume, blended toward its neighbours away from the wall
    dv = step * np.where(area > 0.0, area, 0.0)
    for _ in range(params.volume_iters):
        nu = params.volume_smoothing * rho
        dv = (1.0 - nu) * dv + nu * _neighbour_mean(dv, periodic)
    dist = np.where(area > 1e-300, dv / np.where(area > 0.0, area, 1.0), step)
    # bound it: small cells next to large ones (sharp convex edges) would
    # otherwise inherit their neighbours' volume as a large marching distance
    dist = np.clip(dist, abs(step) / params.distance_clip, abs(step) * params.distance_clip)
    rz0 = dist[..., None] * nhat  # local solution r_zeta0

    # grid-convergence sensor: neighbour spacing of the local-solution layer
    # relative to the current one (< 1 where fronts converge)
    trial = layer + rz0
    t_ti, t_tj = _tangents(trial, periodic, {})
    shrink = np.minimum(
        np.linalg.norm(t_ti, axis=-1) / np.maximum(np.linalg.norm(ti, axis=-1), 1e-300),
        np.linalg.norm(t_tj, axis=-1) / np.maximum(np.linalg.norm(tj, axis=-1), 1e-300),
    )
    converge = np.clip(1.0 - shrink, 0.0, 1.0) / np.maximum(rho, 1e-12)
    boost = 1.0 + params.convergence_gain * np.minimum(converge, 1.0)
    eps_i = params.eps_implicit * rho * boost
    eps_e = params.eps_explicit * rho * rho * boost  # tilt ~ eps * rho * |lap|/h

    # explicit tangential dissipation of the current layer. Not across a seam:
    # grid lines kink where blocks meet (on the cubed sphere, the great-circle
    # grid lines bend at every cube edge), a topological feature that the
    # Laplacian would read as roughness and straighten by sliding the seam
    # nodes sideways, tilting the k-lines off the wall normal (7-8 deg)
    lap = _laplacian(layer, periodic, {})
    lap_t = lap - np.einsum("...j,...j->...", lap, nhat)[..., None] * nhat
    rhs_field = rz0 + eps_e[..., None] * lap_t

    # Edge rows are unknowns where a neighbour block's ghost row exists (the
    # neighbour's increment is extrapolated from the local solution), and
    # fixed (Dirichlet) where there is none or the caller pins them.
    dr = rz0.copy()
    fixed = np.zeros((ni, nj), dtype=bool)
    for axis, (lo, hi) in ((0, ("i_lo", "i_hi")), (1, ("j_lo", "j_hi"))):
        if periodic[axis]:
            continue
        for side, key in ((0, lo), (-1, hi)):
            sl = [slice(None)] * 2
            sl[axis] = side
            if boundary_target is not None or key not in ghost:
                fixed[tuple(sl)] = True
    if boundary_target is not None:
        dr[fixed] = boundary_target[fixed] - layer[fixed]
    interior = ~fixed
    if not interior.any():
        return layer + dr

    # local matrices: C = [r_xi; r_eta; n], A, B linearised about r_zeta0
    c_mat = np.stack([ti, tj, n], axis=-2)
    a_mat = np.stack([rz0, np.zeros_like(rz0), np.cross(tj, rz0)], axis=-2)
    b_mat = np.stack([np.zeros_like(rz0), rz0, np.cross(rz0, ti)], axis=-2)
    det = np.linalg.det(c_mat)
    safe = np.abs(det) > 1e-300
    c_inv = np.zeros_like(c_mat)
    c_inv[safe] = np.linalg.inv(c_mat[safe])
    ma = c_inv @ a_mat  # (ni, nj, 3, 3)
    mb = c_inv @ b_mat

    # The volume row holds over the outer face, but r_zeta0 = dist * n sets
    # the inner face's volume: on a convex wall each layer would march only
    # ~dist * (1 - 2 * dist * curvature). Add the outer-face excess, with the
    # system's own differences, so that dr = r_zeta0 solves the undissipated
    # system. Only where the curvature is resolved: the excess is bounded to
    # excess_bound * dV, and none is added on concave layers. At a convex kink
    # (excess ~ dV) and on converging fronts the shorter or longer step of the
    # inner-face target is what keeps the cells from folding.
    excess = np.einsum("...j,...j->...", a_mat[..., 2, :], _central(dr, 0, periodic[0]))
    excess += np.einsum("...j,...j->...", b_mat[..., 2, :], _central(dr, 1, periodic[1]))
    excess = np.clip(excess, 0.0, params.excess_bound * dist * area)
    rhs_field = rhs_field + c_inv[..., :, 2] * excess[..., None]

    idx = -np.ones((ni, nj), dtype=np.int64)
    ii, jj = np.nonzero(interior)
    idx[ii, jj] = np.arange(len(ii))
    m = len(ii)
    eye = np.eye(3)
    rows, cols, vals = [], [], []
    rhs = rhs_field[ii, jj].copy()

    def couple(blocks: np.ndarray, di: int, dj: int) -> None:
        ni_ = ii + di
        nj_ = jj + dj
        if periodic[0]:
            ni_ %= ni
        if periodic[1]:
            nj_ %= nj
        inside = (ni_ >= 0) & (ni_ < ni) & (nj_ >= 0) & (nj_ < nj)
        # across a seam (ghost row): the neighbour's increment ~ this node's
        # local solution (zero-gradient extrapolation), a known value
        po = np.nonzero(~inside)[0]
        if len(po):
            rhs[po] -= np.einsum("pab,pb->pa", blocks[po], rz0[ii[po], jj[po]])
        p = np.nonzero(inside)[0]
        q = idx[ni_[p], nj_[p]]
        unknown = q >= 0
        # neighbour is an unknown: matrix entry
        pu, qu = p[unknown], q[unknown]
        if len(pu):
            blk = blocks[pu]
            r = (3 * pu[:, None, None] + np.arange(3)[None, :, None]).repeat(3, axis=2)
            c = (3 * qu[:, None, None] + np.arange(3)[None, None, :]).repeat(3, axis=1)
            rows.append(r.ravel())
            cols.append(c.ravel())
            vals.append(blk.ravel())
        # neighbour is a fixed boundary node: move to the right-hand side
        pf = p[~unknown]
        if len(pf):
            fixed = dr[ni_[pf], nj_[pf]]
            rhs[pf] -= np.einsum("pab,pb->pa", blocks[pf], fixed)

    e_i = eps_i[ii, jj][:, None, None] * eye
    for di, dj, mat, sign in (
        (1, 0, ma, 0.5),
        (-1, 0, ma, -0.5),
        (0, 1, mb, 0.5),
        (0, -1, mb, -0.5),
    ):
        couple(sign * mat[ii, jj] - e_i, di, dj)
    diag = eye + 4.0 * e_i
    r = (3 * np.arange(m)[:, None, None] + np.arange(3)[None, :, None]).repeat(3, axis=2)
    c = (3 * np.arange(m)[:, None, None] + np.arange(3)[None, None, :]).repeat(3, axis=1)
    rows.append(r.ravel())
    cols.append(c.ravel())
    vals.append(diag.ravel())

    matrix = sp.csr_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * m, 3 * m)
    )
    solution = spla.spsolve(matrix.tocsc(), rhs.ravel())
    dr[ii, jj] = solution.reshape(m, 3)
    return layer + dr
