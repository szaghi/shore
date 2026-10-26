"""Surface grids on one parametric face: TFI and elliptic smoothing in (u, v) (ADR 0007).

``surface_grid`` fills a face from its four boundary rows:

1. the boundary nodes are projected once onto the SURFACE entity, which
   gives their parameters (u, v);
2. transfinite interpolation (TFI) of those parameters;
3. with ``smoothing="elliptic"``, Winslow smoothing in (u, v):

       alpha (r_xixi + P r_xi) - 2 beta r_xieta + gamma (r_etaeta + Q r_eta) + c = 0

   for ``r = (u, v)``. ``alpha = |x_eta|^2``, ``beta = x_xi . x_eta`` and
   ``gamma = |x_xi|^2`` are physical (``x_xi = S_u u_xi + S_v v_xi``), and
   ``c = G^-1 [S_u, S_v]^T (alpha S(r_xi, r_xi) - 2 beta S(r_xi, r_eta) +
   gamma S(r_eta, r_eta))`` holds the parametrisation's second derivatives
   (``S(p, w) = S_uu p_u w_u + S_uv (p_u w_v + p_v w_u) + S_vv p_v w_v``,
   ``G`` the first fundamental form): the physical Winslow equations
   projected on the tangent plane, so the result does not depend on how
   the face is parametrised. ``P`` and ``Q`` are Thomas-Middlecoff control
   functions, ``-(x_xi . x_xixi) / |x_xi|^2`` along the two j-boundaries
   (``P``) and the two i-boundaries (``Q``) of the physical boundary,
   blended linearly across the face: on tensor-product boundary data they
   make the TFI grid an exact solution, so the boundary clustering is kept;
4. one evaluation on the surface at the end; the boundary rows are the
   input.

The surface metric comes from a parameter lattice evaluated once (its
finite-difference derivatives are interpolated bilinearly), so the number
of geometry queries is bounded per node, never per sweep (ADR 0002). The
nonlinear system is solved by Picard iteration, one SuperLU factorisation
per outer step for both u and v (deterministic).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from shore.exceptions import ParameterError
from shore.geometry.base import EntityKind, EntityRef, Geometry

SMOOTHING = ("none", "elliptic")


@dataclass(frozen=True)
class EllipticParams:
    """Defaults of the elliptic smoother (not per-case tuning knobs)."""

    lattice_refine: int = 2  # metric lattice spacing = grid parameter spacing / refine
    max_iterations: int = 400
    tolerance: float = 1e-12  # max parameter update, relative to the (u, v) extent
    control_limit: float = 1.5  # |P|, |Q| bound: keeps the 5-point operator an M-matrix


DEFAULT = EllipticParams()


@dataclass(frozen=True)
class SurfaceGrid:
    """A structured grid on one face.

    Parameters
    ----------
    points : np.ndarray, shape (ni, nj, 3)
        Node positions, on the surface; the boundary rows are the input.
    params : np.ndarray, shape (ni, nj, 2)
        The nodes' (u, v) on the surface.
    """

    points: np.ndarray
    params: np.ndarray


def surface_grid(
    geometry: Geometry,
    surface: EntityRef,
    boundary: np.ndarray,
    smoothing: str = "elliptic",
    params: EllipticParams = DEFAULT,
) -> SurfaceGrid:
    """Fill the face of *geometry*'s SURFACE entity *surface* from *boundary*.

    Parameters
    ----------
    geometry :
        A backend with parametric SURFACE entities (``OcctGeometry``).
    surface :
        The SURFACE entity the face lies on.
    boundary : np.ndarray, shape (ni, nj, 3)
        Only its four boundary rows are read; they must lie on the surface.
    smoothing :
        ``"none"`` (plain TFI of the parameters) or ``"elliptic"``.
    params :
        Elliptic smoother settings.

    Returns
    -------
    SurfaceGrid

    Raises
    ------
    ParameterError
        Unknown *smoothing*; *boundary* not ``(ni, nj, 3)`` with ni, nj >= 3;
        boundary nodes off the surface (beyond 1e-6 of the model diagonal);
        an entity without (u, v) parameters.
    """
    if smoothing not in SMOOTHING:
        raise ParameterError(f"smoothing must be one of {SMOOTHING}, got {smoothing!r}.")
    b = np.asarray(boundary, dtype=np.float64)
    if b.ndim != 3 or b.shape[2] != 3 or b.shape[0] < 3 or b.shape[1] < 3:
        raise ParameterError(
            f"boundary must have shape (ni, nj, 3) with ni, nj >= 3, got {b.shape}."
        )
    if not isinstance(surface, EntityRef) or surface.kind is not EntityKind.SURFACE:
        raise ParameterError(f"surface must be a SURFACE entity with (u, v), got {surface!r}.")
    ni, nj, _ = b.shape

    rows = _boundary_mask(ni, nj)
    proj = geometry.project(b[rows], surface, exact=True)
    if proj.params is None:
        raise ParameterError(f"{surface!r} has no (u, v) parametrisation.")
    bounds = np.asarray(geometry.bounds, dtype=np.float64)
    tol = 1e-6 * float(np.linalg.norm(bounds[1] - bounds[0]))
    worst = float(np.max(proj.distances))
    if worst > tol:
        raise ParameterError(
            f"boundary nodes are off the surface by up to {worst:.3e} (tolerance {tol:.3e})."
        )

    uv = np.zeros((ni, nj, 2))
    uv[rows] = proj.params
    uv = _tfi(uv)
    if smoothing == "elliptic":
        uv = _winslow(geometry, surface, b, uv, params)

    points = np.asarray(geometry.evaluate(surface, uv), dtype=np.float64)
    points[rows] = b[rows]
    return SurfaceGrid(points=points, params=uv)


# ---------------------------------------------------------------------------
# TFI
# ---------------------------------------------------------------------------


def _boundary_mask(ni: int, nj: int) -> np.ndarray:
    mask = np.zeros((ni, nj), dtype=bool)
    mask[0] = mask[-1] = mask[:, 0] = mask[:, -1] = True
    return mask


def _tfi(f: np.ndarray) -> np.ndarray:
    """Linear transfinite interpolation of *f*'s boundary rows (index-uniform blending)."""
    ni, nj = f.shape[:2]
    s = np.linspace(0.0, 1.0, ni)[:, None, None]
    t = np.linspace(0.0, 1.0, nj)[None, :, None]
    out = (
        (1 - s) * f[0][None] + s * f[-1][None] + (1 - t) * f[:, 0][:, None] + t * f[:, -1][:, None]
        - ((1 - s) * (1 - t) * f[0, 0] + s * (1 - t) * f[-1, 0]
           + (1 - s) * t * f[0, -1] + s * t * f[-1, -1])
    )  # fmt: skip
    out[_boundary_mask(ni, nj)] = f[_boundary_mask(ni, nj)]
    return out


# ---------------------------------------------------------------------------
# Elliptic smoothing in (u, v)
# ---------------------------------------------------------------------------


class _Metric:
    """``S_u, S_v, S_uu, S_uv, S_vv`` from one lattice of surface evaluations."""

    def __init__(
        self, geometry: Geometry, surface: EntityRef, lo: np.ndarray, hi: np.ndarray, shape
    ) -> None:
        mu, mv = shape
        span = np.where(hi - lo > 0.0, hi - lo, 1.0)
        # a margin of one lattice cell, so nodes on the boundary interpolate inside
        lo = lo - span / (np.array(shape) - 3)
        hi = hi + span / (np.array(shape) - 3)
        self.lo, self.du = lo, (hi - lo) / (np.array(shape) - 1)
        u = lo[0] + self.du[0] * np.arange(mu)
        v = lo[1] + self.du[1] * np.arange(mv)
        uu, vv = np.meshgrid(u, v, indexing="ij")
        x = np.asarray(geometry.evaluate(surface, np.stack([uu, vv], axis=-1)), dtype=np.float64)
        hu, hv = self.du
        xu = np.gradient(x, hu, axis=0, edge_order=2)
        xv = np.gradient(x, hv, axis=1, edge_order=2)
        self.fields = np.concatenate(
            [
                xu,
                xv,
                np.gradient(xu, hu, axis=0, edge_order=2),
                np.gradient(xu, hv, axis=1, edge_order=2),
                np.gradient(xv, hv, axis=1, edge_order=2),
            ],
            axis=-1,
        )  # (mu, mv, 15)
        self.shape = np.array(shape)

    def at(self, uv: np.ndarray) -> tuple[np.ndarray, ...]:
        """Bilinear interpolation of the five derivative fields at *uv* (..., 2)."""
        g = (uv - self.lo) / self.du
        k = np.clip(np.floor(g).astype(np.int64), 0, self.shape - 2)
        w = np.clip(g - k, 0.0, 1.0)
        i, j = k[..., 0], k[..., 1]
        a, c = w[..., 0:1], w[..., 1:2]
        f = self.fields
        out = (
            (1 - a) * (1 - c) * f[i, j] + a * (1 - c) * f[i + 1, j]
            + (1 - a) * c * f[i, j + 1] + a * c * f[i + 1, j + 1]
        )  # fmt: skip
        return tuple(out[..., 3 * n : 3 * n + 3] for n in range(5))


def _control(edge: np.ndarray, limit: float) -> np.ndarray:
    """Thomas-Middlecoff ``-(x_s . x_ss) / |x_s|^2`` along a boundary row (ends: 0)."""
    d1 = np.zeros_like(edge)
    d2 = np.zeros_like(edge)
    d1[1:-1] = 0.5 * (edge[2:] - edge[:-2])
    d2[1:-1] = edge[2:] + edge[:-2] - 2.0 * edge[1:-1]
    num = np.einsum("...j,...j->...", d1, d2)
    den = np.einsum("...j,...j->...", d1, d1)
    out = np.where(den > 0.0, -num / np.where(den > 0.0, den, 1.0), 0.0)
    return np.clip(out, -limit, limit)


def _second(hessian: tuple[np.ndarray, ...], p: np.ndarray, w: np.ndarray) -> np.ndarray:
    """``S_uu p_u w_u + S_uv (p_u w_v + p_v w_u) + S_vv p_v w_v``."""
    suu, suv, svv = hessian
    return (
        suu * (p[..., :1] * w[..., :1])
        + suv * (p[..., :1] * w[..., 1:] + p[..., 1:] * w[..., :1])
        + svv * (p[..., 1:] * w[..., 1:])
    )


def _winslow(
    geometry: Geometry,
    surface: EntityRef,
    boundary: np.ndarray,
    uv: np.ndarray,
    params: EllipticParams,
) -> np.ndarray:
    ni, nj = uv.shape[:2]
    rows = _boundary_mask(ni, nj)
    lo, hi = uv[rows].min(axis=0), uv[rows].max(axis=0)
    refine = params.lattice_refine
    metric = _Metric(geometry, surface, lo, hi, (refine * (ni - 1) + 3, refine * (nj - 1) + 3))
    extent = float(np.max(hi - lo)) or 1.0

    # control functions from the physical boundary, blended linearly in index
    t = np.linspace(0.0, 1.0, nj)[None, :]
    s = np.linspace(0.0, 1.0, ni)[:, None]
    p_fn = (1 - t) * _control(boundary[:, 0], params.control_limit)[:, None]
    p_fn = p_fn + t * _control(boundary[:, -1], params.control_limit)[:, None]
    q_fn = (1 - s) * _control(boundary[0], params.control_limit)[None, :]
    q_fn = q_fn + s * _control(boundary[-1], params.control_limit)[None, :]
    p_in, q_in = p_fn[1:-1, 1:-1], q_fn[1:-1, 1:-1]

    mi, mj = ni - 2, nj - 2
    idx = np.arange(mi * mj).reshape(mi, mj)
    r = uv.copy()
    for _ in range(params.max_iterations):
        ru = r[1:-1, 1:-1]
        r_xi = 0.5 * (r[2:, 1:-1] - r[:-2, 1:-1])
        r_eta = 0.5 * (r[1:-1, 2:] - r[1:-1, :-2])
        r_xe = 0.25 * (r[2:, 2:] - r[2:, :-2] - r[:-2, 2:] + r[:-2, :-2])
        su, sv, suu, suv, svv = metric.at(ru)
        x_xi = su * r_xi[..., :1] + sv * r_xi[..., 1:]
        x_eta = su * r_eta[..., :1] + sv * r_eta[..., 1:]
        alpha = np.einsum("...j,...j->...", x_eta, x_eta)
        beta = np.einsum("...j,...j->...", x_xi, x_eta)
        gamma = np.einsum("...j,...j->...", x_xi, x_xi)

        hessian = (suu, suv, svv)
        curv = (
            alpha[..., None] * _second(hessian, r_xi, r_xi)
            - 2.0 * beta[..., None] * _second(hessian, r_xi, r_eta)
            + gamma[..., None] * _second(hessian, r_eta, r_eta)
        )
        g11 = np.einsum("...j,...j->...", su, su)
        g12 = np.einsum("...j,...j->...", su, sv)
        g22 = np.einsum("...j,...j->...", sv, sv)
        b1 = np.einsum("...j,...j->...", su, curv)
        b2 = np.einsum("...j,...j->...", sv, curv)
        det = g11 * g22 - g12 * g12
        corr = np.stack([g22 * b1 - g12 * b2, g11 * b2 - g12 * b1], axis=-1) / det[..., None]

        # alpha (r_E + r_W - 2 r) + alpha P (r_E - r_W) / 2 + gamma (...) = 2 beta r_xe - c
        a_e = alpha * (1.0 + 0.5 * p_in)
        a_w = alpha * (1.0 - 0.5 * p_in)
        a_n = gamma * (1.0 + 0.5 * q_in)
        a_s = gamma * (1.0 - 0.5 * q_in)
        rhs = 2.0 * beta[..., None] * r_xe - corr
        rows_, cols_, vals_ = [idx.ravel()], [idx.ravel()], [(-2.0 * (alpha + gamma)).ravel()]
        for coef, di, dj in ((a_e, 1, 0), (a_w, -1, 0), (a_n, 0, 1), (a_s, 0, -1)):
            ii = np.arange(mi)[:, None] + di
            jj = np.arange(mj)[None, :] + dj
            inside = (ii >= 0) & (ii < mi) & (jj >= 0) & (jj < mj)
            inside = np.broadcast_to(inside, (mi, mj))
            src = idx[inside]
            dst = idx[np.clip(ii, 0, mi - 1), np.clip(jj, 0, mj - 1)][inside]
            rows_.append(src)
            cols_.append(dst)
            vals_.append(coef[inside])
            # neighbours on the boundary: known, to the right-hand side
            edge = ~inside
            if edge.any():
                bi = (np.arange(mi)[:, None] + 1 + di) * np.ones((1, mj), dtype=np.int64)
                bj = (np.arange(mj)[None, :] + 1 + dj) * np.ones((mi, 1), dtype=np.int64)
                rhs[edge] -= coef[edge][:, None] * r[bi[edge], bj[edge]]
        matrix = sp.csc_matrix(
            (np.concatenate(vals_), (np.concatenate(rows_), np.concatenate(cols_))),
            shape=(mi * mj, mi * mj),
        )
        solution = spla.splu(matrix).solve(rhs.reshape(-1, 2))
        new = solution.reshape(mi, mj, 2)
        change = float(np.max(np.abs(new - ru)))
        r[1:-1, 1:-1] = new
        if change <= params.tolerance * extent:
            break
    return r
