"""Synthetic multi-block assembly for the GUI spike (ADR 0006).

A spherical shell split into ``n_phi * n_theta`` curvilinear blocks, each a
``(ni, nj, nk, 3)`` float64 node array. The k direction is wall-normal with a
geometric law whose first step is ``ds_wall``, so the float32 precision trap
(small wall spacing at a large world offset) can be exercised.

Throwaway spike code: not part of the ``shore`` package.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Block:
    """One structured block of the synthetic assembly."""

    name: str
    nodes: np.ndarray  # (ni, nj, nk, 3) float64

    @property
    def n_nodes(self) -> int:
        """Number of nodes in the block."""
        ni, nj, nk, _ = self.nodes.shape
        return ni * nj * nk


def _geometric_steps(first: float, total: float, n: int) -> np.ndarray:
    """``n`` geometric steps starting at ``first`` and summing to ``total``.

    The growth ratio is solved by bisection; if ``first * n >= total`` the
    steps are uniform.
    """
    if first * n >= total:
        return np.full(n, total / n)
    lo, hi = 1.0 + 1e-12, 10.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        series = first * (r**n - 1.0) / (r - 1.0)
        lo, hi = (r, hi) if series < total else (lo, r)
    steps = first * r ** np.arange(n)
    return steps * (total / steps.sum())


def _geometric_radii(radius: float, thickness: float, ds_wall: float, nk: int) -> np.ndarray:
    """Radii of ``nk`` layers: first step ``ds_wall``, total ``thickness``."""
    steps = _geometric_steps(ds_wall, thickness, nk - 1)
    return radius + np.concatenate(([0.0], np.cumsum(steps)))


def _clustered(t0: float, t1: float, n: int, width: float) -> np.ndarray:
    """``n`` nodes on ``[t0, t1]`` clustered about an interior centre node.

    The two cells adjacent to the centre node have width ``width``; cells grow
    geometrically towards both ends.
    """
    left = (n - 1) // 2
    right = n - 1 - left
    tc = t0 + (t1 - t0) * left / (n - 1)
    sl = _geometric_steps(width, tc - t0, left)
    sr = _geometric_steps(width, t1 - tc, right)
    return np.concatenate((tc - np.cumsum(sl)[::-1], [tc], tc + np.cumsum(sr)))


def make_assembly(
    n_nodes: int = 10_000_000,
    n_phi: int = 4,
    n_theta: int = 5,
    radius: float = 1.0,
    thickness: float = 2.0,
    ds_wall: float = 1e-6,
    offset: tuple[float, float, float] = (1.0e3, 0.0, 0.0),
    aspect_ratio: float | None = None,
) -> list[Block]:
    """Build a shell assembly with about ``n_nodes`` nodes in total.

    Parameters
    ----------
    n_nodes : int
        Target total node count; each block gets ``n_nodes / n_blocks`` nodes
        as a near-cubic ``ni = nj = nk`` (k slightly favoured).
    n_phi, n_theta : int
        Blocks in azimuth and polar angle. Polar angle spans
        ``[0.1 pi, 0.9 pi]`` so no block touches a pole (no singular rows).
    radius, thickness : float
        Wall radius and shell thickness.
    ds_wall : float
        First wall-normal step.
    offset : tuple of float
        World-space translation of the whole assembly.
    aspect_ratio : float, optional
        If given, polar angle is clustered about each block's centre node so
        that the two wall cells there have tangential width
        ``aspect_ratio * ds_wall`` (realistic boundary-layer aspect ratio,
        independent of the node count). Otherwise polar angle is uniform.

    Returns
    -------
    list of Block
    """
    n_blocks = n_phi * n_theta
    side = max(round((n_nodes / n_blocks) ** (1.0 / 3.0)), 3)
    ni = nj = nk = side
    radii = _geometric_radii(radius, thickness, ds_wall, nk)
    origin = np.asarray(offset, dtype=np.float64)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, n_phi + 1)
    theta_edges = np.linspace(0.1 * np.pi, 0.9 * np.pi, n_theta + 1)
    blocks: list[Block] = []
    for a in range(n_phi):
        phi = np.linspace(phi_edges[a], phi_edges[a + 1], ni)
        for b in range(n_theta):
            t0, t1 = theta_edges[b], theta_edges[b + 1]
            if aspect_ratio is None:
                theta = np.linspace(t0, t1, nj)
            else:
                theta = _clustered(t0, t1, nj, aspect_ratio * ds_wall / radius)
            ph, th, rr = np.meshgrid(phi, theta, radii, indexing="ij")
            nodes = np.empty((ni, nj, nk, 3), dtype=np.float64)
            nodes[..., 0] = rr * np.sin(th) * np.cos(ph)
            nodes[..., 1] = rr * np.sin(th) * np.sin(ph)
            nodes[..., 2] = rr * np.cos(th)
            nodes += origin
            blocks.append(Block(name=f"blk_{a}_{b}", nodes=nodes))
    return blocks
