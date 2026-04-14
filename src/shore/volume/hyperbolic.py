"""Hyperbolic marching: extrude a structured surface layer into a volume.

Algorithm
---------
Given the k=0 surface layer (ni, nj, 3), each subsequent layer k+1 is
computed by:

    P[i,j,k+1] = P[i,j,k] + ds * r^k * n̂[i,j,k]

where:
- ds   : first-layer thickness
- r    : geometric growth ratio
- n̂    : outward unit normal at node (i,j,k), estimated from the four
         in-plane neighbours in the structured grid

After each layer, optional Laplacian smoothing is applied to suppress
high-frequency wiggles.  A per-cell Jacobian check detects inverted cells.

The j-direction is treated as periodic (O-grid topology).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from shore.volume.jacobian import min_jacobian_two_layers


# ---------------------------------------------------------------------------
# Normal estimation
# ---------------------------------------------------------------------------

def _compute_normals(layer: np.ndarray) -> np.ndarray:
    """Estimate outward unit normals for every node of a structured layer.

    Central differences in i (one-sided at boundaries) and j (periodic).

    Parameters
    ----------
    layer :
        Node coordinates, shape (ni, nj, 3).

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
    """
    nj = layer.shape[1]

    ti = np.empty_like(layer)
    ti[1:-1] = layer[2:] - layer[:-2]
    ti[0]    = layer[1]  - layer[0]
    ti[-1]   = layer[-1] - layer[-2]

    jp1 = (np.arange(nj) + 1) % nj
    jm1 = (np.arange(nj) - 1) % nj
    tj  = layer[:, jp1] - layer[:, jm1]

    normals = np.cross(ti, tj)
    mag = np.linalg.norm(normals, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    return normals / mag


# ---------------------------------------------------------------------------
# Laplacian smoothing
# ---------------------------------------------------------------------------

def _laplacian_smooth(layer: np.ndarray, strength: float) -> np.ndarray:
    """One step of in-plane Laplacian smoothing. Periodic in j, clamped at i.

    The k=0 wall row (index 0) is never moved.

    Parameters
    ----------
    layer :
        Shape (ni, nj, 3).
    strength :
        Blend factor in [0, 1].

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
    """
    nj = layer.shape[1]
    jp1 = (np.arange(nj) + 1) % nj
    jm1 = (np.arange(nj) - 1) % nj

    laplacian = (
        layer[:, jp1]
        + layer[:, jm1]
        + np.roll(layer, -1, axis=0)
        + np.roll(layer,  1, axis=0)
    ) / 4.0

    smoothed      = layer + strength * (laplacian - layer)
    smoothed[0]   = layer[0]
    smoothed[-1]  = layer[-1]
    return smoothed


# ---------------------------------------------------------------------------
# Main marching routine
# ---------------------------------------------------------------------------

def march(
    surface: np.ndarray,
    nk: int,
    ds: float,
    growth: float = 1.15,
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    layer_callback: Callable[[int, float, float], None] | None = None,
) -> np.ndarray:
    """Extrude a structured surface layer into a volume by hyperbolic marching.

    Parameters
    ----------
    surface :
        k=0 wall layer, shape (ni, nj, 3).
    nk :
        Total number of layers in the wall-normal direction (including k=0).
    ds :
        First-layer thickness (physical length units matching the STL).
    growth :
        Geometric growth ratio between successive layers.
    smooth_strength :
        Laplacian blend factor applied after each layer (0 = off).
    smooth_iters :
        Number of Laplacian sweeps per layer.
    layer_callback :
        Optional callable invoked after each accepted layer with signature
        ``(k, jmin, step)`` where k is 1-based, jmin is the minimum cell
        Jacobian, and step is the physical layer thickness.

    Returns
    -------
    numpy.ndarray, shape (ni, nj, nk, 3)
        Full structured volume grid.

    Raises
    ------
    RuntimeError
        If a negative-Jacobian cell is detected during marching.
    """
    ni, nj, _ = surface.shape
    grid = np.empty((ni, nj, nk, 3), dtype=np.float64)
    grid[:, :, 0, :] = surface

    for k in range(nk - 1):
        layer     = grid[:, :, k, :]
        normals   = _compute_normals(layer)
        step      = ds * (growth ** k)
        new_layer = layer + step * normals

        for _ in range(smooth_iters):
            new_layer = _laplacian_smooth(new_layer, smooth_strength)

        jmin = min_jacobian_two_layers(layer, new_layer)
        if jmin <= 0.0:
            raise RuntimeError(
                f"Negative Jacobian ({jmin:.3e}) detected at layer k={k+1}.  "
                "Reduce ds, growth ratio, or increase smooth_strength."
            )

        grid[:, :, k + 1, :] = new_layer

        if layer_callback is not None:
            layer_callback(k + 1, jmin, step)

    return grid
