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

from shore.exceptions import MarchingError
from shore.volume.jacobian import min_jacobian_two_layers


# ---------------------------------------------------------------------------
# Normal estimation
# ---------------------------------------------------------------------------

def _compute_normals(
    layer: np.ndarray,
    periodic: tuple[bool, bool] = (False, True),
) -> np.ndarray:
    """Estimate outward unit normals for every node of a structured layer.

    Central differences in both directions; one-sided at non-periodic boundaries.

    Parameters
    ----------
    layer :
        Node coordinates, shape (ni, nj, 3).
    periodic :
        ``(periodic_i, periodic_j)``.  When ``True`` the corresponding
        direction uses modular (wrap-around) indexing for the end nodes;
        when ``False`` one-sided differences are used at the boundaries.
        Default ``(False, True)`` preserves the original O-grid behaviour
        (i open, j periodic).

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
    """
    ni, nj = layer.shape[0], layer.shape[1]
    periodic_i, periodic_j = periodic

    # --- i-direction tangent ---
    ti = np.empty_like(layer)
    if periodic_i:
        ip1 = (np.arange(ni) + 1) % ni
        im1 = (np.arange(ni) - 1) % ni
        ti  = layer[ip1] - layer[im1]
    else:
        ti[1:-1] = layer[2:] - layer[:-2]
        ti[0]    = layer[1]  - layer[0]
        ti[-1]   = layer[-1] - layer[-2]

    # --- j-direction tangent ---
    if periodic_j:
        jp1 = (np.arange(nj) + 1) % nj
        jm1 = (np.arange(nj) - 1) % nj
        tj  = layer[:, jp1] - layer[:, jm1]
    else:
        tj = np.empty_like(layer)
        tj[:, 1:-1] = layer[:, 2:] - layer[:, :-2]
        tj[:, 0]    = layer[:, 1]  - layer[:, 0]
        tj[:, -1]   = layer[:, -1] - layer[:, -2]

    normals = np.cross(ti, tj)
    mag = np.linalg.norm(normals, axis=-1, keepdims=True)
    mag = np.where(mag < 1e-15, 1.0, mag)
    return normals / mag


# ---------------------------------------------------------------------------
# Laplacian smoothing
# ---------------------------------------------------------------------------

def _laplacian_smooth(
    layer: np.ndarray,
    strength: float,
    periodic: tuple[bool, bool] = (False, True),
) -> np.ndarray:
    """One step of in-plane Laplacian smoothing.

    Boundary rows/columns in non-periodic directions are clamped (never moved).
    The k=0 wall row (index 0) is always clamped regardless of periodicity.

    Parameters
    ----------
    layer :
        Shape (ni, nj, 3).
    strength :
        Blend factor in [0, 1].
    periodic :
        ``(periodic_i, periodic_j)``.  Non-periodic boundaries are clamped.
        Default ``(False, True)`` preserves the original O-grid behaviour.

    Returns
    -------
    numpy.ndarray, shape (ni, nj, 3)
    """
    ni, nj = layer.shape[0], layer.shape[1]
    periodic_i, periodic_j = periodic

    # j neighbours
    if periodic_j:
        jp1 = (np.arange(nj) + 1) % nj
        jm1 = (np.arange(nj) - 1) % nj
        j_fwd = layer[:, jp1]
        j_bwd = layer[:, jm1]
    else:
        j_fwd = np.empty_like(layer)
        j_bwd = np.empty_like(layer)
        j_fwd[:, :-1] = layer[:, 1:]
        j_fwd[:, -1]  = layer[:, -1]   # clamp: boundary contributes itself
        j_bwd[:, 1:]  = layer[:, :-1]
        j_bwd[:, 0]   = layer[:, 0]

    # i neighbours
    if periodic_i:
        i_fwd = np.roll(layer, -1, axis=0)
        i_bwd = np.roll(layer,  1, axis=0)
    else:
        i_fwd = np.empty_like(layer)
        i_bwd = np.empty_like(layer)
        i_fwd[:-1] = layer[1:]
        i_fwd[-1]  = layer[-1]
        i_bwd[1:]  = layer[:-1]
        i_bwd[0]   = layer[0]

    laplacian = (j_fwd + j_bwd + i_fwd + i_bwd) / 4.0
    smoothed  = layer + strength * (laplacian - layer)

    # always clamp the wall row
    smoothed[0] = layer[0]
    # clamp non-periodic boundaries
    if not periodic_i:
        smoothed[-1] = layer[-1]
    if not periodic_j:
        smoothed[:, 0]  = layer[:, 0]
        smoothed[:, -1] = layer[:, -1]

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
    periodic: tuple[bool, bool] = (False, True),
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
    periodic :
        ``(periodic_i, periodic_j)`` flags passed to normal estimation and
        Laplacian smoothing.  Default ``(False, True)`` is the standard
        O-grid topology (i open, j periodic).  Use ``(False, False)`` for
        polar cap blocks.
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
        normals   = _compute_normals(layer, periodic=periodic)
        step      = ds * (growth ** k)
        new_layer = layer + step * normals

        for _ in range(smooth_iters):
            new_layer = _laplacian_smooth(new_layer, smooth_strength, periodic=periodic)

        jmin = min_jacobian_two_layers(layer, new_layer, periodic_j=periodic[1])
        if jmin <= 0.0:
            raise MarchingError(
                f"Negative Jacobian ({jmin:.3e}) detected at layer k={k+1}.  "
                "Reduce ds, growth ratio, or increase smooth_strength."
            )

        grid[:, :, k + 1, :] = new_layer

        if layer_callback is not None:
            layer_callback(k + 1, jmin, step)

    return grid
