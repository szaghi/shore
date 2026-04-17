"""Jacobian field computation for structured volume grids.

All public functions operate on grids of shape (ni, nj, nk, 3).
The j-direction is treated as periodic (O-grid topology).
"""

from __future__ import annotations

import numpy as np


def jacobian_field(grid: np.ndarray) -> np.ndarray:
    """Compute the per-cell Jacobian scalar field for a structured volume.

    Each cell Jacobian is the scalar triple product of the three edge
    vectors from the cell's (i, j, k) corner:

        J = (e_i × e_j) · e_k

    where e_i, e_j, e_k are the edge vectors in the i-, j-, and k-directions
    respectively. J > 0 for a valid cell; J ≤ 0 for an inverted cell.

    The j-direction is periodic: the cell at j = nj-1 uses j = 0 as its
    j+1 neighbour. The i-direction is not periodic: cells only exist for
    i ∈ [0, ni-2].

    Parameters
    ----------
    grid :
        Structured volume, shape (ni, nj, nk, 3).

    Returns
    -------
    numpy.ndarray, shape (ni-1, nj, nk-1)
        Per-cell Jacobian values. There are nj cells in the j-direction
        (one per column, with the last column wrapping to j=0).
    """
    nj = grid.shape[1]
    jp1 = (np.arange(nj) + 1) % nj   # periodic j+1 index

    lk  = grid[:-1, :,   :-1, :]          # (ni-1, nj, nk-1, 3) — corner
    lk1 = grid[:-1, :,    1:, :]          # k+1 neighbour
    e1  = grid[:-1, jp1, :-1, :] - lk    # j-tangent  (periodic)
    e2  = grid[1:,  :,   :-1, :] - lk    # i-tangent
    e3  = lk1 - lk                        # k-tangent

    return np.einsum("ijkl,ijkl->ijk", np.cross(e2, e1), e3)


def jacobian_per_layer(grid: np.ndarray) -> np.ndarray:
    """Return the minimum cell Jacobian for each wall-normal layer interval.

    Parameters
    ----------
    grid :
        Structured volume, shape (ni, nj, nk, 3).

    Returns
    -------
    numpy.ndarray, shape (nk-1,)
        ``result[k]`` is the minimum Jacobian among all cells between
        layer k and layer k+1.
    """
    return jacobian_field(grid).min(axis=(0, 1))


def min_jacobian_two_layers(
    layer_k: np.ndarray,
    layer_k1: np.ndarray,
    periodic_j: bool = True,
) -> float:
    """Return the minimum Jacobian between two consecutive structured layers.

    Convenience wrapper used during hyperbolic marching.

    Parameters
    ----------
    layer_k, layer_k1 :
        Consecutive layers, shape (ni, nj, 3).
    periodic_j :
        When ``True`` (default) the j-direction wraps around (O-grid).
        When ``False`` the last j-column does not wrap; the cell count in j
        is ``nj-1`` instead of ``nj``.

    Returns
    -------
    float
        Minimum cell Jacobian. Negative → inverted cell.
    """
    two_layer = np.stack([layer_k, layer_k1], axis=2)  # (ni, nj, 2, 3)
    if periodic_j:
        return float(jacobian_field(two_layer).min())
    return float(_jacobian_field_open(two_layer).min())


def _jacobian_field_open(grid: np.ndarray) -> np.ndarray:
    """Jacobian field for a non-periodic (open) structured volume.

    Identical to ``jacobian_field`` but the j-direction is not wrapped:
    there are ``nj-1`` cells in j rather than ``nj``.

    Parameters
    ----------
    grid :
        Structured volume, shape (ni, nj, nk, 3).

    Returns
    -------
    numpy.ndarray, shape (ni-1, nj-1, nk-1)
    """
    lk  = grid[:-1, :-1, :-1, :]           # (ni-1, nj-1, nk-1, 3)
    lk1 = grid[:-1, :-1,  1:, :]
    e1  = grid[:-1, 1:,  :-1, :] - lk      # j-tangent  (no wrap)
    e2  = grid[1:,  :-1, :-1, :] - lk      # i-tangent
    e3  = lk1 - lk                          # k-tangent

    return np.einsum("ijkl,ijkl->ijk", np.cross(e2, e1), e3)
