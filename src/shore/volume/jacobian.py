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

        J = (e_i x e_j) . e_k

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
    jp1 = (np.arange(nj) + 1) % nj  # periodic j+1 index

    lk = grid[:-1, :, :-1, :]  # (ni-1, nj, nk-1, 3) — corner
    lk1 = grid[:-1, :, 1:, :]  # k+1 neighbour
    e1 = grid[:-1, jp1, :-1, :] - lk  # j-tangent  (periodic)
    e2 = grid[1:, :, :-1, :] - lk  # i-tangent
    e3 = lk1 - lk  # k-tangent

    return np.einsum("ijkl,ijkl->ijk", np.cross(e2, e1), e3)


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
    lk = grid[:-1, :-1, :-1, :]  # (ni-1, nj-1, nk-1, 3)
    lk1 = grid[:-1, :-1, 1:, :]
    e1 = grid[:-1, 1:, :-1, :] - lk  # j-tangent  (no wrap)
    e2 = grid[1:, :-1, :-1, :] - lk  # i-tangent
    e3 = lk1 - lk  # k-tangent

    return np.einsum("ijkl,ijkl->ijk", np.cross(e2, e1), e3)
