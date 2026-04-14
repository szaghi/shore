"""VTK export for ParaView.

Writes structured grids to VTK XML StructuredGrid format (.vts), which
ParaView reads natively.

O-grid periodicity
------------------
VTK StructuredGrid has no concept of periodic connectivity. The standard
approach for O-grids is to close the j-direction by repeating the first
j-column at the end, producing a grid of shape (ni, nj+1, nk). This
creates a visible seam in ParaView that is geometrically correct.

Jacobian cell data
------------------
The full per-cell Jacobian field is attached as ``CellData`` under the
array name ``"Jacobian"``, enabling direct quality visualisation in
ParaView without any post-processing filter.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import pyvista as pv
except ImportError as e:
    raise ImportError(
        "pyvista is required for VTK export.  "
        "Install it with: pip install shore[vis]  or  pip install pyvista"
    ) from e

from shore.volume.jacobian import jacobian_field


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def _close_j(grid: np.ndarray) -> np.ndarray:
    """Close the periodic j-direction by appending j=0 at the end.

    Parameters
    ----------
    grid :
        Shape (ni, nj, nk, 3).

    Returns
    -------
    numpy.ndarray, shape (ni, nj+1, nk, 3)
    """
    return np.concatenate([grid, grid[:, 0:1, :, :]], axis=1)


def _grid_points(grid: np.ndarray) -> np.ndarray:
    """Flatten a structured grid to (N, 3) in PyVista point order (i-fastest).

    Parameters
    ----------
    grid :
        Shape (ni, nj, nk, 3).

    Returns
    -------
    numpy.ndarray, shape (ni*nj*nk, 3)
    """
    return np.column_stack([
        grid[:, :, :, 0].ravel(order="F"),
        grid[:, :, :, 1].ravel(order="F"),
        grid[:, :, :, 2].ravel(order="F"),
    ])


def geo_to_pyvista(grid: np.ndarray, with_jacobian: bool = True) -> "pv.StructuredGrid":
    """Convert a structured volume grid to a PyVista StructuredGrid.

    The j-direction is closed (first column repeated at end) so that
    the O-grid forms a topologically complete shell.

    Parameters
    ----------
    grid :
        Volume grid, shape (ni, nj, nk, 3).
    with_jacobian :
        If True (default), attach the per-cell Jacobian field as
        ``CellData["Jacobian"]``.

    Returns
    -------
    pyvista.StructuredGrid
    """
    closed = _close_j(grid)              # (ni, nj+1, nk, 3)
    ni, njc, nk, _ = closed.shape       # njc = nj+1

    sg = pv.StructuredGrid()
    sg.dimensions = (ni, njc, nk)
    sg.points = _grid_points(closed)

    if with_jacobian:
        # Jacobian on the original (periodic) grid: shape (ni-1, nj, nk-1)
        jac = jacobian_field(grid).ravel(order="F")  # (ni-1)*nj*(nk-1)
        sg.cell_data["Jacobian"] = jac

    return sg


def surface_to_pyvista(k0: np.ndarray) -> "pv.StructuredGrid":
    """Convert a k=0 surface layer to a PyVista StructuredGrid (nk=1 slab).

    Parameters
    ----------
    k0 :
        Surface layer, shape (ni, nj, 3).

    Returns
    -------
    pyvista.StructuredGrid
    """
    ni, nj, _ = k0.shape
    closed = np.concatenate([k0, k0[:, 0:1, :]], axis=1)  # (ni, nj+1, 3)
    nijc = nj + 1

    sg = pv.StructuredGrid()
    sg.dimensions = (ni, nijc, 1)
    sg.points = _grid_points(closed[:, :, np.newaxis, :])
    return sg


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_vts(path: str | Path, grid: np.ndarray, with_jacobian: bool = True) -> None:
    """Export a volume grid to VTK XML StructuredGrid format (.vts).

    Parameters
    ----------
    path :
        Output file path. The ``.vts`` extension is appended if absent.
    grid :
        Volume grid, shape (ni, nj, nk, 3).
    with_jacobian :
        Attach per-cell Jacobian as CellData (default True).
    """
    path = Path(path)
    if path.suffix != ".vts":
        path = path.with_suffix(".vts")
    sg = geo_to_pyvista(grid, with_jacobian=with_jacobian)
    sg.save(str(path))


def write_surface_vts(path: str | Path, k0: np.ndarray) -> None:
    """Export a k=0 surface layer to VTK XML StructuredGrid format (.vts).

    Parameters
    ----------
    path :
        Output file path. The ``.vts`` extension is appended if absent.
    k0 :
        Surface layer, shape (ni, nj, 3).
    """
    path = Path(path)
    if path.suffix != ".vts":
        path = path.with_suffix(".vts")
    sg = surface_to_pyvista(k0)
    sg.save(str(path))
