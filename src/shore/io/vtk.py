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
        "Install it with: pip install shore-mesh[vis]  or  pip install pyvista"
    ) from e

from shore.volume.jacobian import _jacobian_field_open, jacobian_field

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
    return np.column_stack(
        [
            grid[:, :, :, 0].ravel(order="F"),
            grid[:, :, :, 1].ravel(order="F"),
            grid[:, :, :, 2].ravel(order="F"),
        ]
    )


def geo_to_pyvista(grid: np.ndarray, with_jacobian: bool = True) -> pv.StructuredGrid:
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
    closed = _close_j(grid)  # (ni, nj+1, nk, 3)
    ni, njc, nk, _ = closed.shape  # njc = nj+1

    sg = pv.StructuredGrid()
    sg.dimensions = (ni, njc, nk)
    sg.points = _grid_points(closed)

    if with_jacobian:
        # Jacobian on the original (periodic) grid: shape (ni-1, nj, nk-1)
        jac = jacobian_field(grid).ravel(order="F")  # (ni-1)*nj*(nk-1)
        sg.cell_data["Jacobian"] = jac

    return sg


def surface_to_pyvista(k0: np.ndarray) -> pv.StructuredGrid:
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


def open_to_pyvista(grid: np.ndarray, with_jacobian: bool = True) -> pv.StructuredGrid:
    """Convert a non-periodic (open) structured block to a PyVista StructuredGrid.

    Unlike ``geo_to_pyvista``, no j-direction closing is performed.  Used for
    polar cap blocks and the equatorial sub-blocks of the cubed-sphere topology.

    Parameters
    ----------
    grid :
        Volume grid, shape (ni, nj, nk, 3).
    with_jacobian :
        If True (default), attach the per-cell Jacobian field as
        ``CellData["Jacobian"]`` using the open (non-periodic) Jacobian.

    Returns
    -------
    pyvista.StructuredGrid
    """
    ni, nj, nk, _ = grid.shape
    sg = pv.StructuredGrid()
    sg.dimensions = (ni, nj, nk)
    sg.points = _grid_points(grid)

    if with_jacobian:
        jac = _jacobian_field_open(grid).ravel(order="F")
        sg.cell_data["Jacobian"] = jac

    return sg


def write_caps_vtm(
    path: str | Path,
    blocks: list[np.ndarray],
    with_jacobian: bool = True,
    labels: list[str] | None = None,
    boxes: list | None = None,
) -> None:
    """Export a list of structured grids as a PyVista MultiBlock (.vtm).

    Originally written for the 6-block cubed-sphere topology — when *labels*
    is None and *blocks* has exactly 6 entries, the cubed-sphere block labels
    (``sub0..sub3``, ``cap_north``, ``cap_south``) are used.  Otherwise
    *labels* must be provided with one entry per block.

    Parameters
    ----------
    path :
        Output path. The ``.vtm`` extension is used if no extension given.
    blocks :
        List of volume grids, each shape ``(ni, nj, nk, 3)``.
    with_jacobian :
        Attach per-cell Jacobian as CellData on each block.
    labels :
        Optional list of block labels (used as the keys of the resulting
        ``MultiBlock``).  Must match ``len(blocks)`` when given.
    boxes :
        Optional list of :class:`shore.io.boxes.Box` instances (the
        ``cc.par`` marker boxes for this topology).  When provided,
        each box is appended as an additional MultiBlock child with
        label ``"box_<n>"`` (1-based) so ParaView shows the meshes
        and the marker boxes together.  See :func:`write_boxes_vtm`
        for the standalone form.
    """
    path = Path(path)
    if path.suffix not in {".vtm", ".vtmb"}:
        path = path.with_suffix(".vtm")
    if labels is None:
        if len(blocks) == 6:
            labels = ["sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south"]
        else:
            labels = [f"block_{i}" for i in range(len(blocks))]
    if len(labels) != len(blocks):
        raise ValueError(f"write_caps_vtm: got {len(blocks)} blocks but {len(labels)} labels")
    mb = pv.MultiBlock()
    for label, block in zip(labels, blocks, strict=True):
        mb[label] = open_to_pyvista(block, with_jacobian=with_jacobian)
    if boxes:
        for q, box in enumerate(boxes, start=1):
            mb[f"box_{q}"] = _box_to_unstructured(box, box_index=q)
    mb.save(str(path))


# ---------------------------------------------------------------------------
# Marker boxes (cc.par boxes section) → VTK
# ---------------------------------------------------------------------------

# SHORE canonical (i, j, k) bit-layout vertex order → VTK hexahedron
# winding (bottom quad CCW, then top quad CCW, looking from +k).  See
# :class:`shore.io.boxes.Box` for the SHORE convention.
_SHORE_TO_VTK_HEX_ORDER: tuple[int, ...] = (0, 1, 3, 2, 4, 5, 7, 6)


def _box_to_unstructured(box, box_index: int) -> pv.UnstructuredGrid:
    """Build a 1-cell ``pv.UnstructuredGrid`` (VTK_HEXAHEDRON) from a Box.

    The single cell carries ``associated_block`` (string) and
    ``box_index`` (int, 1-based) as cell data so ParaView users can
    colour by block or pick a specific box.
    """
    # Reorder vertices SHORE → VTK and emit a single 8-node connectivity.
    points = np.asarray(box.vertices, dtype=np.float64)[list(_SHORE_TO_VTK_HEX_ORDER)]
    # PyVista cells: [n_pts_in_cell, p0, p1, ..., p_{n-1}].
    cells = np.array([8, 0, 1, 2, 3, 4, 5, 6, 7], dtype=np.int64)
    cell_types = np.array([pv.CellType.HEXAHEDRON], dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells, cell_types, points)
    grid.cell_data["associated_block"] = np.array([box.associated_block], dtype=object)
    grid.cell_data["box_index"] = np.array([box_index], dtype=np.int64)
    return grid


def write_boxes_vtm(
    boxes: list,
    path: str | Path,
) -> Path:
    """Export ``cc.par`` marker boxes as a PyVista MultiBlock (.vtm).

    Each :class:`shore.io.boxes.Box` becomes one ``vtkUnstructuredGrid``
    with a single ``VTK_HEXAHEDRON`` cell, exposed under the
    MultiBlock key ``"box_<n>"`` (1-based).  Cell data fields:
    ``associated_block`` and ``box_index``.

    Parameters
    ----------
    boxes :
        Non-empty list of :class:`Box` instances.
    path :
        Output path; ``.vtm`` is used when no extension is given.

    Returns
    -------
    Path
        The path actually written.

    Raises
    ------
    ValueError
        On an empty ``boxes`` list.
    """
    if not boxes:
        raise ValueError("write_boxes_vtm: boxes list must not be empty")
    out = Path(path)
    if out.suffix not in {".vtm", ".vtmb"}:
        out = out.with_suffix(".vtm")
    mb = pv.MultiBlock()
    for q, box in enumerate(boxes, start=1):
        mb[f"box_{q}"] = _box_to_unstructured(box, box_index=q)
    mb.save(str(out))
    return out
