"""File-based inspection utilities: info, check, export, view.

These functions operate on already-written ``.geo`` / ``.npz`` / ``.vtm``
files without requiring a live ``Mesh`` object in memory.  They mirror the
in-memory ``Mesh`` methods for offline or post-hoc inspection.
"""

from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

def info_geo(path: str | Path) -> dict:
    """Return metadata for a single-block ``.geo`` file.

    Parameters
    ----------
    path :
        Path to a ``.geo`` file.

    Returns
    -------
    dict with keys ``shape``, ``n_points``, ``n_cells``,
    ``xyz_min``, ``xyz_max``.
    """
    from shore.io.geo import read_geo
    grid = read_geo(path)
    ni, nj, nk, _ = grid.shape
    pts = grid.reshape(-1, 3)
    return {
        "shape": (ni, nj, nk),
        "n_points": ni * nj * nk,
        "n_cells": (ni - 1) * (nj - 1) * (nk - 1),
        "xyz_min": pts.min(axis=0).tolist(),
        "xyz_max": pts.max(axis=0).tolist(),
    }


def info_surface(path: str | Path) -> dict:
    """Return metadata for a ``.npz`` surface file.

    Parameters
    ----------
    path :
        Path to a ``.npz`` file with a ``'k0'`` key.

    Returns
    -------
    dict with keys ``shape``, ``n_points``, ``xyz_min``, ``xyz_max``.
    """
    from shore.mesh.body import Body
    k0 = Body.load_surface(path)
    ni, nj, _ = k0.shape
    pts = k0.reshape(-1, 3)
    return {
        "shape": (ni, nj),
        "n_points": ni * nj,
        "xyz_min": pts.min(axis=0).tolist(),
        "xyz_max": pts.max(axis=0).tolist(),
    }


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------

def check_geo(
    path: str | Path,
    warn_threshold: float = 1e-6,
    periodic_j: bool = False,
) -> dict:
    """Check cell Jacobians layer-by-layer for a single-block ``.geo`` file.

    Parameters
    ----------
    path :
        Path to a ``.geo`` file.
    warn_threshold :
        Warn if jmin drops below this (but remains positive).
    periodic_j :
        Use the periodic-j Jacobian (``True``) or the open (non-periodic)
        Jacobian (``False``, default).  Caps-topology blocks are always
        open; use ``True`` only for legacy single-block O-grids.

    Returns
    -------
    dict with keys ``ok`` (bool), ``jmins`` (list of per-layer minimums),
    ``any_negative``, ``any_warn``.
    """
    from shore.io.geo import read_geo
    from shore.volume.jacobian import _jacobian_field_open, jacobian_field
    grid  = read_geo(path)
    nk    = grid.shape[2]
    jfunc = jacobian_field if periodic_j else _jacobian_field_open
    jac   = jfunc(grid)   # (ni-1, nj[-1], nk-1)
    jmins = [float(jac[:, :, k].min()) for k in range(nk - 1)]
    any_neg  = any(j <= 0.0 for j in jmins)
    any_warn = any(0.0 < j < warn_threshold for j in jmins)
    return {
        "ok": not any_neg,
        "jmins": jmins,
        "any_negative": any_neg,
        "any_warn": any_warn,
    }


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

def export_geo(
    path: str | Path,
    output: str | Path | None = None,
    with_jacobian: bool = True,
) -> Path:
    """Export a ``.geo`` file to VTK StructuredGrid ``.vts`` for ParaView.

    Parameters
    ----------
    path :
        Input ``.geo`` file.
    output :
        Output ``.vts`` path.  Defaults to same stem with ``.vts`` extension.
    with_jacobian :
        Attach per-cell Jacobian as CellData.

    Returns
    -------
    Path
        The written ``.vts`` file path.
    """
    from shore.io.geo import read_geo
    from shore.io.vtk import write_vts
    path   = Path(path)
    out    = Path(output) if output else path.with_suffix(".vts")
    grid   = read_geo(path)
    write_vts(out, grid, with_jacobian=with_jacobian)
    return out


def export_surface(
    path: str | Path,
    output: str | Path | None = None,
) -> Path:
    """Export a ``.npz`` surface file to VTK StructuredGrid ``.vts``.

    Parameters
    ----------
    path :
        Input ``.npz`` file with ``'k0'`` key.
    output :
        Output ``.vts`` path.  Defaults to same stem with ``.vts`` extension.

    Returns
    -------
    Path
        The written ``.vts`` file path.
    """
    from shore.io.vtk import write_surface_vts
    from shore.mesh.body import Body
    path = Path(path)
    out  = Path(output) if output else path.with_suffix(".vts")
    k0   = Body.load_surface(path)
    write_surface_vts(out, k0)
    return out


# ---------------------------------------------------------------------------
# view
# ---------------------------------------------------------------------------

def view_geo_file(
    path: str | Path,
    *,
    colorby: str = "Jacobian",
    show_edges: bool = True,
    opacity: float = 0.6,
    jupyter: bool = False,
    off_screen: bool = False,
) -> object:
    """Open an interactive PyVista window for a ``.geo`` file.

    Parameters
    ----------
    path :
        Input ``.geo`` file.
    colorby :
        ``"Jacobian"`` (default) or ``"none"``.
    show_edges :
        Draw cell edges.
    opacity :
        Grid opacity [0-1].
    jupyter :
        Use HTML backend for inline Jupyter display.
    off_screen :
        Render off-screen (for testing).

    Returns
    -------
    pyvista.Plotter
    """
    from shore.io.geo import read_geo
    from shore.vis import view_geo
    grid = read_geo(path)
    return view_geo(
        grid,
        jupyter=jupyter,
        off_screen=off_screen,
        show_edges=show_edges,
        opacity=opacity,
        colorby=colorby,
    )


def view_surface_file(
    path: str | Path,
    *,
    colorby: str = "radius",
    show_edges: bool = True,
    jupyter: bool = False,
    off_screen: bool = False,
) -> object:
    """Open an interactive PyVista window for a ``.npz`` surface file.

    Parameters
    ----------
    path :
        Input ``.npz`` file with ``'k0'`` key.
    colorby :
        ``"radius"`` (default) or ``"none"``.
    show_edges :
        Draw cell edges.
    jupyter :
        Use HTML backend for inline Jupyter display.
    off_screen :
        Render off-screen (for testing).

    Returns
    -------
    pyvista.Plotter
    """
    from shore.mesh.body import Body
    from shore.vis import view_surface
    k0 = Body.load_surface(path)
    return view_surface(
        k0,
        jupyter=jupyter,
        off_screen=off_screen,
        show_edges=show_edges,
        colorby=colorby,
    )
