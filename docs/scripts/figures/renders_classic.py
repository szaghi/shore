"""3-D renders of the classic SHORE grids: cubed sphere, O-/C-grids, hybrids, primitives, Chimera.

Every grid is built with SHORE's own API (``Mesh``, the topologies, ``shore.primitive``,
``shore.boxes``) at a small size, cached per run so the light and dark variants share it,
and rendered off-screen with PyVista through ``theme.plotter`` / ``theme.shot``.
Temporary STL files go to a temporary directory, never into the repository.
"""

from __future__ import annotations

import tempfile
from functools import cache
from pathlib import Path

import numpy as np
import pyvista as pv
from theme import PALETTES, plotter, shot

# ---------------------------------------------------------------------------
# Builders (cached: light and dark reuse one grid)
# ---------------------------------------------------------------------------


@cache
def _tmp() -> Path:
    return Path(tempfile.mkdtemp(prefix="shore-figs-"))


@cache
def _sphere_stl(radius: float = 1.0) -> Path:
    from shore.shapes import sphere

    path = _tmp() / f"sphere-{radius:g}.stl"
    sphere(radius=radius, subdivisions=4).export(str(path))
    return path


@cache
def _wing_stl() -> Path:
    from shore.profiles import naca
    from shore.shapes import extrude_section

    path = _tmp() / "naca0012.stl"
    section = naca("0012", n_pts=60, closed_te=True)
    extrude_section(section, span=0.6, n_stations=12, end_caps=True).export(str(path))
    return path


@cache
def _sphere_mesh():
    """Cubed sphere on the unit sphere, marched (6 blocks)."""
    from shore.mesh import CubedSphereTopology, Mesh, Spacing1D

    m = Mesh.from_stl(
        _sphere_stl(),
        topology=CubedSphereTopology(),
        ni=24,
        nj=64,
        nk=22,
        theta_cap_deg=30.0,
        spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
    )
    m.march(ds=0.02, growth=1.1)
    return m


@cache
def _ogrid_mesh():
    """Single-block O-grid around a prolate spheroid (2:1, axis z)."""
    import trimesh

    from shore.mesh import Mesh, OGridTopology, Spacing1D

    path = _tmp() / "spheroid.stl"
    body = trimesh.creation.icosphere(subdivisions=4, radius=1.0)
    body.apply_scale((0.5, 0.5, 1.0))
    body.export(str(path))
    m = Mesh.from_stl(
        path,
        topology=OGridTopology(),
        ni=33,
        nj=48,
        nk=20,
        theta_cap_deg=8.0,
        spacing_k=Spacing1D(law="geometric", ds=0.01, growth=1.12),
    )
    m.march(ds=0.01, growth=1.12)
    return m


@cache
def _cgrid_mesh(cap_tips: str = "none"):
    """Three-block C-grid around the NACA 0012 wing (span along z)."""
    from shore.mesh import Body, CGridTopology, Mesh, Spacing1D

    body = Body.from_stl(_wing_stl(), skip_hygiene=True)
    topo = CGridTopology()
    blocks = topo.build(
        surface=None,
        nk=28,
        mesh=body.mesh,
        spacing_k=Spacing1D(law="geometric", ds=2e-3, growth=1.1),
        ni_body=81,
        ni_wake=20,
        n_stations=9,
        wake_length=2.0,
        wake_growth=1.15,
        span_axis="z",
        blocks=3,
        cap_tips=cap_tips,
    )
    m = Mesh(body=body, blocks=blocks, topology=topo)
    m.march(ds=2e-3, growth=1.1, smooth_strength=0.3, smooth_iters=3)
    return m


@cache
def _ch_mesh():
    """CH hybrid: C-grid wrap + H-channel to a chord/normal frame."""
    from shore.mesh import Body, CHTopology, Mesh, Spacing1D

    body = Body.from_stl(_wing_stl(), skip_hygiene=True)
    topo = CHTopology(
        chord_min=-1.0,
        chord_max=3.0,
        normal_min=-1.2,
        normal_max=1.2,
        nk_channel=10,
        k_law_channel="uniform",
    )
    blocks = topo.build(
        surface=None,
        nk=14,
        mesh=body.mesh,
        spacing_k=Spacing1D(law="geometric", ds=2e-3, growth=1.2),
        ni_body=61,
        ni_wake=16,
        n_stations=7,
        wake_length=2.0,
    )
    m = Mesh(body=body, blocks=blocks, topology=topo)
    m.march(ds=2e-3, growth=1.2, smooth_strength=0.0, smooth_iters=0)
    return m


@cache
def _oh_mesh():
    """OH hybrid: cubed-sphere wrap + H-channel to a box."""
    from shore.mesh import Mesh, OHTopology, Spacing1D

    m = Mesh.from_stl(
        _sphere_stl(),
        topology=OHTopology(
            box_vmin=(-2.5, -2.5, -2.5),
            box_vmax=(2.5, 2.5, 2.5),
            nk_channel=10,
            k_law_channel="uniform",
        ),
        ni=20,
        nj=48,
        nk=14,
        spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
    )
    m.march(ds=0.02, growth=1.1)
    return m


@cache
def _chimera_meshes():
    """Near-body cubed sphere (radius 0.5) and a flat-caps cylinder background along x."""
    from shore.mesh import CubedSphereTopology, Mesh, Spacing1D
    from shore.primitive import FlatCapsCylinder

    near = Mesh.from_stl(
        _sphere_stl(0.5),
        topology=CubedSphereTopology(),
        ni=16,
        nj=40,
        nk=13,
        theta_cap_deg=30.0,
        spacing_k=Spacing1D(law="geometric", ds=0.01, growth=1.15),
    )
    near.march(ds=0.01, growth=1.15)
    background = FlatCapsCylinder.build(
        r_out=1.6, z0=-1.8, z1=3.6, ni=8, nc=13, nk=34, inner_frac=0.45, axis="x"
    )
    return near, background


@cache
def _h_grids():
    from shore.primitive import HGridRectMesh, HGridTubeMesh

    duct = HGridRectMesh.build(
        length=0.6,
        width=0.12,
        height=0.2,
        ni=30,
        nj=16,
        nk=22,
        streamwise_axis="x",
        j_law="tanh2",
        j_beta=4.0,
        k_law="tanh2",
        k_beta=4.0,
    )
    tube = HGridTubeMesh.build(
        length=0.6,
        radius=0.1,
        n_axial=30,
        n_circ=40,
        n_radial=12,
        streamwise_axis="x",
        k_law="tanh",
        k_beta=3.0,
    )
    return duct, tube


@cache
def _primitives():
    from shore.primitive import AnnularCylinderMesh, BoxMesh, FlatCapsCylinder

    box = BoxMesh.build(
        vmin=(0.0, 0.0, 0.0), vmax=(1.2, 0.8, 0.6), ni=20, nj=12, nk=10, i_law="tanh", i_beta=3.0
    )
    annulus = AnnularCylinderMesh.build(
        r_in=0.22,
        r_out=0.6,
        z0=0.0,
        z1=0.6,
        ni=9,
        nj=37,
        nk=10,
        theta_start=0.0,
        theta_end=1.5 * np.pi,
        center=(2.2, 0.4, 0.0),
        i_law="tanh",
        i_beta=2.0,
    )
    caps = FlatCapsCylinder.build(
        r_out=0.6, z0=0.0, z1=0.6, ni=7, nc=11, nk=10, inner_frac=0.5, center=(3.9, 0.4, 0.0)
    )
    return box, annulus, caps


@cache
def _marker_boxes():
    from shore.boxes import boxes_from_stl
    from shore.shapes import capsule

    path = _tmp() / "capsule.stl"
    body = capsule(radius=0.5, height=1.0, sections=48)
    body.apply_transform([[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]])
    body.export(str(path))
    boxes = boxes_from_stl(path, "body", max_outward_gap=0.08, method="voxel-fill")
    return path, boxes


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------


def _grid(nodes: np.ndarray) -> pv.StructuredGrid:
    """A StructuredGrid over a ``(ni, nj, nk, 3)`` or ``(ni, nj, 3)`` node array."""
    a = np.asarray(nodes, dtype=float)
    if a.ndim == 3:
        a = a[:, :, None, :]
    return pv.StructuredGrid(a[..., 0], a[..., 1], a[..., 2])


def _cell_index(nodes: np.ndarray, axis: int) -> np.ndarray:
    """Per-cell index along *axis* in VTK (i-fastest) cell order."""
    shape = tuple(s - 1 for s in nodes.shape[:3])
    return np.indices(shape)[axis].ravel(order="F")


def _keep(grid, mask: np.ndarray):
    """The cells of *grid* where *mask* (one entry per cell) is true."""
    return grid.extract_cells(np.flatnonzero(mask))


def _wedge(grid, centre_deg: float, half_deg: float = 45.0):
    """*grid* less the cells whose centres lie in the azimuthal wedge ``centre +- half``."""
    c = grid.cell_centers().points
    phi = np.degrees(np.arctan2(c[:, 1], c[:, 0]))
    d = (phi - centre_deg + 180.0) % 360.0 - 180.0
    return _keep(grid, np.abs(d) > half_deg)


def _close_j(nodes: np.ndarray) -> np.ndarray:
    """A j-periodic node array with its first j row repeated at the end (for display)."""
    return np.concatenate([nodes, nodes[:, :1]], axis=1)


def _add(pl, mesh, pal, color=None, *, edges=True, opacity=1.0, lw=0.6, smooth=True, **kw):
    """Add *mesh* with the shared material: thin edges, modest specular, ambient 0.3."""
    if mesh.n_cells == 0:
        return None
    return pl.add_mesh(
        mesh,
        color=color,
        show_edges=edges,
        edge_color=pal.edge,
        line_width=lw,
        smooth_shading=smooth,
        specular=0.1,
        specular_power=12.0,
        ambient=0.3,
        diffuse=0.8,
        opacity=opacity,
        show_scalar_bar=False,
        **kw,
    )


def _camera(pl, azim: float = -60.0, elev: float = 25.0, zoom: float = 1.3, parallel=False):
    """Elevated 3/4 view, z up, fitted to the scene and zoomed by *zoom*."""
    a, e = np.radians(azim), np.radians(elev)
    direction = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    bounds = np.array(pl.bounds).reshape(3, 2)
    focus = bounds.mean(axis=1)
    pl.camera_position = [tuple(focus + 10.0 * direction), tuple(focus), (0.0, 0.0, 1.0)]
    if parallel:
        pl.enable_parallel_projection()
    pl.reset_camera()
    pl.camera.zoom(zoom)


def _bar(pl, actor, pal, title: str, fmt: str = "%.2f", n_labels: int = 5) -> None:
    """A small vertical colour bar for *actor* in the ink colour."""
    pl.add_scalar_bar(
        title=title,
        mapper=actor.mapper,
        vertical=True,
        position_x=0.86,
        position_y=0.2,
        width=0.035,
        height=0.5,
        color=pal.ink,
        title_font_size=20,
        label_font_size=17,
        n_labels=n_labels,
        fmt=fmt,
        font_family="arial",
        shadow=False,
    )


def _wing_view(nodes: np.ndarray) -> np.ndarray:
    """Wing coordinates for display: span (z) into -y, thickness (y) up (rotation about x)."""
    a = np.asarray(nodes, dtype=float)
    return np.stack([a[..., 0], -a[..., 2], a[..., 1]], axis=-1)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def sphere_blocks(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    for k, blk in enumerate(_sphere_mesh().blocks):
        _add(pl, _grid(blk.nodes[:, :, 0]), pal, pal.block(k))
    _camera(pl, azim=-60.0, elev=38.0, zoom=1.35)
    shot(pl, "sphere-blocks", mode)


def sphere_topdown(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode, (1000, 1000))
    for k, blk in enumerate(_sphere_mesh().blocks):
        surf = blk.nodes[:, :, 0]
        if surf[..., 2].max() > 0.0:
            _add(pl, _grid(surf), pal, pal.block(k), lw=0.7)
    pl.camera_position = [(0.0, 0.0, 10.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    pl.enable_parallel_projection()
    pl.reset_camera()
    pl.camera.zoom(1.35)
    shot(pl, "sphere-topdown", mode)


def marching_layers(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    m = _sphere_mesh()
    nk = m.blocks[0].nodes.shape[2]
    for blk in m.blocks:
        grid = _grid(blk.nodes)
        grid.cell_data["k"] = _cell_index(blk.nodes, 2)
        part = _wedge(grid, -90.0)
        actor = _add(
            pl, part, pal, scalars="k", cmap=pal.depth, clim=(0, nk - 2), lw=0.5, smooth=False
        )
    body = pv.read(str(_sphere_stl())).scale(0.999)
    _add(pl, body, pal, pal.tint(pal.slate, 0.15), edges=False)
    _camera(pl, azim=-100.0, elev=26.0, zoom=1.35)
    _bar(pl, actor, pal, "layer k", fmt="%.0f", n_labels=3)
    shot(pl, "marching-layers", mode)


def seam_continuity(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    for k, blk in enumerate(_sphere_mesh().blocks):
        _add(pl, _grid(blk.nodes[:, :, 0]), pal, pal.block(k), lw=0.5)
        outer = _wedge(_grid(blk.nodes[:, :, -1]), -90.0)
        _add(pl, outer, pal, pal.block(k), lw=0.6)
    _camera(pl, azim=-100.0, elev=26.0, zoom=1.2)
    shot(pl, "seam-continuity", mode)


def quality_jacobian(mode: str) -> None:
    from shore.quality import scaled_jacobian

    pal = PALETTES[mode]
    pl = plotter(mode)
    actors = []
    for blk in _oh_mesh().blocks:
        grid = _grid(blk.nodes)
        grid.cell_data["scaled Jacobian"] = scaled_jacobian(blk.nodes).ravel(order="F")
        actors.append(
            _add(
                pl,
                _wedge(grid, -45.0),
                pal,
                scalars="scaled Jacobian",
                cmap=pal.quality,
                clim=(0.5, 1.0),
                lw=0.5,
                smooth=False,
            )
        )
    _add(pl, pv.read(str(_sphere_stl())), pal, pal.tint(pal.slate, 0.15), edges=False)
    _camera(pl, azim=-55.0, elev=28.0, zoom=1.0)
    actor = next(a for a in actors if a is not None)
    _bar(pl, actor, pal, "scaled Jacobian", fmt="%.1f", n_labels=6)
    shot(pl, "quality-jacobian", mode)


def ogrid_single(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    nodes = _close_j(_ogrid_mesh().blocks[0].nodes)
    nj = nodes.shape[1] - 1
    _add(pl, _grid(nodes[:, :, 0]), pal, pal.block(0), lw=0.5)
    for j in (0, nj // 2):
        _add(pl, _grid(nodes[:, j : j + 1, :]), pal, pal.sand, lw=0.5, smooth=False)
    _camera(pl, azim=-70.0, elev=20.0, zoom=1.15)
    shot(pl, "ogrid-single", mode)


def _wall(pl, pal, blk) -> None:
    """The k = 0 wall of a wing block in the neutral body colour."""
    _add(pl, _grid(_wing_view(blk.nodes[:, :, 0])), pal, pal.tint(pal.slate, 0.1), lw=0.4)


def _root_slice(pl, pal, blocks) -> None:
    """Root-station (j = 0) sheet of every block, coloured per block, and the body wall."""
    for k, blk in enumerate(blocks):
        sheet = _wing_view(blk.nodes[:, :1, :])
        _add(pl, _grid(sheet), pal, pal.block(k), lw=0.5, smooth=False)
    _wall(pl, pal, next(b for b in blocks if b.label == "body"))


def cgrid_wing(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    _root_slice(pl, pal, _cgrid_mesh().blocks)
    _camera(pl, azim=-108.0, elev=18.0, zoom=1.1)
    shot(pl, "cgrid-wing", mode)


def cgrid_tip_cap(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    blocks = {b.label: b for b in _cgrid_mesh("both").blocks}
    _wall(pl, pal, blocks["body"])
    for k, label in enumerate(("upper", "body", "lower")):
        nodes = blocks[label].nodes[:, -1:, :20]
        nodes = nodes[nodes[:, 0, 0, 0] < 1.35]
        _add(pl, _grid(_wing_view(nodes)), pal, pal.tint(pal.block(k), 0.5), lw=0.4, smooth=False)
    caps = [b for label, b in blocks.items() if label.startswith("cap_jhi")]
    for k, blk in enumerate(caps):
        tip = blk.nodes[:, :, 0]
        tip = _wing_view(tip[tip[:, 0, 0] < 1.35])
        _add(pl, _grid(tip), pal, pal.block(k + 2), lw=0.6, smooth=False)
    _camera(pl, azim=-112.0, elev=26.0, zoom=1.2)
    shot(pl, "cgrid-tip-cap", mode)


def ch_hybrid(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    _root_slice(pl, pal, _ch_mesh().blocks)
    _camera(pl, azim=-120.0, elev=20.0, zoom=1.1)
    shot(pl, "ch-hybrid", mode)


def oh_hybrid(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    for k, blk in enumerate(_oh_mesh().blocks):
        colour = pal.block(k) if k < 6 else pal.tint(pal.block(k - 6), 0.4)
        _add(pl, _wedge(_grid(blk.nodes), -45.0), pal, colour, lw=0.5, smooth=False)
    _add(pl, pv.read(str(_sphere_stl())), pal, pal.tint(pal.slate, 0.15), edges=False)
    _camera(pl, azim=-55.0, elev=28.0, zoom=1.0)
    shot(pl, "oh-hybrid", mode)


def h_grid(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    duct, tube = _h_grids()
    _add(pl, _grid(duct.blocks[0].nodes), pal, pal.block(0), lw=0.5, smooth=False)
    nodes = _close_j(tube.blocks[0].nodes) + np.array([0.0, 0.4, 0.1])
    _add(pl, _grid(nodes), pal, pal.block(1), lw=0.5, smooth=False)
    _camera(pl, azim=-145.0, elev=24.0, zoom=1.1)
    shot(pl, "h-grid", mode)


def primitives(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode, (1500, 800))
    box, annulus, caps = _primitives()
    _add(pl, _grid(box.blocks[0].nodes), pal, pal.block(0), lw=0.5, smooth=False)
    _add(pl, _grid(annulus.blocks[0].nodes), pal, pal.block(1), lw=0.5, smooth=False)
    for k, blk in enumerate(caps.blocks):
        _add(pl, _grid(blk.nodes), pal, pal.block(2 + k), lw=0.5, smooth=False)
    _camera(pl, azim=-70.0, elev=28.0, zoom=1.2)
    shot(pl, "primitives", mode)


def chimera_assembly(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    near, background = _chimera_meshes()
    for blk in background.blocks:
        grid = _grid(blk.nodes)
        part = _keep(grid, grid.cell_centers().points[:, 1] > 0.0)
        grey = pal.tint(pal.slate, 0.55 if mode == "light" else 0.1)
        _add(pl, part, pal, grey, lw=0.5, smooth=False)
    for k, blk in enumerate(near.blocks):
        _add(pl, _wedge(_grid(blk.nodes), -90.0), pal, pal.block(k), lw=0.5, smooth=False)
    _add(pl, pv.read(str(_sphere_stl(0.5))), pal, pal.tint(pal.slate, 0.15), edges=False)
    _camera(pl, azim=-65.0, elev=24.0, zoom=1.1)
    shot(pl, "chimera-assembly", mode)


def marker_boxes(mode: str) -> None:
    pal = PALETTES[mode]
    pl = plotter(mode)
    path, boxes = _marker_boxes()
    cubes = pv.merge([_grid(b.vertices.reshape(2, 2, 2, 3, order="F")) for b in boxes])
    _add(pl, pv.read(str(path)), pal, pal.tint(pal.slate, 0.1), edges=False)
    _add(pl, cubes, pal, pal.sunset, opacity=0.25, lw=0.8, smooth=False)
    _camera(pl, azim=-60.0, elev=28.0, zoom=1.0)
    shot(pl, "marker-boxes", mode)


FIGURES = {
    "sphere-blocks": sphere_blocks,
    "sphere-topdown": sphere_topdown,
    "marching-layers": marching_layers,
    "seam-continuity": seam_continuity,
    "quality-jacobian": quality_jacobian,
    "ogrid-single": ogrid_single,
    "cgrid-wing": cgrid_wing,
    "cgrid-tip-cap": cgrid_tip_cap,
    "ch-hybrid": ch_hybrid,
    "oh-hybrid": oh_hybrid,
    "h-grid": h_grid,
    "primitives": primitives,
    "chimera-assembly": chimera_assembly,
    "marker-boxes": marker_boxes,
}
