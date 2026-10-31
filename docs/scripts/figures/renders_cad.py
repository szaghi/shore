"""3-D renders of the CAD, collar and blocking figures (PyVista, off-screen).

The wing-body figures reproduce the phase-7 acceptance (ADR 0051,
``tests/spec/test_wing_body.py``): the half model is written as STEP with
build123d and blocked through ``shore.project.Project`` with the same
journaled commands. The model is built once per run (``functools.cache``).
Needs ``shore[cad]`` and ``shore[vis]``.
"""

from __future__ import annotations

import functools
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from theme import PALETTES, plotter, shot

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")
FUSELAGE_BLOCKS = [(x, y, 0.4) for x in (-1.5, 0.5, 2.0) for y in (-0.55, 0.0, 0.55)]
WING_BLOCK = (0.5, 0.0, 2.0)
BASE_EDGES = ("j_lo_k_lo", "j_hi_k_lo", "i_lo_k_lo", "i_hi_k_lo")

# -- the wing-body model of ADR 0051 -------------------------------------------------


def wing_body(path: Path, span: float = 3.0, radius: float = 0.5, length: float = 5.0) -> Path:
    """Write the half model as STEP to *path* (as ``tests/spec/test_wing_body.py``)."""
    from build123d import (
        Box,
        Cylinder,
        Line,
        Pos,
        Sphere,
        Spline,
        Wire,
        export_step,
        extrude,
        make_face,
    )

    from shore.profiles import naca

    body = Cylinder(radius, length, rotation=(0, 90, 0))
    body = body + Pos(-length / 2, 0, 0) * Sphere(radius) + Pos(length / 2, 0, 0) * Sphere(radius)
    sec = naca("0012", 40, closed_te=False)
    le = int(np.argmin(sec[:, 0]))
    upper = [(float(x), float(y), 0.0) for x, y in sec[: le + 1]]
    lower = [(float(x), float(y), 0.0) for x, y in sec[le:]]
    wire = Wire([Spline(*upper), Spline(*lower), Line(lower[-1], upper[0])])
    wing = extrude(make_face(wire), amount=span)
    export_step((body + wing) & (Pos(0, 0, 10) * Box(40, 40, 20)), str(path))
    return path


def junction_edges(geometry, radius: float) -> list[int]:
    """Indices of the CAD edges lying on the fuselage and on the wing."""
    from shore.geometry.base import EntityKind

    out = []
    for k in range(len(geometry.entities(EntityKind.CURVE))):
        pts = geometry.display().curves[k]
        on_fuselage = np.allclose(np.hypot(pts[:, 1], pts[:, 2]), radius, atol=1e-3)
        at_wing = np.all(np.abs(pts[:, 1]) < 0.07) and np.all(
            (pts[:, 0] > -0.01) & (pts[:, 0] < 1.01)
        )
        if on_fuselage and at_wing:
            out.append(k)
    return out


def _block_at(p, point) -> str:
    b = p.blocking
    for label in b.labels():
        pts = b.points[list(b.block(label))]
        if np.all(pts.min(axis=0) < point) and np.all(point < pts.max(axis=0)):
            return label
    raise AssertionError(f"no block holds {point}")


def _cut(p, axis: str, c: float, other) -> None:
    d = "ijk".index(axis)
    point = np.array(other, float)
    point[d] = c
    label = _block_at(p, point)
    pts = p.blocking.points[list(p.blocking.block(label))]
    a, z = pts.min(axis=0)[d], pts.max(axis=0)[d]
    p.split_blocking(label, axis, (c - a) / (z - a))


@dataclass
class WingBody:
    """The wing-body stages: the geometry, blockings and grids along the way."""

    step: Path
    geometry: object
    junction: list[int]
    split: object  # blocking after the cuts (hole still there)
    holed: object  # blocking after the delete (hole open)
    ogrid: object  # blocking after the O-grid
    premesh: list  # [(label, {face: bc}, nodes)] after premesh(exact=True)
    final: list  # the same after orthogonalize(0.5) and wall_layer(2e-3)
    hole: list[str]


def _grid(p) -> list:
    return [
        (b.label, {f.name: f.bc.name for f in b.faces()}, b.nodes.copy()) for b in p.mesh.blocks
    ]


@functools.cache
def wing_body_stages() -> WingBody:
    """Build ADR 0051's blocking once, keeping every stage."""
    from shore.project import Project

    step = wing_body(Path(tempfile.mkdtemp(prefix="shore-fig-")) / "wb.step")
    p = Project()
    p.load_geometry(step, name="wb")
    p.create_blocking(53, 33, 25, vmin=(-6, -4, 0), vmax=(7, 4, 6))
    p.symmetry_plane((0, 0, 0), (0, 0, 1))  # ADR 0053
    corner = (-5.99, -3.99, 0.01)
    for axis, cuts in (("i", (-3.0, -0.3, 1.3, 3.0)), ("j", (-0.8, -0.3, 0.3, 0.8)),
                       ("k", (0.8, 3.3))):  # fmt: skip
        for c in cuts:
            _cut(p, axis, c, corner)
    split = p.blocking.copy()
    hole = [_block_at(p, np.array(c)) for c in [*FUSELAGE_BLOCKS, WING_BLOCK]]
    b = p.blocking
    owners: dict = {}
    for label in hole:
        for f in FACES:
            owners.setdefault(b.entity(label, f)[0], []).append((label, f))
    walls = {k: v[0] for k, v in owners.items()
             if len(v) == 1 and not np.allclose(b.points[list(k[1:])][:, 2], 0.0)}  # fmt: skip
    for label, f in walls.values():
        p.associate_blocking(label, f, "wb", "body")
    geometry = p.geometries["wb"].geometry
    je = junction_edges(geometry, 0.5)
    for e in BASE_EDGES:
        p.associate_blocking(hole[-1], e, "wb", "curve", je)
    p.delete_blocks(hole)
    holed = p.blocking.copy()
    layered = [(label, f) for label in p.blocking.labels() for f in FACES
               if p.blocking.entity(label, f)[0] in walls]  # fmt: skip
    p.ogrid(faces=layered, factor=0.7, nodes=17)
    ogrid = p.blocking.copy()
    p.premesh(exact=True)
    premesh = _grid(p)
    p.orthogonalize(0.5)
    p.wall_layer(2e-3)
    return WingBody(step, geometry, je, split, holed, ogrid, premesh, _grid(p), hole)


# -- drawing helpers -----------------------------------------------------------------


def _up(points: np.ndarray) -> np.ndarray:
    """The display frame: the model's y (the section's thickness) up.

    A proper rotation, ``(x, y, z) -> (x, -z, y)``: the span (+z) runs to -y,
    the symmetry plane z = 0 becomes y = 0, and the theme's lights (set for z
    up) light the wing and fuselage from above.
    """
    p = np.asarray(points, dtype=np.float64)
    return np.stack([p[..., 0], -p[..., 2], p[..., 1]], axis=-1)


def _sheet(grid: np.ndarray):
    """A PyVista structured surface of an (na, nb, 3) node grid."""
    import pyvista as pv

    g = np.ascontiguousarray(_up(grid))
    return pv.StructuredGrid(g[..., 0], g[..., 1], g[..., 2])


def _face(nodes: np.ndarray, name: str) -> np.ndarray:
    return np.take(nodes, 0 if name.endswith("lo") else -1, axis="ijk".index(name[0]))


def _poly(points: np.ndarray, closed: bool = False):
    import pyvista as pv

    pts = _up(points)
    if closed:
        pts = np.vstack([pts, pts[:1]])
    return pv.lines_from_points(pts)


def _tube(points: np.ndarray, radius: float, closed: bool = False):
    pts = np.asarray(points, dtype=np.float64)
    keep = np.r_[True, np.linalg.norm(np.diff(pts, axis=0), axis=1) > 1e-9 * radius]
    return _poly(pts[keep], closed).tube(radius=radius, n_sides=16)


def _surface_kw(colour, mode: str, *, edges: bool = True, width: float = 0.6, opacity=1.0):
    pal = PALETTES[mode]
    return {
        "color": colour,
        "show_edges": edges,
        "edge_color": pal.edge,
        "line_width": width,
        "smooth_shading": True,
        "specular": 0.25,
        "specular_power": 20.0,
        "ambient": 0.18 if mode == "light" else 0.25,
        "diffuse": 0.85,
        "opacity": opacity,
    }


def _trimesh_pv(points: np.ndarray, triangles: np.ndarray):
    import pyvista as pv

    tris = np.asarray(triangles, dtype=np.int64)
    faces = np.hstack([np.full((len(tris), 1), 3), tris]).ravel()
    return pv.PolyData(_up(points), faces)


def _camera(pl, focal, direction, distance: float, zoom: float = 1.0):
    """Look at *focal* (model frame) from *direction* (display frame, z up)."""
    d = np.asarray(direction, float)
    d /= np.linalg.norm(d)
    f = _up(focal)
    up = (0.0, 0.0, 1.0)
    pl.camera_position = [tuple(f + distance * d), tuple(f), up]
    pl.camera.view_angle = 30.0
    pl.camera.zoom(zoom)


# the elevated 3/4 view in the display frame: from the front, the wing side, above
VIEW = (-0.9, -1.25, 0.95)


def _wall_faces(grid: list):
    """(block index, label, face name, node grid) of every WALL face."""
    out = []
    for n, (label, bcs, nodes) in enumerate(grid):
        for name, bc in bcs.items():
            if bc == "WALL":
                out.append((n, label, name, _face(nodes, name)))
    return out


def _plane_faces(grid: list, box=None):
    """(block index, node grid) of every face on the symmetry plane z = 0."""
    out = []
    for n, (_, _, nodes) in enumerate(grid):
        for name in FACES:
            g = _face(nodes, name)
            # the faces on the declared symmetry plane z = 0 (ADR 0053)
            if np.abs(g[..., 2]).max() > 1e-9:
                continue
            if box is not None:
                lo, hi = g.reshape(-1, 3).min(axis=0), g.reshape(-1, 3).max(axis=0)
                if np.any(lo[:2] < box[0]) or np.any(hi[:2] > box[1]):
                    continue
            out.append((n, g))
    return out


def _plane_window(pl, mode: str, sheet, centre, radii, opacity: float = 0.35) -> None:
    """Draw *sheet* (on the symmetry plane) inside an ellipse of the model's
    (x, y): a translucent fill and its grid lines, both cut smoothly."""
    pal = PALETTES[mode]
    pts = np.asarray(sheet.points)  # display frame: model x = x, model y = z
    sheet["r"] = ((pts[:, 0] - centre[0]) / radii[0]) ** 2 + (
        (pts[:, 2] - centre[1]) / radii[1]
    ) ** 2
    fill = sheet.clip_scalar(scalars="r", value=1.0)
    if fill.n_cells == 0:
        return
    lines = sheet.extract_all_edges().clip_scalar(scalars="r", value=1.0)
    pl.add_mesh(
        fill, color=pal.tint(pal.sky, 0.35), opacity=opacity, ambient=0.8, diffuse=0.2,
        scalars=None, show_scalar_bar=False,
    )  # fmt: skip
    if lines.n_cells:
        pl.add_mesh(
            lines, color=pal.tint(pal.ocean, 0.25), line_width=0.7, opacity=0.6,
            scalars=None, show_scalar_bar=False,
        )  # fmt: skip


def _volume(nodes: np.ndarray):
    import pyvista as pv

    g = np.ascontiguousarray(_up(nodes))
    return pv.StructuredGrid(g[..., 0], g[..., 1], g[..., 2])


def _station(nodes: np.ndarray, z: float, xr, yr):
    """The cut of a block by the plane z = *z* (model frame), inside x in *xr*, y in *yr*."""
    if not (nodes[..., 2].min() < z < nodes[..., 2].max()):
        return None
    if nodes[..., 0].max() < xr[0] or nodes[..., 0].min() > xr[1]:
        return None
    cut = _volume(nodes).slice(normal=(0.0, 1.0, 0.0), origin=(0.0, -z, 0.0))
    if cut.n_cells == 0:
        return None
    ctr = np.asarray(cut.cell_centers().points)
    keep = np.nonzero((ctr[:, 0] > xr[0]) & (ctr[:, 0] < xr[1])
                      & (ctr[:, 2] > yr[0]) & (ctr[:, 2] < yr[1]))[0]  # fmt: skip
    return (
        cut.extract_cells(keep).extract_surface(algorithm="dataset_surface") if len(keep) else None
    )


def _soft(pal, colour: str, amount: float):
    """*colour* lightened for a fill: towards white in light mode, only a little
    towards the page in dark mode (a dark fill would hide its grid lines)."""
    return pal.tint(colour, amount if pal.mode == "light" else 0.4 * amount)


def _colour_of(labels: list[str]) -> dict[str, int]:
    """A palette index per block label, stable across figures."""
    return {label: k for k, label in enumerate(sorted(labels))}


# -- wing-body figures ---------------------------------------------------------------


def hero_wing_body(mode: str) -> None:
    pal = PALETTES[mode]
    wb = wing_body_stages()
    pl = plotter(mode, (1800, 1150))
    walls = _wall_faces(wb.final)
    for k, (_, _, _, g) in enumerate(walls):
        pl.add_mesh(_sheet(g), **_surface_kw(pal.block(k), mode, width=0.5))
    # the symmetry plane round the body: the grid as a thin translucent sheet
    # cut to an ellipse, so it fades into the page
    for _, g in _plane_faces(wb.final):
        _plane_window(pl, mode, _sheet(g), (0.5, 0.0), (4.6, 2.3))
    # a spanwise station of the volume through the wing: the O-grid round the section
    for _, bcs, nodes in wb.final:
        if "WALL" not in bcs.values():
            continue  # the O-grid layer only: the blocks on the wall
        cut = _station(nodes, 2.1, (-5.0, 5.0), (-5.0, 5.0))
        if cut is not None:
            pl.add_mesh(
                cut,
                color=_soft(pal, pal.sand, 0.35),
                show_edges=True,
                edge_color=pal.edge,
                line_width=0.6,
                ambient=0.55,
                diffuse=0.5,
            )
    _camera(pl, (0.2, 0.0, 1.3), VIEW, 11.5)
    shot(pl, "hero-wing-body", mode)


def blocking_premesh(mode: str) -> None:
    pal = PALETTES[mode]
    wb = wing_body_stages()
    pl = plotter(mode, (2100, 1250))
    labels = _colour_of([label for label, _, _ in wb.premesh])
    for _, label, _, g in _wall_faces(wb.premesh):
        pl.add_mesh(_sheet(g), **_surface_kw(pal.block(labels[label]), mode, width=0.6))
    for _, g in _plane_faces(wb.premesh):
        _plane_window(pl, mode, _sheet(g), (0.5, 0.0), (5.6, 3.2), opacity=0.5)
    _camera(pl, (0.3, 0.0, 1.0), VIEW, 12.0)
    shot(pl, "blocking-premesh", mode)


def blocking_skeleton(mode: str) -> None:
    import pyvista as pv

    pal = PALETTES[mode]
    wb = wing_body_stages()
    b = wb.ogrid
    pl = plotter(mode, (1700, 1100))
    d = wb.geometry.display()
    pl.add_mesh(
        _trimesh_pv(d.points, d.triangles),
        color=_soft(pal, pal.slate, 0.35),
        smooth_shading=True,
        specular=0.2,
        ambient=0.25,
        opacity=0.6,
    )
    # skeleton edges: the O-grid layer in sea, the rest in ocean
    layer = set(b.labels()) - set(wb.holed.labels())
    segs = {"layer": [], "outer": []}
    for label, _, pts in b.polylines():
        segs["layer" if label in layer else "outer"].append(pts[[0, -1]])
    for key, colour, width in (("outer", pal.tint(pal.muted, 0.3), 1.2), ("layer", pal.ocean, 2.6)):
        lines = np.concatenate(segs[key])
        n = len(lines) // 2
        cells = np.hstack([np.full((n, 1), 2), np.arange(2 * n).reshape(n, 2)]).ravel()
        pl.add_mesh(pv.PolyData(_up(lines), lines=cells), color=colour, line_width=width)
    centres = np.array([c for _, _, c in b.centres("block")])
    balls = pv.PolyData(_up(centres)).glyph(geom=pv.Sphere(radius=0.05), scale=False)
    pl.add_mesh(balls, color=pal.coral, smooth_shading=True, specular=0.3)
    _camera(pl, (0.3, 0.0, 1.4), VIEW, 12.5)
    shot(pl, "blocking-skeleton", mode)


def _xcut(nodes: np.ndarray, x: float, yr, zr):
    """The cut of a block by the plane x = *x* (model frame), inside y in *yr*, z in *zr*."""
    if not (nodes[..., 0].min() < x < nodes[..., 0].max()):
        return None
    pts = nodes.reshape(-1, 3)
    if pts[:, 1].max() < yr[0] or pts[:, 1].min() > yr[1]:
        return None
    if pts[:, 2].max() < zr[0] or pts[:, 2].min() > zr[1]:
        return None
    cut = _volume(nodes).slice(normal=(1.0, 0.0, 0.0), origin=(x, 0.0, 0.0))
    if cut.n_cells == 0:
        return None
    box = (x - 1.0, x + 1.0, -zr[1], -zr[0], yr[0], yr[1])  # display frame
    return cut.clip_box(box, invert=False, crinkle=True).extract_surface(
        algorithm="dataset_surface"
    )


def wall_layer(mode: str) -> None:
    pal = PALETTES[mode]
    wb = wing_body_stages()
    pl = plotter(mode, (1600, 1250))
    labels = _colour_of([label for label, _, _ in wb.final])
    x = 0.3
    yr, zr = (-0.6, 0.6), (0.2, 1.2)
    for label, _, nodes in wb.final:
        cut = _xcut(nodes, x, yr, zr)
        if cut is not None and cut.n_cells:
            pl.add_mesh(
                cut,
                color=_soft(pal, pal.block(labels[label]), 0.6),
                show_edges=True,
                edge_color=pal.edge,
                line_width=0.7,
                ambient=0.7,
                diffuse=0.3,
            )
    # the wall in the cut plane: the body's section, drawn over the grid
    d = wb.geometry.display()
    body = _trimesh_pv(d.points, d.triangles)
    section = body.slice(normal=(1.0, 0.0, 0.0), origin=(x, 0, 0))
    inside = body.clip(normal=(-1.0, 0.0, 0.0), origin=(x, 0.0, 0.0))  # the body behind the cut
    pl.add_mesh(inside, color=_soft(pal, pal.slate, 0.4), ambient=0.6, diffuse=0.3)
    section.points[:, 0] -= 1e-3  # towards the camera
    pl.add_mesh(section, color=pal.ink, line_width=3.0)
    _camera(pl, (x, 0.02, 0.62), (-1.0, 0.0, 0.0), 2.0)
    pl.enable_parallel_projection()
    pl.camera.parallel_scale = 0.34
    shot(pl, "wall-layer", mode)


def cad_entities(mode: str) -> None:
    import pyvista as pv

    pal = PALETTES[mode]
    wb = wing_body_stages()
    d = wb.geometry.display()
    pl = plotter(mode, (1700, 1100))
    mesh = _trimesh_pv(d.points, d.triangles)
    for k in np.unique(d.triangle_entity):
        part = mesh.extract_cells(np.nonzero(d.triangle_entity == k)[0]).extract_surface(
            algorithm="dataset_surface"
        )
        pl.add_mesh(
            part.compute_normals(split_vertices=True, feature_angle=40.0),
            color=pal.block(int(k)),
            smooth_shading=True,
            specular=0.3,
            specular_power=20.0,
            ambient=0.2,
        )
    for k, c in enumerate(d.curves):
        if np.linalg.norm(np.ptp(c, axis=0)) > 1e-6:
            pl.add_mesh(_tube(_dense_curve(wb.geometry, k, c), 0.012), color=pal.ink)
    if len(d.vertices):
        pl.add_mesh(
            pv.PolyData(_up(d.vertices)),
            color=pal.ink,
            point_size=14.0,
            render_points_as_spheres=True,
        )
    _camera(pl, (0.1, 0.0, 1.1), VIEW, 10.5)
    shot(pl, "cad-entities", mode)


def _dense_curve(geometry, k: int, c: np.ndarray, per: int = 12) -> np.ndarray:
    """Polyline *c* of CAD curve *k* refined *per* times, projected back onto the curve."""
    from shore.geometry.base import EntityKind, EntityRef

    t = np.linspace(0.0, 1.0, per, endpoint=False)[:, None, None]
    pts = ((1 - t) * c[None, :-1] + t * c[None, 1:]).transpose(1, 0, 2).reshape(-1, 3)
    pts = np.vstack([pts, c[-1:]])
    return geometry.project(pts, EntityRef(EntityKind.CURVE, k), exact=True).points


# -- STL feature curves (ADR 0049) ---------------------------------------------------


def stl_features(mode: str) -> None:
    import pyvista as pv

    from shore.geometry.stl import StlGeometry
    from shore.profiles import naca
    from shore.shapes import extrude_section

    pal = PALETTES[mode]
    sec = naca("0012", n_pts=60, closed_te=False)
    wing = extrude_section(
        sec, span=3.0, n_stations=7, taper=0.5, sweep_deg=20.0, span_axis="z", end_caps=True
    )
    g = StlGeometry(wing)
    d = g.display()
    pl = plotter(mode, (1700, 1100))
    pl.add_mesh(
        _trimesh_pv(d.points, d.triangles),
        color=_soft(pal, pal.slate, 0.55),
        smooth_shading=False,
        specular=0.2,
        ambient=0.25,
    )
    for k, c in enumerate(d.curves):
        pl.add_mesh(
            _tube(c, 0.012), color=pal.block([0, 1, 3, 4, 5, 6][k % 6]), smooth_shading=True
        )
    balls = pv.PolyData(_up(d.vertices)).glyph(geom=pv.Sphere(radius=0.035), scale=False)
    pl.add_mesh(balls, color=pal.ink, smooth_shading=True)
    _camera(pl, (1.0, 0.0, 1.5), VIEW, 6.6)
    shot(pl, "stl-features", mode)


# -- intersection, surface march and collar (ADR 0008-0010, 0024) -------------------


@functools.cache
def _fuselage_wing(sharp: bool):
    """A cylinder fuselage along x (radius 0.5) and a NACA 0012 wing along +z
    from inside it, as STL geometries (the junction loop near z = 0.5)."""
    import trimesh

    from shore.geometry.stl import StlGeometry
    from shore.profiles import naca
    from shore.shapes import extrude_section

    fus = trimesh.creation.cylinder(radius=0.5, height=4.0, sections=256)
    fus.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [0, 1, 0]))
    # half a facet off: the trailing edge does not land on a facet edge
    fus.apply_transform(trimesh.transformations.rotation_matrix(0.37 * np.pi / 256, [1, 0, 0]))
    wing = extrude_section(
        naca("0012", n_pts=80, closed_te=sharp), span=1.6, n_stations=2, span_axis="z",
        end_caps=True,
    )  # fmt: skip
    return fus, wing, StlGeometry(fus), StlGeometry(wing)


SEED = np.array([0.5, 0.06, 0.5])


@functools.cache
def _loop() -> np.ndarray:
    from shore.surface.intersect import intersection_loop

    _, _, a, b = _fuselage_wing(False)
    return intersection_loop(a, b, SEED, 96, curvature=1.0)


def _bodies(pl, mode: str, sharp: bool, opacity: float = 1.0, clip: float | None = None) -> None:
    pal = PALETTES[mode]
    fus, wing, _, _ = _fuselage_wing(sharp)
    f = _trimesh_pv(fus.vertices, fus.faces)
    if clip is not None:
        f = f.clip_box((-clip, clip, -1, 1, -1, 1), invert=False).extract_surface(
            algorithm="dataset_surface"
        )
    w = _trimesh_pv(wing.vertices, wing.faces)
    for mesh, colour in ((f, _soft(pal, pal.slate, 0.45)), (w, pal.tint(pal.sky, 0.25))):
        pl.add_mesh(
            mesh.compute_normals(split_vertices=True, feature_angle=50.0),
            color=colour,
            smooth_shading=True,
            specular=0.25,
            specular_power=20.0,
            ambient=0.2,
            opacity=opacity,
        )


def intersection_curve(mode: str) -> None:
    pal = PALETTES[mode]
    loop = _loop()
    pl = plotter(mode, (1600, 1100))
    _bodies(pl, mode, False, clip=1.8)
    pl.add_mesh(_tube(loop, 0.012, closed=True), color=pal.coral, smooth_shading=True)
    _camera(pl, (0.5, 0.1, 0.55), (-0.9, -1.0, 0.95), 4.2)
    shot(pl, "intersection-curve", mode)


@functools.cache
def _marched() -> tuple[np.ndarray, np.ndarray]:
    from shore.surface.march import march_surface

    _, _, a, b = _fuselage_wing(False)
    loop = _loop()
    steps = 0.02 * 1.12 ** np.arange(15)
    on_a = march_surface(a, loop, steps, periodic=True)
    on_b = march_surface(b, loop[::-1], 0.02 * 1.12 ** np.arange(10), periodic=True)
    return on_a, on_b


def _closed(grid: np.ndarray) -> np.ndarray:
    return np.concatenate([grid, grid[:1]], axis=0)


def _lift(grid: np.ndarray, eps: float, *geometries) -> np.ndarray:
    """*grid* moved *eps* off the nearest of *geometries* along its outward
    normal, so a grid on an STL does not z-fight with it."""
    from shore.geometry.base import BODY

    pts = grid.reshape(-1, 3)
    best = [g.project(pts, BODY) for g in geometries]
    k = np.argmin(np.stack([b.distances for b in best]), axis=0)
    normals = np.stack([b.normals for b in best])[k, np.arange(len(pts))]
    return (pts + eps * normals).reshape(grid.shape)


def surface_march(mode: str) -> None:
    pal = PALETTES[mode]
    _, _, a, b = _fuselage_wing(False)
    on_a, on_b = _marched()
    pl = plotter(mode, (1600, 1100))
    _bodies(pl, mode, False, clip=1.8)
    pl.add_mesh(_sheet(_closed(_lift(on_a, 5e-3, a))), **_surface_kw(pal.sea, mode, width=0.7))
    pl.add_mesh(_sheet(_closed(_lift(on_b, 2e-3, b))), **_surface_kw(pal.sand, mode, width=0.7))
    pl.add_mesh(_tube(_loop(), 0.007, closed=True), color=pal.coral, smooth_shading=True)
    _camera(pl, (0.5, 0.1, 0.6), (-0.9, -1.0, 0.95), 3.6)
    shot(pl, "surface-march", mode)


@functools.cache
def _collar_blocks():
    from shore.collar import split_collar

    _, _, a, b = _fuselage_wing(True)
    steps = 0.005 * 1.12 ** np.arange(15)
    return split_collar(
        a, b, SEED, n_loop=96, steps_a=steps, steps_b=steps,
        steps_k=5e-5 * 1.2 ** np.arange(25), loop_curvature=1.0,
    )  # fmt: skip


def collar(mode: str) -> None:
    pal = PALETTES[mode]
    _, _, a, b = _fuselage_wing(True)
    blocks = _collar_blocks()
    pl = plotter(mode, (1600, 1150))
    _bodies(pl, mode, True, clip=1.8)
    for k, blk in enumerate(blocks):
        if blk.face("k_lo").bc.name != "WALL":
            continue  # the wake's virtual sheet: no wall
        wall = _lift(blk.nodes[:, :, 0], 1e-3, a, b)
        pl.add_mesh(_sheet(wall), **_surface_kw(pal.block(k), mode, width=0.5))
    # a j-slice of the volume at the leading edge, both halves
    for label in ("body_b", "body_a"):
        blk = next(x for x in blocks if x.label.endswith(label))
        row = blk.nodes[:, :, 0, 0].min(axis=0)
        j = int(np.argmin(row))
        pl.add_mesh(
            _sheet(blk.nodes[:, j]),
            color=_soft(pal, pal.sand, 0.15),
            show_edges=True,
            edge_color=pal.edge,
            line_width=0.6,
            ambient=0.5,
        )
    _camera(pl, (0.05, 0.0, 0.52), (-0.8, -0.9, 0.8), 0.75)
    shot(pl, "collar", mode)


# -- O-grid insertion (ADR 0043) -----------------------------------------------------


def blocking_ogrid(mode: str) -> None:
    from shore.mesh.blocking import box, ogrid

    pal = PALETTES[mode]
    before = box((-1, -1, -1), (1, 1, 1), 9, 9, 9)
    after = ogrid(before, ["b0"], exclude=[("b0", "k_lo")], factor=0.5, nodes=5)
    pl = plotter(mode, (2000, 900))
    shift = np.array([3.2, 0.0, 0.0])
    for n, (b, dx) in enumerate(((before, np.zeros(3)), (after, shift))):
        for k, (label, verts) in enumerate(b.blocks):
            core = label == "b0"
            colour = pal.slate if n == 0 else (pal.coral if core else pal.block(k - 1))
            pl.add_mesh(
                _hex_display(b.points[list(verts)] + dx),
                color=colour,
                opacity=0.8 if core and n == 1 else 0.22,
                smooth_shading=False,
                ambient=0.35,
                specular=0.2,
            )
        for _, _, pts in b.polylines():
            pl.add_mesh(_tube(pts[[0, -1]] + dx, 0.012), color=pal.ink, smooth_shading=True)
    _camera(pl, (1.6, -0.1, 0.0), (-0.45, -1.0, 0.6), 9.0)
    shot(pl, "blocking-ogrid", mode)


def _hex_display(corners: np.ndarray):
    """The hexahedron of *corners* in the display frame (the model's y up)."""
    return _hex(_up(corners))


def _hex(corners: np.ndarray):
    """An unstructured hexahedron of 8 skeleton corners (vertex index 4i + 2j + k)."""
    import pyvista as pv

    order = [0, 4, 6, 2, 1, 5, 7, 3]  # VTK hexahedron from the i-j-k bit order
    cells = np.array([8, *range(8)])
    return pv.UnstructuredGrid(cells, np.array([pv.CellType.HEXAHEDRON]), corners[order])


FIGURES = {
    "hero-wing-body": hero_wing_body,
    "stl-features": stl_features,
    "cad-entities": cad_entities,
    "intersection-curve": intersection_curve,
    "surface-march": surface_march,
    "collar": collar,
    "blocking-skeleton": blocking_skeleton,
    "blocking-ogrid": blocking_ogrid,
    "blocking-premesh": blocking_premesh,
    "wall-layer": wall_layer,
}
