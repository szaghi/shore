"""pygfx/wgpu implementation of the GUI viewport (ADR 0006, 0037, 0038).

The only module of SHORE that imports pygfx, wgpu or rendercanvas.

Float precision: positions are kept in float64 on the CPU and uploaded as
``float32(p - origin)``; the origin follows the orbit target and moves
(re-uploading everything, "rebase") when the target drifts more than
``rebase_factor`` view widths from it, or on ``look_at_point``. Coordinates
reported by ``pick`` come from the float64 arrays (ADR 0006, spike item 3).

Drawn objects: per block its six boundary faces (quads) and twelve edges;
per slice plane its quads; per geometry its triangles, curves and points.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pygfx as gfx

from shore.gui.viewport import (
    FACES,
    GeometryPick,
    PickResult,
    SliceKey,
    colormap,
    edge_segments,
    face_grids,
    quad_cell,
    quad_centre,
    quad_indices,
    slice_cell,
    slice_grid,
)


def create_canvas(
    parent: Any = None, *, offscreen: bool = False, size: tuple[int, int] = (800, 600)
) -> Any:
    """A rendercanvas canvas: a Qt widget, or an offscreen canvas for tests."""
    if offscreen:
        from rendercanvas.offscreen import RenderCanvas

        return RenderCanvas(size=size)
    from rendercanvas.qt import QRenderWidget

    return QRenderWidget(parent, update_mode="ondemand")


def _normals(grid: np.ndarray) -> np.ndarray:
    da = np.gradient(grid, axis=0)
    db = np.gradient(grid, axis=1)
    n = np.cross(da, db)
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    return n / np.where(norm > 0.0, norm, 1.0)


def _quad_mesh(grids: list[np.ndarray], colour: tuple[float, ...]) -> tuple[gfx.Mesh, np.ndarray]:
    """One pickable quad mesh of node grids; its float64 positions."""
    pos, nrm, idx, base = [], [], [], 0
    for g in grids:
        na, nb, _ = g.shape
        pos.append(g.reshape(-1, 3))
        nrm.append(_normals(g).reshape(-1, 3))
        idx.append(quad_indices(na, nb, base))
        base += na * nb
    pos64 = np.concatenate(pos).astype(np.float64)
    n_quads = sum(len(i) for i in idx)
    geometry = gfx.Geometry(
        positions=pos64.astype(np.float32),
        normals=np.concatenate(nrm).astype(np.float32),
        indices=np.concatenate(idx),
        colors=np.tile(np.asarray(colour, np.float32), (n_quads, 1)),
    )
    mesh = gfx.Mesh(geometry, gfx.MeshPhongMaterial(color_mode="face", pick_write=True))
    return mesh, pos64


@dataclass
class _Item:
    """A drawn object: its pygfx objects and their float64 positions."""

    key: Any
    kind: str  # "block", "slice", "geometry"
    objects: list[gfx.WorldObject]
    pos64: list[np.ndarray]  # positions of each object, same order
    nodes: np.ndarray | None = None  # block / slice: the block's nodes
    colour: tuple[float, float, float, float] = (0.8, 0.8, 0.8, 1.0)
    display: Any = None  # geometry: its Display
    curve_offsets: np.ndarray = field(default_factory=lambda: np.zeros(1, np.int64))


class PygfxViewport:
    """pygfx/wgpu implementation of ``Viewport``."""

    def __init__(self, canvas: Any, *, rebase_factor: float = 100.0) -> None:
        self.canvas = canvas
        self.rebase_factor = rebase_factor
        self.rebase_count = 0
        self.renderer = gfx.WgpuRenderer(canvas)
        self.scene = gfx.Scene()
        self.scene.add(gfx.Background.from_color("#f4f5f7", "#d9dce1"))
        self.scene.add(gfx.AmbientLight(intensity=0.5))
        self.camera = gfx.PerspectiveCamera(45, 4 / 3)
        self.camera.add(gfx.DirectionalLight(intensity=2.0))
        self.scene.add(self.camera)
        self.groups = {kind: gfx.Group() for kind in ("block", "slice", "geometry")}
        self.content = gfx.Group()
        self.content.add(*self.groups.values())
        self.scene.add(self.content)
        self.controller = gfx.OrbitController(self.camera, register_events=self.renderer)
        self._items: list[_Item] = []
        self._by_object: dict[int, _Item] = {}
        self._origin = np.zeros(3)
        self._highlight = gfx.Line(
            gfx.Geometry(positions=np.zeros((5, 3), np.float32)),
            gfx.LineMaterial(thickness=3.0, color="#ff3b30", depth_test=False),
        )
        self._highlight.visible = False
        self.scene.add(self._highlight)
        canvas.request_draw(self._draw)

    # -- content -------------------------------------------------------------------

    def _remove(self, kind: str) -> None:
        self.groups[kind].clear()
        for item in [i for i in self._items if i.kind == kind]:
            self._items.remove(item)
            for obj in item.objects:
                self._by_object.pop(id(obj), None)
        self._highlight.visible = False

    def _add(self, item: _Item) -> None:
        for obj, pos in zip(item.objects, item.pos64, strict=True):
            obj.geometry.positions.data[:] = (pos - self._origin).astype(np.float32)
            self._by_object[id(obj)] = item
            self.groups[item.kind].add(obj)
        self._items.append(item)

    def set_blocks(self, blocks: Sequence[tuple[Any, np.ndarray]]) -> None:
        self._remove("block")
        for index, (key, nodes) in enumerate(blocks):
            hue = gfx.Color.from_hsv((index * 0.137) % 1.0, 0.35, 0.85)
            colour = (hue.r, hue.g, hue.b, 1.0)
            mesh, pos64 = _quad_mesh(face_grids(nodes), colour)
            edges64 = edge_segments(nodes).astype(np.float64)
            edges = gfx.Line(
                gfx.Geometry(positions=edges64.astype(np.float32)),
                gfx.LineSegmentMaterial(thickness=1.5, color="#202020"),
            )
            self._add(_Item(key, "block", [mesh, edges], [pos64, edges64], nodes, colour))
        self.request_draw()

    def set_slices(self, slices: Sequence[tuple[SliceKey, np.ndarray]]) -> None:
        """Replace the slice planes: ``(SliceKey, the block's nodes)``."""
        self._remove("slice")
        for key, nodes in slices:
            colour = (0.95, 0.85, 0.55, 1.0)
            mesh, pos64 = _quad_mesh([slice_grid(nodes, key.axis, key.index)], colour)
            self._add(_Item(key, "slice", [mesh], [pos64], nodes, colour))
        self.request_draw()

    def set_geometries(self, geometries: Sequence[tuple[str, Any]]) -> None:
        """Replace the geometries: ``(name, shore.geometry.base.Display)``."""
        self._remove("geometry")
        for name, display in geometries:
            objects: list[gfx.WorldObject] = []
            positions: list[np.ndarray] = []
            pts = np.asarray(display.points, np.float64)
            if len(display.triangles):
                tris = gfx.Mesh(
                    gfx.Geometry(
                        positions=pts.astype(np.float32),
                        indices=np.asarray(display.triangles, np.int32),
                    ),
                    gfx.MeshPhongMaterial(color="#9aa3ad", pick_write=True, flat_shading=True),
                )
                objects.append(tris)
                positions.append(pts)
            offsets = np.zeros(1, np.int64)
            if display.curves:
                # polylines separated by NaN rows: one Line, vertex index -> curve;
                # drawn over the faces they bound (CAD viewers do), so they stay pickable
                chunks, starts = [], [0]
                for c in display.curves:
                    chunks.append(np.vstack([c, np.full((1, 3), np.nan)]))
                    starts.append(starts[-1] + len(c) + 1)
                curves64 = np.concatenate(chunks)
                line = gfx.Line(
                    gfx.Geometry(positions=curves64.astype(np.float32)),
                    gfx.LineMaterial(
                        thickness=3.0, color="#1f4e79", pick_write=True, depth_test=False
                    ),
                )
                objects.append(line)
                positions.append(curves64)
                offsets = np.asarray(starts, np.int64)
            if len(display.vertices):
                verts = np.asarray(display.vertices, np.float64)
                points = gfx.Points(
                    gfx.Geometry(positions=verts.astype(np.float32)),
                    gfx.PointsMaterial(
                        size=9.0, color="#c0392b", pick_write=True, depth_test=False
                    ),
                )
                objects.append(points)
                positions.append(verts)
            item = _Item(name, "geometry", objects, positions, display=display)
            item.curve_offsets = offsets
            self._add(item)
        self.request_draw()

    def show(self, kind: str, visible: bool) -> None:
        """Show or hide every ``"block"``, ``"slice"`` or ``"geometry"`` object."""
        self.groups[kind].visible = visible
        self.request_draw()

    def visible(self, kind: str) -> bool:
        return bool(self.groups[kind].visible)

    def set_face_scalars(
        self, values: dict[Any, np.ndarray] | None, clim: tuple[float, float] = (0.0, 1.0)
    ) -> None:
        """Colour block and slice quads by value (keys: block keys, ``SliceKey``), or clear."""
        for item in self._items:
            if item.kind == "geometry":
                continue
            colours = item.objects[0].geometry.colors
            v = None if values is None else values.get(item.key)
            if v is None:
                colours.data[:] = np.asarray(item.colour, np.float32)
            else:
                colours.data[:] = colormap(np.asarray(v, np.float64), clim)
            colours.update_full()
        self.request_draw()

    def face_colours(self, key: Any) -> np.ndarray:
        """The uploaded per-quad RGBA of a block or slice (for tests)."""
        item = next(i for i in self._items if i.key == key and i.kind != "geometry")
        return np.asarray(item.objects[0].geometry.colors.data)

    @property
    def keys(self) -> list[Any]:
        return [i.key for i in self._items if i.kind == "block"]

    @property
    def slice_keys(self) -> list[SliceKey]:
        return [i.key for i in self._items if i.kind == "slice"]

    @property
    def geometry_names(self) -> list[str]:
        return [i.key for i in self._items if i.kind == "geometry"]

    def uploaded_vertices(self, key: Any) -> int:
        item = next(i for i in self._items if i.key == key)
        return len(item.objects[0].geometry.positions.data)

    def uploaded_positions(self, key: Any) -> np.ndarray:
        """float64 world positions of what was uploaded for *key*'s first object."""
        item = next(i for i in self._items if i.key == key)
        return np.asarray(item.objects[0].geometry.positions.data, np.float64) + self._origin

    # -- origin (float32 rebase) ------------------------------------------------------

    def _set_origin(self, origin: np.ndarray, *, move_camera: bool = True) -> None:
        delta = np.asarray(origin, np.float64) - self._origin
        self._origin = np.asarray(origin, np.float64)
        for item in self._items:
            for obj, pos in zip(item.objects, item.pos64, strict=True):
                buf = obj.geometry.positions
                buf.data[:] = (pos - self._origin).astype(np.float32)
                buf.update_full()
        if move_camera:
            # the controller's target setter calls camera.look_at: keep the pose
            self.camera.local.position = np.asarray(self.camera.local.position) - delta
            self.controller.target = tuple(np.asarray(self.controller.target) - delta)
        self.rebase_count += 1

    def view_width(self) -> float:
        tgt = np.asarray(self.controller.target, np.float64)
        dist = float(np.linalg.norm(np.asarray(self.camera.local.position) - tgt))
        return 2.0 * dist * np.tan(np.radians(float(self.camera.fov)) / 2.0)

    def maybe_rebase(self) -> bool:
        if not self._items:
            return False
        tgt = np.asarray(self.controller.target, np.float64)
        width = self.view_width()
        if width <= 0.0 or float(np.linalg.norm(tgt)) <= self.rebase_factor * width:
            return False
        self._set_origin(self._origin + tgt)
        return True

    # -- camera --------------------------------------------------------------------

    def show_all(self) -> None:
        pts = [p[np.isfinite(p).all(axis=1)] for i in self._items for p in i.pos64]
        pts = [p for p in pts if len(p)]
        if pts:
            allp = np.concatenate(pts)
            centre = 0.5 * (allp.min(axis=0) + allp.max(axis=0))
            self._set_origin(centre, move_camera=False)
        self.camera.depth_range = None
        self.camera.show_object(self.content, view_dir=(-1.0, -1.0, -1.0), up=(0.0, 0.0, 1.0))
        if pts:
            self.controller.target = (0.0, 0.0, 0.0)
        self.request_draw()

    def look_at_point(
        self, xyz: np.ndarray, distance: float, view_dir: np.ndarray | None = None
    ) -> None:
        self._set_origin(np.asarray(xyz, np.float64), move_camera=False)
        d = np.array([-0.3, -0.4, -0.866]) if view_dir is None else np.asarray(view_dir, float)
        d = d / np.linalg.norm(d)
        up = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
        self.camera.world.reference_up = up
        self.camera.local.position = -distance * d
        self.camera.look_at((0.0, 0.0, 0.0))
        self.camera.depth_range = (distance * 1e-3, distance * 1e3)
        self.controller.target = (0.0, 0.0, 0.0)
        self.request_draw()

    def camera_position(self) -> np.ndarray:
        """The camera's float64 world position."""
        return np.asarray(self.camera.local.position, np.float64) + self._origin

    def project(self, xyz: np.ndarray) -> tuple[float, float]:
        p = np.append(np.asarray(xyz, np.float64) - self._origin, 1.0)
        clip = np.asarray(self.camera.camera_matrix, np.float64) @ p
        ndc = clip[:3] / clip[3]
        w, h = self.renderer.logical_size
        return (ndc[0] + 1.0) * 0.5 * w, (1.0 - ndc[1]) * 0.5 * h

    def pixel_size(self, xyz: np.ndarray) -> float:
        """World size of one pixel at *xyz* (perspective)."""
        dist = float(np.linalg.norm(np.asarray(xyz, np.float64) - self.camera_position()))
        _, h = self.renderer.logical_size
        return 2.0 * dist * np.tan(np.radians(float(self.camera.fov)) / 2.0) / float(h)

    # -- picking -------------------------------------------------------------------

    def pick(self, x: float, y: float) -> PickResult | GeometryPick | None:
        info = self.renderer.get_pick_info((x, y))
        wo = info.get("world_object")
        item = self._by_object.get(id(wo)) if wo is not None else None
        if item is None:
            self._highlight.visible = False
            self.request_draw()
            return None
        if item.kind == "geometry":
            return self._pick_geometry(item, wo, info)
        if "face_index" not in info:
            return None
        quad = int(info["face_index"])
        shape = tuple(int(n) for n in item.nodes.shape[:3])
        if item.kind == "block":
            face, cell = quad_cell(shape, quad)
            xyz = quad_centre(item.nodes, face, cell)
            self._outline(face_grids(item.nodes)[FACES.index(face)], cell, face)
            return PickResult(item.key, face, cell, xyz)
        key: SliceKey = item.key
        cell = slice_cell(shape, key.axis, key.index, quad)
        grid = slice_grid(item.nodes, key.axis, key.index)
        a = "ijk".index(key.axis)
        p, q = (c for d, c in enumerate(cell) if d != a)
        xyz = grid[p : p + 2, q : q + 2].reshape(-1, 3).mean(axis=0)
        self._outline_quad(grid[p : p + 2, q : q + 2])
        return PickResult(key, f"{key.axis}={key.index}", cell, xyz)

    def _pick_geometry(self, item: _Item, wo: Any, info: dict) -> GeometryPick | None:
        from shore.geometry.base import EntityKind, EntityRef

        d = item.display
        self._highlight.visible = False
        if isinstance(wo, gfx.Mesh):
            t = int(info["face_index"])
            coord = np.asarray(info.get("face_coord", (1 / 3, 1 / 3, 1 / 3)), np.float64)
            corners = np.asarray(d.points, np.float64)[np.asarray(d.triangles)[t]]
            xyz = (coord[:3, None] * corners).sum(axis=0) / coord[:3].sum()
            entity = EntityRef(d.surface_kind, int(d.triangle_entity[t]))
            return GeometryPick(item.key, entity, xyz)
        if isinstance(wo, gfx.Points):
            v = int(info["vertex_index"])
            return GeometryPick(item.key, EntityRef(EntityKind.POINT, v), d.vertices[v].copy())
        v = int(info["vertex_index"])
        c = int(np.searchsorted(item.curve_offsets, v, side="right") - 1)
        curve = d.curves[c]
        local = min(v - int(item.curve_offsets[c]), len(curve) - 1)
        return GeometryPick(item.key, EntityRef(EntityKind.CURVE, c), curve[local].copy())

    def _outline(self, grid: np.ndarray, cell: tuple[int, int, int], face: str) -> None:
        f = FACES.index(face)
        i, j, k = cell
        a, b = {0: (j, k), 1: (j, k), 2: (i, k), 3: (i, k), 4: (i, j), 5: (i, j)}[f]
        self._outline_quad(grid[a : a + 2, b : b + 2])

    def highlight_polyline(self, points: np.ndarray) -> None:
        """Outline *points* (float64, (n, 3)): a selected edge or vertex."""
        pts = np.asarray(points, np.float64).reshape(-1, 3)
        if len(pts) == 1:  # a vertex: a small cross of the view's pixel size
            h = 6.0 * self.pixel_size(pts[0])
            pts = pts[0] + np.array([[-h, 0, 0], [h, 0, 0], [0, 0, 0], [0, -h, 0], [0, h, 0]])
        geometry = gfx.Geometry(positions=(pts - self._origin).astype(np.float32))
        self._highlight.geometry = geometry
        self._highlight.visible = True
        self.request_draw()

    def _outline_quad(self, q: np.ndarray) -> None:
        quad = np.array([q[0, 0], q[1, 0], q[1, 1], q[0, 1]])
        loop = np.vstack([quad, quad[:1]]) - self._origin
        self._highlight.geometry = gfx.Geometry(positions=loop.astype(np.float32))
        self._highlight.visible = True
        self.request_draw()

    def on_pick(self, callback: Any) -> None:
        """Call ``callback(pick result or None, shift)`` on a double click in the view."""

        def handler(event: Any) -> None:
            shift = "Shift" in (getattr(event, "modifiers", None) or ())
            callback(self.pick(event.x, event.y), shift)

        self.renderer.add_event_handler(handler, "double_click")

    # -- drawing -------------------------------------------------------------------

    def _draw(self) -> None:
        self.maybe_rebase()
        self.renderer.render(self.scene, self.camera)

    def request_draw(self) -> None:
        self.canvas.request_draw()

    def draw_now(self) -> np.ndarray:
        """Render synchronously and return the image (offscreen canvases)."""
        return np.asarray(self.canvas.draw())
