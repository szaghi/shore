"""pygfx/wgpu implementation of the GUI viewport (ADR 0006, 0037).

The only module of SHORE that imports pygfx, wgpu or rendercanvas.

Float precision: positions are kept in float64 on the CPU and uploaded as
``float32(p - origin)``; the origin follows the orbit target and moves
(re-uploading every block, "rebase") when the target drifts more than
``rebase_factor`` view widths from it, or on ``look_at_point``. Coordinates
reported by ``pick`` come from the float64 arrays (ADR 0006, spike item 3).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pygfx as gfx

from shore.gui.viewport import (
    PickResult,
    colormap,
    edge_segments,
    face_grids,
    quad_cell,
    quad_centre,
    quad_indices,
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


@dataclass
class _Block:
    key: Any
    nodes: np.ndarray
    mesh: gfx.Mesh
    edges: gfx.Line
    pos64: np.ndarray
    edges64: np.ndarray
    colour: tuple[float, float, float, float]


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
        self.content = gfx.Group()
        self.scene.add(self.content)
        self.controller = gfx.OrbitController(self.camera, register_events=self.renderer)
        self._blocks: list[_Block] = []
        self._by_mesh: dict[int, _Block] = {}
        self._origin = np.zeros(3)
        self._highlight = gfx.Line(
            gfx.Geometry(positions=np.zeros((5, 3), np.float32)),
            gfx.LineMaterial(thickness=3.0, color="#ff3b30", depth_test=False),
        )
        self._highlight.visible = False
        self.scene.add(self._highlight)
        canvas.request_draw(self._draw)

    # -- content -------------------------------------------------------------------

    def set_blocks(self, blocks: Sequence[tuple[Any, np.ndarray]]) -> None:
        self.content.clear()
        self._blocks.clear()
        self._by_mesh.clear()
        self._highlight.visible = False
        for index, (key, nodes) in enumerate(blocks):
            grids = face_grids(nodes)
            pos, nrm, idx, base = [], [], [], 0
            for g in grids:
                na, nb, _ = g.shape
                pos.append(g.reshape(-1, 3))
                nrm.append(_normals(g).reshape(-1, 3))
                idx.append(quad_indices(na, nb, base))
                base += na * nb
            pos64 = np.concatenate(pos).astype(np.float64)
            hue = gfx.Color.from_hsv((index * 0.137) % 1.0, 0.35, 0.85)
            colour = (hue.r, hue.g, hue.b, 1.0)
            n_quads = sum(len(i) for i in idx)
            geometry = gfx.Geometry(
                positions=(pos64 - self._origin).astype(np.float32),
                normals=np.concatenate(nrm).astype(np.float32),
                indices=np.concatenate(idx),
                colors=np.tile(np.asarray(colour, np.float32), (n_quads, 1)),
            )
            mesh = gfx.Mesh(geometry, gfx.MeshPhongMaterial(color_mode="face", pick_write=True))
            edges64 = edge_segments(nodes).astype(np.float64)
            edges = gfx.Line(
                gfx.Geometry(positions=(edges64 - self._origin).astype(np.float32)),
                gfx.LineSegmentMaterial(thickness=1.5, color="#202020"),
            )
            block = _Block(key, nodes, mesh, edges, pos64, edges64, colour)
            self._blocks.append(block)
            self._by_mesh[id(mesh)] = block
            self.content.add(mesh, edges)
        self.request_draw()

    def set_face_scalars(
        self, values: dict[Any, np.ndarray] | None, clim: tuple[float, float] = (0.0, 1.0)
    ) -> None:
        for block in self._blocks:
            colours = block.mesh.geometry.colors
            v = None if values is None else values.get(block.key)
            if v is None:
                colours.data[:] = np.asarray(block.colour, np.float32)
            else:
                colours.data[:] = colormap(np.asarray(v, np.float64), clim)
            colours.update_full()
        self.request_draw()

    def face_colours(self, key: Any) -> np.ndarray:
        """The uploaded per-quad RGBA of a block (for tests)."""
        block = next(b for b in self._blocks if b.key == key)
        return np.asarray(block.mesh.geometry.colors.data)

    @property
    def keys(self) -> list[Any]:
        return [b.key for b in self._blocks]

    def uploaded_vertices(self, key: Any) -> int:
        block = next(b for b in self._blocks if b.key == key)
        return len(block.mesh.geometry.positions.data)

    # -- origin (float32 rebase) ------------------------------------------------------

    def _set_origin(self, origin: np.ndarray, *, move_camera: bool = True) -> None:
        delta = np.asarray(origin, np.float64) - self._origin
        self._origin = np.asarray(origin, np.float64)
        for block in self._blocks:
            buf = block.mesh.geometry.positions
            buf.data[:] = (block.pos64 - self._origin).astype(np.float32)
            buf.update_full()
            ebuf = block.edges.geometry.positions
            ebuf.data[:] = (block.edges64 - self._origin).astype(np.float32)
            ebuf.update_full()
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
        if not self._blocks:
            return False
        tgt = np.asarray(self.controller.target, np.float64)
        width = self.view_width()
        if width <= 0.0 or float(np.linalg.norm(tgt)) <= self.rebase_factor * width:
            return False
        self._set_origin(self._origin + tgt)
        return True

    # -- camera --------------------------------------------------------------------

    def show_all(self) -> None:
        if self._blocks:
            pts = np.concatenate([b.pos64 for b in self._blocks])
            centre = 0.5 * (pts.min(axis=0) + pts.max(axis=0))
            self._set_origin(centre, move_camera=False)
        self.camera.depth_range = None
        self.camera.show_object(self.content, view_dir=(-1.0, -1.0, -1.0), up=(0.0, 0.0, 1.0))
        if self._blocks:
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

    def project(self, xyz: np.ndarray) -> tuple[float, float]:
        p = np.append(np.asarray(xyz, np.float64) - self._origin, 1.0)
        clip = np.asarray(self.camera.camera_matrix, np.float64) @ p
        ndc = clip[:3] / clip[3]
        w, h = self.renderer.logical_size
        return (ndc[0] + 1.0) * 0.5 * w, (1.0 - ndc[1]) * 0.5 * h

    # -- picking -------------------------------------------------------------------

    def pick(self, x: float, y: float) -> PickResult | None:
        info = self.renderer.get_pick_info((x, y))
        wo = info.get("world_object")
        block = self._by_mesh.get(id(wo)) if wo is not None else None
        if block is None or "face_index" not in info:
            self._highlight.visible = False
            self.request_draw()
            return None
        shape = tuple(int(n) for n in block.nodes.shape[:3])
        face, cell = quad_cell(shape, int(info["face_index"]))
        xyz = quad_centre(block.nodes, face, cell)
        self._show_highlight(block.nodes, face, cell)
        return PickResult(block.key, face, cell, xyz)

    def _show_highlight(self, nodes: np.ndarray, face: str, cell: tuple[int, int, int]) -> None:
        from shore.gui.viewport import FACES

        f = FACES.index(face)
        grid = face_grids(nodes)[f]
        i, j, k = cell
        a, b = {0: (j, k), 1: (j, k), 2: (i, k), 3: (i, k), 4: (i, j), 5: (i, j)}[f]
        quad = np.array([grid[a, b], grid[a + 1, b], grid[a + 1, b + 1], grid[a, b + 1]])
        loop = np.vstack([quad, quad[:1]]) - self._origin
        self._highlight.geometry.positions.data[:] = loop.astype(np.float32)
        self._highlight.geometry.positions.update_full()
        self._highlight.visible = True
        self.request_draw()

    def on_pick(self, callback: Any) -> None:
        """Call ``callback(PickResult | None)`` on a double click in the view."""

        def handler(event: Any) -> None:
            callback(self.pick(event.x, event.y))

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
