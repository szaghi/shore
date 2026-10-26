"""Viewport interface and its pygfx backend (ADR 0006 spike).

Only this module imports pygfx/wgpu. The rest of the spike talks to the
``Viewport`` protocol, which is the isolation boundary ADR 0006 prescribes.

Origin modes (the float32 trap test):

- ``absolute``: vertex positions uploaded as absolute float32 world coords.
- ``block``: positions relative to each block centroid; the centroid goes in
  the object's world transform. pygfx applies ``world_transform`` and then
  ``cam_transform`` on the GPU in float32, so the large translation still
  cancels in float32. Expected to fail the trap.
- ``floating``: like ``block``, but every object translation is taken
  relative to a movable world origin, computed in float64 on the CPU. The
  origin follows the camera target, so GPU-side numbers stay small.
- ``tiled``: like ``floating``, but each block face is split by a quadtree
  into tiles with their own local origin, until every tile satisfies
  ``min cell edge >= safety * eps32 * max|p - tile origin|``.
- ``rebase``: one mesh per block, but on every origin move all positions are
  re-uploaded as ``float32(p - origin)`` from float64 on the CPU, so the
  vertices near the camera target carry full float32 precision.

Rebase trigger (all floating-origin modes): before each draw, if the orbit
target has drifted from the origin by more than ``rebase_factor`` view widths
``V``, the origin moves to the target. After a rebase, visible vertices
(within ~V of the target) are in error by ~``eps32 * (drift + V)``; for that
to stay under 0.1 px at ~1000 px, ``drift <~ 1e-4 V / eps32 ~ 840 V``. The
default factor 100 keeps a ~8x margin; a rebase costs ~4 ms for 1e7 nodes.

Throwaway spike code: not part of the ``shore`` package.
"""

from __future__ import annotations

import time
from collections import deque
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal, Protocol

import numpy as np
import pygfx as gfx
from synthetic import Block

OriginMode = Literal["absolute", "block", "floating", "tiled", "rebase"]
EPS32 = float(np.finfo(np.float32).eps)
FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")


@dataclass
class PickResult:
    """Structured-index answer to a screen pick."""

    block: str
    face: str
    cell: tuple[int, int, int]
    xyz: np.ndarray  # float64 centre of the picked boundary quad
    elapsed_ms: float


@dataclass
class FrameStats:
    """Frame interval statistics in milliseconds."""

    n: int
    mean_ms: float
    p95_ms: float
    max_ms: float


class Viewport(Protocol):
    """What the GUI needs from a renderer."""

    def add_block(self, block: Block) -> None: ...
    def clear(self) -> None: ...
    def show_all(self) -> None: ...
    def look_at_point(
        self, xyz: np.ndarray, distance: float, view_dir: np.ndarray | None = None
    ) -> None: ...
    def pick(self, x: float, y: float) -> PickResult | None: ...
    def project(self, xyz: np.ndarray) -> tuple[float, float]: ...
    def frame_stats(self) -> FrameStats: ...
    def request_draw(self) -> None: ...


def _face_grids(nodes: np.ndarray) -> list[np.ndarray]:
    """The six boundary node grids ``(na, nb, 3)`` of a block, in FACES order."""
    return [
        nodes[0, :, :],
        nodes[-1, :, :],
        nodes[:, 0, :],
        nodes[:, -1, :],
        nodes[:, :, 0],
        nodes[:, :, -1],
    ]


def _quad_indices(na: int, nb: int, base: int) -> np.ndarray:
    """Quad connectivity ``(na-1)*(nb-1), 4`` of an ``na x nb`` node grid."""
    p, q = np.meshgrid(np.arange(na - 1), np.arange(nb - 1), indexing="ij")
    v00 = (p * nb + q).ravel()
    return base + np.stack([v00, v00 + nb, v00 + nb + 1, v00 + 1], axis=1).astype(np.int32)


def _grid_normals(grid: np.ndarray) -> np.ndarray:
    """Unit vertex normals of a node grid from central differences."""
    da = np.gradient(grid, axis=0)
    db = np.gradient(grid, axis=1)
    n = np.cross(da, db)
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    return n / np.where(norm > 0.0, norm, 1.0)


def _edge_segments(nodes: np.ndarray) -> np.ndarray:
    """The 12 block edges as line-segment vertex pairs ``(2m, 3)``."""
    e = [nodes[:, i, k] for i in (0, -1) for k in (0, -1)]
    e += [nodes[i, :, k] for i in (0, -1) for k in (0, -1)]
    e += [nodes[i, j, :] for i in (0, -1) for j in (0, -1)]
    return np.concatenate([np.stack([c[:-1], c[1:]], axis=1).reshape(-1, 3) for c in e])


def _tiles(
    g: np.ndarray, safety: float, max_tiles: int
) -> tuple[list[tuple[int, int, int, int]], bool]:
    """Quadtree tiles ``(a0, a1, b0, b1)`` (inclusive node ranges) of a face grid.

    A tile is accepted when its smallest cell edge is at least
    ``safety * eps32 * max|p - tile centroid|`` (component-wise, as float32
    rounds per component), or when it is a single cell. Returns the tiles and
    whether ``max_tiles`` capped the refinement.
    """
    na, nb, _ = g.shape
    da = np.linalg.norm(np.diff(g, axis=0), axis=-1)  # (na-1, nb)
    db = np.linalg.norm(np.diff(g, axis=1), axis=-1)  # (na, nb-1)
    h = np.minimum(np.minimum(da[:, :-1], da[:, 1:]), np.minimum(db[:-1], db[1:]))
    stack = [(0, na - 1, 0, nb - 1)]
    out: list[tuple[int, int, int, int]] = []
    capped = False
    while stack:
        a0, a1, b0, b1 = stack.pop()
        sub = g[a0 : a1 + 1, b0 : b1 + 1].reshape(-1, 3)
        r = float(np.abs(sub - sub.mean(axis=0)).max())
        single = a1 - a0 == 1 and b1 - b0 == 1
        if single or h[a0:a1, b0:b1].min() >= safety * EPS32 * r:
            out.append((a0, a1, b0, b1))
            continue
        if len(out) + len(stack) + 2 > max_tiles:
            capped = True
            out.append((a0, a1, b0, b1))
            continue
        am, bm = (a0 + a1) // 2, (b0 + b1) // 2
        la = float(np.linalg.norm(g[a1, bm] - g[a0, bm]))
        lb = float(np.linalg.norm(g[am, b1] - g[am, b0]))
        if b1 - b0 == 1 or (a1 - a0 > 1 and la >= lb):
            stack += [(a0, am, b0, b1), (am, a1, b0, b1)]
        else:
            stack += [(a0, a1, b0, bm), (a0, a1, bm, b1)]
    return out, capped


@dataclass
class _Piece:
    """One pygfx mesh: a whole block surface (6 segments) or one face tile."""

    block_index: int
    mesh: gfx.Mesh
    origin: np.ndarray  # float64 origin the uploaded float32 positions are relative to
    pos64: np.ndarray | None  # float64 absolute positions, kept in ``rebase`` mode
    quad_offsets: np.ndarray  # (n_segments + 1,)
    seg_face: list[int]
    seg_a0: list[int]
    seg_b0: list[int]
    seg_nb: list[int]  # cells along b in the segment


@dataclass
class _BlockGfx:
    block: Block
    centroid: np.ndarray
    edges: gfx.Line
    edges_origin: np.ndarray
    edges_pos64: np.ndarray
    pieces: list[_Piece] = field(default_factory=list)


class PygfxViewport:
    """pygfx/wgpu implementation of ``Viewport``."""

    def __init__(
        self,
        canvas: object,
        origin_mode: OriginMode = "floating",
        tile_safety: float = 4.0,
        max_tiles_per_face: int = 4096,
        rebase_factor: float = 100.0,
    ) -> None:
        self.canvas = canvas
        self.origin_mode: OriginMode = origin_mode
        self.tile_safety = tile_safety
        self.max_tiles_per_face = max_tiles_per_face
        self.tiles_capped = False
        self.last_rebase_ms = 0.0
        self.rebase_factor = rebase_factor
        self.rebase_count = 0
        self.on_rebase: Callable[[dict], None] | None = None
        self.renderer = gfx.WgpuRenderer(canvas)
        self.scene = gfx.Scene()
        self.scene.add(gfx.AmbientLight(intensity=0.4))
        self.camera = gfx.PerspectiveCamera(45, 16 / 9)
        self.camera.add(gfx.DirectionalLight(intensity=2.5))
        self.scene.add(self.camera)
        self.content = gfx.Group()  # block objects only: what "fit" frames
        self.scene.add(self.content)
        self.controller = gfx.OrbitController(self.camera, register_events=self.renderer)
        self._blocks: list[_BlockGfx] = []
        self._by_mesh: dict[int, _Piece] = {}
        self._origin = np.zeros(3)
        self.target = np.zeros(3)  # orbit centre, scene coordinates
        self._intervals: deque[float] = deque(maxlen=600)
        self._last = 0.0
        self._highlight = gfx.Line(
            gfx.Geometry(positions=np.zeros((5, 3), np.float32)),
            gfx.LineMaterial(thickness=3.0, color="#ff3b30", depth_test=False),
        )
        self._highlight.visible = False
        self.scene.add(self._highlight)
        self._background = gfx.Background.from_color("#f4f5f7", "#d9dce1")
        self.scene.add(self._background)
        self.before_draw: Callable[[], None] | None = None
        canvas.request_draw(self._draw)

    def set_background(self, top: str, bottom: str) -> None:
        """Vertical gradient background (theme hook)."""
        self.scene.remove(self._background)
        self._background = gfx.Background.from_color(top, bottom)
        self.scene.add(self._background)
        self.request_draw()

    # -- scene content -------------------------------------------------------

    def _make_piece(
        self,
        index: int,
        grids: list[np.ndarray],
        normals: list[np.ndarray],
        segs: list[tuple[int, int, int, int, int]],
        origin: np.ndarray,
        color: gfx.Color,
    ) -> _Piece:
        """One mesh from face sub-grids ``(face, a0, a1, b0, b1)``."""
        pos, nrm, idx, offsets = [], [], [], [0]
        base = 0
        for f, a0, a1, b0, b1 in segs:
            sub = grids[f][a0 : a1 + 1, b0 : b1 + 1]
            na, nb, _ = sub.shape
            pos.append(sub.reshape(-1, 3))
            nrm.append(normals[f][a0 : a1 + 1, b0 : b1 + 1].reshape(-1, 3))
            idx.append(_quad_indices(na, nb, base))
            offsets.append(offsets[-1] + (na - 1) * (nb - 1))
            base += na * nb
        pos64 = np.concatenate(pos)
        geometry = gfx.Geometry(
            positions=(pos64 - origin).astype(np.float32),
            normals=np.concatenate(nrm).astype(np.float32),
            indices=np.concatenate(idx),
        )
        mesh = gfx.Mesh(geometry, gfx.MeshPhongMaterial(color=color, pick_write=True))
        piece = _Piece(
            block_index=index,
            mesh=mesh,
            origin=np.array(origin, dtype=np.float64),
            pos64=pos64 if self.origin_mode == "rebase" else None,
            quad_offsets=np.asarray(offsets),
            seg_face=[s[0] for s in segs],
            seg_a0=[s[1] for s in segs],
            seg_b0=[s[3] for s in segs],
            seg_nb=[s[4] - s[3] for s in segs],
        )
        self._by_mesh[id(mesh)] = piece
        return piece

    def add_block(self, block: Block) -> None:
        nodes = block.nodes
        index = len(self._blocks)
        mode = self.origin_mode
        centroid = nodes.reshape(-1, 3).mean(axis=0)
        grids = _face_grids(nodes)
        normals = [_grid_normals(g) for g in grids]
        color = gfx.Color.from_hsv((index * 0.137) % 1.0, 0.35, 0.85)
        if mode == "absolute":
            origin = np.zeros(3)
        elif mode == "rebase":
            origin = self._origin.copy()
        else:
            origin = centroid
        edges_pos64 = _edge_segments(nodes)
        edges = gfx.Line(
            gfx.Geometry(positions=(edges_pos64 - origin).astype(np.float32)),
            gfx.LineSegmentMaterial(thickness=1.5, color="#202020"),
        )
        bg = _BlockGfx(block, centroid, edges, np.array(origin), edges_pos64)
        if mode == "tiled":
            for f, g in enumerate(grids):
                tiles, capped = _tiles(g, self.tile_safety, self.max_tiles_per_face)
                self.tiles_capped |= capped
                for a0, a1, b0, b1 in tiles:
                    tile_origin = g[a0 : a1 + 1, b0 : b1 + 1].reshape(-1, 3).mean(axis=0)
                    segs = [(f, a0, a1, b0, b1)]
                    bg.pieces.append(
                        self._make_piece(index, grids, normals, segs, tile_origin, color)
                    )
        else:
            segs = [(f, 0, g.shape[0] - 1, 0, g.shape[1] - 1) for f, g in enumerate(grids)]
            bg.pieces.append(self._make_piece(index, grids, normals, segs, origin, color))
        self._blocks.append(bg)
        self.content.add(*(p.mesh for p in bg.pieces), edges)
        self._place(bg)

    def clear(self) -> None:
        for bg in self._blocks:
            self.content.remove(*(p.mesh for p in bg.pieces), bg.edges)
        self._blocks.clear()
        self._by_mesh.clear()
        self._highlight.visible = False

    @property
    def n_nodes(self) -> int:
        return sum(bg.block.n_nodes for bg in self._blocks)

    @property
    def n_pieces(self) -> int:
        """Number of pygfx mesh objects (tiles) in the scene."""
        return sum(len(bg.pieces) for bg in self._blocks)

    # -- origin handling (float32 trap) ------------------------------------------

    def _place(self, bg: _BlockGfx) -> None:
        for p in bg.pieces:
            p.mesh.local.position = p.origin - self._origin
        bg.edges.local.position = bg.edges_origin - self._origin

    def _rebase(self) -> None:
        """Re-upload every position as float32(p - origin) (``rebase`` mode)."""
        t0 = time.perf_counter()
        for bg in self._blocks:
            for p in bg.pieces:
                buf = p.mesh.geometry.positions
                buf.data[:] = (p.pos64 - self._origin).astype(np.float32)
                buf.update_full()
                p.origin = self._origin.copy()
            buf = bg.edges.geometry.positions
            buf.data[:] = (bg.edges_pos64 - self._origin).astype(np.float32)
            buf.update_full()
            bg.edges_origin = self._origin.copy()
        self.last_rebase_ms = 1e3 * (time.perf_counter() - t0)

    def set_origin(self, origin: np.ndarray) -> None:
        """Move the world origin (float64); camera keeps its world pose."""
        if self.origin_mode in ("absolute", "block"):
            return
        delta = np.asarray(origin, dtype=np.float64) - self._origin
        self._origin = np.asarray(origin, dtype=np.float64)
        if self.origin_mode == "rebase":
            self._rebase()
        for bg in self._blocks:
            self._place(bg)
        self.camera.local.position = np.asarray(self.camera.local.position) - delta
        self._set_target(self.target - delta)

    def view_width(self) -> float:
        """Width of the view at the orbit target, in world units."""
        tgt = np.asarray(self.controller.target, dtype=np.float64)
        dist = float(np.linalg.norm(np.asarray(self.camera.local.position) - tgt))
        fov = float(self.camera.fov)
        if fov > 0.0:
            return 2.0 * dist * np.tan(np.radians(fov) / 2.0)
        return float(self.camera.width)

    def maybe_rebase(self) -> bool:
        """Apply the rebase trigger; returns True if the origin moved."""
        if self.origin_mode in ("absolute", "block") or not self._blocks:
            return False
        tgt = np.asarray(self.controller.target, dtype=np.float64)  # scene coords
        drift = float(np.linalg.norm(tgt))
        width = self.view_width()
        if width <= 0.0 or drift <= self.rebase_factor * width:
            return False
        self.target = tgt
        self.set_origin(self._origin + tgt)
        self.rebase_count += 1
        if self.on_rebase is not None:
            self.on_rebase(
                {
                    "origin": self._origin.copy(),
                    "drift_views": drift / width,
                    "view_width": width,
                    "ms": self.last_rebase_ms,
                }
            )
        return True

    def to_scene(self, xyz: np.ndarray) -> np.ndarray:
        """World float64 coordinates to scene (GPU-side) coordinates."""
        return np.asarray(xyz, dtype=np.float64) - self._origin

    # -- camera ------------------------------------------------------------------

    def show_all(self) -> None:
        if self._blocks:
            centre = np.mean([bg.centroid for bg in self._blocks], axis=0)
            self.set_origin(centre)
            self._set_target(self.to_scene(centre))
        self.camera.depth_range = None
        self.camera.show_object(self.content, view_dir=(-1.0, -1.0, -1.0), up=(0.0, 0.0, 1.0))
        self.request_draw()

    def look_at_point(
        self, xyz: np.ndarray, distance: float, view_dir: np.ndarray | None = None
    ) -> None:
        """Place the camera ``distance`` from ``xyz``, looking along ``view_dir``."""
        self.set_origin(xyz)
        target = self.to_scene(xyz)
        d = np.array([-0.3, -0.4, -0.866]) if view_dir is None else np.asarray(view_dir, float)
        d = d / np.linalg.norm(d)
        eye = target - distance * d
        # look_at is degenerate when the view direction is parallel to "up"
        up = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
        self.camera.world.reference_up = up
        self.camera.local.position = eye
        self.camera.look_at(target)
        self.camera.depth_range = (distance * 1e-3, distance * 1e3)
        self._set_target(target)
        self.request_draw()

    def _set_target(self, target: np.ndarray) -> None:
        """Keep our orbit centre and the controller's in sync."""
        self.target = np.asarray(target, dtype=np.float64)
        self.controller.target = tuple(self.target)

    def project(self, xyz: np.ndarray) -> tuple[float, float]:
        """World point to logical pixel coordinates (float64 on the CPU)."""
        p = np.append(self.to_scene(xyz), 1.0)
        clip = np.asarray(self.camera.camera_matrix, dtype=np.float64) @ p
        ndc = clip[:3] / clip[3]
        w, h = self.renderer.logical_size
        return (ndc[0] + 1.0) * 0.5 * w, (1.0 - ndc[1]) * 0.5 * h

    # -- picking -----------------------------------------------------------------

    def pick(self, x: float, y: float) -> PickResult | None:
        t0 = time.perf_counter()
        info = self.renderer.get_pick_info((x, y))
        wo = info.get("world_object")
        piece = self._by_mesh.get(id(wo)) if wo is not None else None
        if piece is None or "face_index" not in info:
            self._highlight.visible = False
            self.request_draw()
            return None
        q = int(info["face_index"])
        s = int(np.searchsorted(piece.quad_offsets, q, side="right") - 1)
        a, b = divmod(q - int(piece.quad_offsets[s]), piece.seg_nb[s])
        a += piece.seg_a0[s]
        b += piece.seg_b0[s]
        f = piece.seg_face[s]
        block = self._blocks[piece.block_index].block
        ni, nj, nk, _ = block.nodes.shape
        cell = {
            "i_lo": (0, a, b),
            "i_hi": (ni - 2, a, b),
            "j_lo": (a, 0, b),
            "j_hi": (a, nj - 2, b),
            "k_lo": (a, b, 0),
            "k_hi": (a, b, nk - 2),
        }[FACES[f]]
        grid = _face_grids(block.nodes)[f]
        quad = np.array([grid[a, b], grid[a + 1, b], grid[a + 1, b + 1], grid[a, b + 1]])
        loop = np.vstack([quad, quad[:1]]) - self._origin
        self._highlight.geometry.positions.data[:] = loop.astype(np.float32)
        self._highlight.geometry.positions.update_full()
        self._highlight.visible = True
        self.request_draw()
        return PickResult(
            block=block.name,
            face=FACES[f],
            cell=cell,
            xyz=quad.mean(axis=0),
            elapsed_ms=1e3 * (time.perf_counter() - t0),
        )

    # -- drawing and timing ------------------------------------------------------

    def _draw(self) -> None:
        now = time.perf_counter()
        if self._last:
            self._intervals.append(1e3 * (now - self._last))
        self._last = now
        self.maybe_rebase()
        if self.before_draw is not None:
            self.before_draw()
        self.renderer.render(self.scene, self.camera)

    def reset_stats(self) -> None:
        self._intervals.clear()
        self._last = 0.0

    def frame_stats(self) -> FrameStats:
        a = np.asarray(self._intervals)
        if a.size == 0:
            return FrameStats(0, 0.0, 0.0, 0.0)
        return FrameStats(a.size, float(a.mean()), float(np.percentile(a, 95)), float(a.max()))

    def request_draw(self) -> None:
        self.canvas.request_draw()
