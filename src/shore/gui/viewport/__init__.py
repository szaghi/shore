"""The GUI's renderer interface (ADR 0006, 0037).

Only :mod:`shore.gui.viewport.pygfx_backend` imports pygfx, wgpu or
rendercanvas. This module holds the interface and the structured-grid
arithmetic the backends share: a block is drawn as its six boundary faces
(one quad per boundary cell face, in ``FACES`` order) and its twelve edges,
and a picked quad maps back to its (block, face, i, j, k) in O(1).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Protocol

import numpy as np

FACES = ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi")


@dataclass(frozen=True)
class SliceKey:
    """A constant-index plane of a block: ``axis`` in ``"ijk"``, node ``index``."""

    block: Any
    axis: str
    index: int


@dataclass(frozen=True)
class GeometryPick:
    """A picked geometry entity and the float64 point on it."""

    geometry: str
    entity: Any  # shore.geometry.base.EntityRef
    xyz: np.ndarray


@dataclass(frozen=True)
class SkeletonPick:
    """A picked blocking edge (ADR 0042): its key, the node nearest the pick
    (of ``nodes``, counted along the drawn polyline) and that node's point."""

    key: Any  # the skeleton key the viewport was given
    node: int
    nodes: int
    xyz: np.ndarray


@dataclass(frozen=True)
class MarkerPick:
    """A picked blocking face or block marker (ADR 0044): its key and centre."""

    key: Any  # SkeletonKey(component, block, face) or (component, block, "")
    xyz: np.ndarray


@dataclass(frozen=True)
class PickResult:
    """A picked boundary cell face: its block, face and cell, and its centre."""

    key: Any  # the block key the viewport was given (a SliceKey for a slice plane)
    face: str  # a face name, or "i=5" for a slice plane
    cell: tuple[int, int, int]  # (i, j, k) of the cell behind the face, 0-based
    xyz: np.ndarray  # float64 centre of the picked quad, from the CPU arrays


class Viewport(Protocol):
    """What the GUI needs from a renderer."""

    def set_blocks(self, blocks: Sequence[tuple[Any, np.ndarray]]) -> None:
        """Replace the scene with these ``(key, nodes (ni, nj, nk, 3))`` blocks."""

    def set_face_scalars(
        self, values: dict[Any, np.ndarray] | None, clim: tuple[float, float]
    ) -> None:
        """Colour each boundary quad by a value (``boundary_values`` order), or clear."""

    def show_all(self) -> None: ...

    def look_at_point(self, xyz: np.ndarray, distance: float) -> None: ...

    def pick(self, x: float, y: float) -> PickResult | None: ...

    def project(self, xyz: np.ndarray) -> tuple[float, float]:
        """A world point (float64) to logical pixel coordinates."""

    def request_draw(self) -> None: ...


def face_grids(nodes: np.ndarray) -> list[np.ndarray]:
    """The six boundary node grids ``(na, nb, 3)`` of a block, in ``FACES`` order."""
    return [
        nodes[0, :, :],
        nodes[-1, :, :],
        nodes[:, 0, :],
        nodes[:, -1, :],
        nodes[:, :, 0],
        nodes[:, :, -1],
    ]


def quad_indices(na: int, nb: int, base: int = 0) -> np.ndarray:
    """Quads ``((na-1)*(nb-1), 4)`` of an ``na x nb`` node grid, a-major."""
    p, q = np.meshgrid(np.arange(na - 1), np.arange(nb - 1), indexing="ij")
    v00 = (p * nb + q).ravel()
    return base + np.stack([v00, v00 + nb, v00 + nb + 1, v00 + 1], axis=1).astype(np.int32)


def quad_offsets(shape: tuple[int, int, int]) -> np.ndarray:
    """Index of each face's first quad in a block's quad list, plus the total."""
    ni, nj, nk = shape
    counts = [(nj - 1) * (nk - 1)] * 2 + [(ni - 1) * (nk - 1)] * 2 + [(ni - 1) * (nj - 1)] * 2
    return np.concatenate([[0], np.cumsum(counts)])


def quad_cell(shape: tuple[int, int, int], quad: int) -> tuple[str, tuple[int, int, int]]:
    """The face and the (i, j, k) cell behind quad number *quad* of a block."""
    ni, nj, nk = shape
    offsets = quad_offsets(shape)
    f = int(np.searchsorted(offsets, quad, side="right") - 1)
    nb = (nk - 1) if f < 4 else (nj - 1)  # cells along the face's second axis
    a, b = divmod(int(quad - offsets[f]), nb)
    face = FACES[f]
    cell = {
        "i_lo": (0, a, b),
        "i_hi": (ni - 2, a, b),
        "j_lo": (a, 0, b),
        "j_hi": (a, nj - 2, b),
        "k_lo": (a, b, 0),
        "k_hi": (a, b, nk - 2),
    }[face]
    return face, cell


def quad_centre(nodes: np.ndarray, face: str, cell: tuple[int, int, int]) -> np.ndarray:
    """float64 centre of the boundary quad of *cell* on *face*."""
    f = FACES.index(face)
    grid = face_grids(nodes)[f]
    i, j, k = cell
    a, b = {0: (j, k), 1: (j, k), 2: (i, k), 3: (i, k), 4: (i, j), 5: (i, j)}[f]
    quad = grid[a : a + 2, b : b + 2].reshape(-1, 3)
    return quad.mean(axis=0)


def boundary_values(cell_values: np.ndarray) -> np.ndarray:
    """Per boundary quad, the value of the cell behind it (``(ni-1, nj-1, nk-1)`` input).

    Quads follow ``FACES`` order and, within a face, ``quad_indices`` order.
    """
    v = cell_values
    faces = [v[0, :, :], v[-1, :, :], v[:, 0, :], v[:, -1, :], v[:, :, 0], v[:, :, -1]]
    return np.concatenate([f.ravel() for f in faces])


def edge_segments(nodes: np.ndarray) -> np.ndarray:
    """The 12 block edges as line-segment vertex pairs ``(2m, 3)``."""
    e = [nodes[:, i, k] for i in (0, -1) for k in (0, -1)]
    e += [nodes[i, :, k] for i in (0, -1) for k in (0, -1)]
    e += [nodes[i, j, :] for i in (0, -1) for j in (0, -1)]
    return np.concatenate([np.stack([c[:-1], c[1:]], axis=1).reshape(-1, 3) for c in e])


def colormap(values: np.ndarray, clim: tuple[float, float]) -> np.ndarray:
    """RGBA float32 of *values* on a blue-white-red ramp over *clim* (NaN: grey)."""
    lo, hi = clim
    t = np.clip((values - lo) / (hi - lo if hi > lo else 1.0), 0.0, 1.0)
    blue, white, red = np.array([0.23, 0.30, 0.75]), np.ones(3), np.array([0.71, 0.02, 0.15])
    lower = (t < 0.5)[:, None]
    s = np.where(t < 0.5, 2.0 * t, 2.0 * t - 1.0)[:, None]
    rgb = np.where(lower, red + (white - red) * s, white + (blue - white) * s)
    rgb[~np.isfinite(values)] = 0.6
    return np.concatenate([rgb, np.ones((len(values), 1))], axis=1).astype(np.float32)


def slice_grid(nodes: np.ndarray, axis: str, index: int) -> np.ndarray:
    """The node grid ``(na, nb, 3)`` of the constant-*axis* plane *index*."""
    a = "ijk".index(axis)
    return np.take(nodes, index, axis=a)


def slice_cell(
    shape: tuple[int, int, int], axis: str, index: int, quad: int
) -> tuple[int, int, int]:
    """The cell of a slice quad: on the plane's + side, or its - side on the last plane."""
    a = "ijk".index(axis)
    others = [n for d, n in enumerate(shape) if d != a]
    p, q = divmod(quad, others[1] - 1)
    layer = min(index, shape[a] - 2)
    cell = [p, q]
    cell.insert(a, layer)
    return (cell[0], cell[1], cell[2])


def slice_values(cell_values: np.ndarray, axis: str, index: int) -> np.ndarray:
    """Per slice quad, the value of its cell (``slice_cell``), in quad order."""
    a = "ijk".index(axis)
    layer = min(index, cell_values.shape[a] - 1)
    return np.take(cell_values, layer, axis=a).ravel()


# the 8 corners (HexBlock.vertices() order, index -> (i, j, k) as 0 = lo, -1 = hi)
VERTEX_IJK = ((0, 0, 0), (0, 0, -1), (0, -1, 0), (0, -1, -1),
              (-1, 0, 0), (-1, 0, -1), (-1, -1, 0), (-1, -1, -1))  # fmt: skip
EDGE_NAMES = (
    "j_lo_k_lo", "j_hi_k_lo", "j_lo_k_hi", "j_hi_k_hi",
    "i_lo_k_lo", "i_hi_k_lo", "i_lo_k_hi", "i_hi_k_hi",
    "i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi",
)  # fmt: skip


def edge_nodes(nodes: np.ndarray, name: str) -> np.ndarray:
    """The node polyline of block edge *name* (``HexBlock`` edge naming)."""
    fixed = {"i": slice(None), "j": slice(None), "k": slice(None)}
    for token in (name[:4], name[5:]):  # e.g. "j_lo", "k_hi"
        fixed[token[0]] = 0 if token.endswith("lo") else -1
    return nodes[fixed["i"], fixed["j"], fixed["k"]]


def face_edges(face: str) -> tuple[str, ...]:
    """The four edges bounding *face*."""
    return tuple(e for e in EDGE_NAMES if face in (e[:4], e[5:]))


def face_vertices(face: str) -> tuple[int, ...]:
    """The four corners (vertex indices) of *face*."""
    a = "ijk".index(face[0])
    side = 0 if face.endswith("lo") else -1
    return tuple(v for v, ijk in enumerate(VERTEX_IJK) if ijk[a] == side)


def _segment_distance(points: np.ndarray, xyz: np.ndarray) -> float:
    a, b = points[:-1], points[1:]
    ab = b - a
    t = np.clip(
        np.einsum("ij,ij->i", xyz - a, ab) / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300), 0, 1
    )
    return float(np.linalg.norm(a + t[:, None] * ab - xyz, axis=1).min())


def nearest_edge(nodes: np.ndarray, face: str, xyz: np.ndarray) -> str:
    """The edge of *face* closest to the point *xyz* (an edge pick)."""
    return min(face_edges(face), key=lambda e: _segment_distance(edge_nodes(nodes, e), xyz))


def nearest_vertex(nodes: np.ndarray, face: str, xyz: np.ndarray) -> int:
    """The corner of *face* closest to the point *xyz* (a vertex pick)."""
    return min(
        face_vertices(face),
        key=lambda v: float(np.linalg.norm(nodes[VERTEX_IJK[v]] - xyz)),
    )
