"""Grid lines leaving the walls along their normal (ADR 0048).

After a premesh (ADR 0042) the lines across a block run from the wall to the
opposite face along the transfinite interpolation, so they leave a curved
wall obliquely. ``orthogonal_walls`` reshapes every line of each block with
exactly one WALL face into a cubic Hermite curve that leaves the wall node
along the wall's normal (with tangent ``strength`` times the line's chord)
and meets the opposite node with the line's old end direction; the nodes
keep their arc-length fractions, the wall and the opposite face stay put.

The normal of a wall node is the mean of the grid normals of every wall face
holding it (each oriented into its block), keyed by the node's exact
coordinates: blocks sharing a wall node use the same normal, and blocks
sharing a face build its lines identically, so seams stay closed.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from shore.exceptions import ParameterError

if TYPE_CHECKING:
    from shore.mesh.block import HexBlock


def _wall(block: HexBlock, name: str) -> tuple[int, int]:
    """The axis of face *name* and its node index (0 or -1)."""
    return "ijk".index(name[0]), (0 if name.endswith("lo") else -1)


def grid_normals(surface: np.ndarray) -> np.ndarray:
    """Unit normals (na, nb, 3) of a structured surface grid (na, nb, 3)."""
    n = np.cross(np.gradient(surface, axis=0), np.gradient(surface, axis=1))
    return n / np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-300)


def wall_normals(blocks: list[HexBlock]) -> dict[bytes, np.ndarray]:
    """Per wall node, by its coordinates' bytes: the mean grid normal of the
    wall faces holding it, each oriented into its block."""
    acc: dict[bytes, list[np.ndarray]] = {}
    for b in blocks:
        for f in b.faces():
            if f.bc.name != "WALL":
                continue
            axis, side = _wall(b, f.name)
            wall = np.take(b.nodes, side, axis=axis)
            inner = np.take(b.nodes, 1 if side == 0 else -2, axis=axis)
            n = grid_normals(wall)
            if np.einsum("...i,...i", n, inner - wall).sum() < 0.0:
                n = -n
            for p, v in zip(wall.reshape(-1, 3), n.reshape(-1, 3), strict=True):
                acc.setdefault(p.tobytes(), []).append(v)
    out = {}
    for key, vs in acc.items():
        m = np.mean(vs, axis=0)
        out[key] = m / max(float(np.linalg.norm(m)), 1e-300)
    return out


def hermite_lines(
    lines: np.ndarray, normals: np.ndarray, strength: float, samples: int = 64
) -> np.ndarray:
    """Lines (m, n, 3), wall node first, reshaped to leave along *normals* (m, 3).

    Each line becomes the cubic Hermite curve from its wall node (tangent
    ``strength * chord * normal``) to its last node (tangent ``chord`` times
    the old last segment's direction); its nodes keep their arc-length
    fractions, measured on the old line and placed on the curve.
    """
    w, o = lines[:, 0], lines[:, -1]
    chord = np.linalg.norm(o - w, axis=-1)[:, None]
    end = o - lines[:, -2]
    t0 = strength * chord * normals
    t1 = chord * end / np.maximum(np.linalg.norm(end, axis=-1, keepdims=True), 1e-300)
    s = np.linspace(0.0, 1.0, samples)[None, :, None]
    curve = (
        (2 * s**3 - 3 * s**2 + 1) * w[:, None]
        + (s**3 - 2 * s**2 + s) * t0[:, None]
        + (-2 * s**3 + 3 * s**2) * o[:, None]
        + (s**3 - s**2) * t1[:, None]
    )
    old = np.concatenate(
        [np.zeros((len(lines), 1)), np.cumsum(np.linalg.norm(np.diff(lines, axis=1), axis=-1), 1)],
        axis=1,
    )
    old = old / np.maximum(old[:, -1:], 1e-300)
    arc = np.concatenate(
        [np.zeros((len(curve), 1)), np.cumsum(np.linalg.norm(np.diff(curve, axis=1), axis=-1), 1)],
        axis=1,
    )
    arc = arc / np.maximum(arc[:, -1:], 1e-300)
    out = np.empty_like(lines)
    for m in range(len(lines)):
        for d in range(3):
            out[m, :, d] = np.interp(old[m], arc[m], curve[m, :, d])
    out[:, 0], out[:, -1] = w, o
    return out


def orthogonal_walls(blocks: list[HexBlock], strength: float = 1.0) -> list[str]:
    """Reshape, in place, the lines of every block with exactly one WALL face
    to leave the wall along its normal; the labels of the blocks changed."""
    if not strength > 0.0:
        raise ParameterError(f"orthogonalize: strength must be > 0, got {strength}.")
    normals = wall_normals(blocks)
    changed = []
    for b in blocks:
        walls = [f.name for f in b.faces() if f.bc.name == "WALL"]
        if len(walls) != 1:
            continue
        axis, side = _wall(b, walls[0])
        x = np.moveaxis(b.nodes, axis, -2)
        if side == -1:
            x = x[..., ::-1, :]
        shape = x.shape
        lines = x.reshape(-1, shape[-2], 3)
        nrm = np.array([normals[p.tobytes()] for p in lines[:, 0]])
        new = hermite_lines(lines, nrm, strength).reshape(shape)
        if side == -1:
            new = new[..., ::-1, :]
        b.nodes[...] = np.moveaxis(new, -2, axis)
        changed.append(b.label)
    return changed
