"""Redistributing the nodes of a block edge (ADR 0041).

The spacing law of one block edge is applied, as a distribution of arc-length
fractions, to every grid line of the block along that axis, and through the
seams to the blocks that share a face containing that axis (reversed across
a flipped seam), transitively. Every node stays on its line's old polyline
and every line keeps its end nodes, so faces across the axis do not move.

The law is taken over the picked edge's arc length: tanh / tanh2 cluster by
``beta`` (``cluster_at`` as a fraction), ``pinned`` holds ``first_step`` and
``last_step`` (lengths along the picked edge) around an ``interior`` law.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.face import FACE_ALONG_AXES, FaceBC
from shore.mesh.spacing import Spacing1D

if TYPE_CHECKING:
    from shore.mesh.block import HexBlock

LAWS = ("uniform", "tanh", "tanh2", "pinned")
_COINCIDENT = (FaceBC.SHARED, FaceBC.DIRICHLET)


def edge_axis(edge: str) -> str:
    """The axis an edge runs along: the letter its name does not fix."""
    fixed = {edge[0], edge[5]}
    (axis,) = {"i", "j", "k"} - fixed
    return axis


def fractions(spacing: Spacing1D, nodes: np.ndarray) -> np.ndarray:
    """Arc-length fractions (n,) of *spacing* over the polyline *nodes* (n, 3)."""
    if spacing.law not in LAWS:
        raise ParameterError(f"redistribute: law must be one of {LAWS}, got {spacing.law!r}.")
    length = float(np.linalg.norm(np.diff(nodes, axis=0), axis=1).sum())
    if length <= 0.0:
        raise ParameterError("redistribute: the edge has zero length.")
    t = np.asarray(spacing.build(len(nodes), 0.0, length), dtype=np.float64) / length
    if not (np.all(np.diff(t) > 0.0) and abs(t[0]) < 1e-12 and abs(t[-1] - 1.0) < 1e-12):
        raise ParameterError("redistribute: the law does not give increasing nodes over the edge.")
    t[0], t[-1] = 0.0, 1.0
    return t


def resample(lines: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Lines (..., n, 3) resampled at arc-length fractions *t* along themselves.

    Collapsed lines (zero length: polar axes) are kept; every line keeps its
    two end nodes exactly.
    """
    shape = lines.shape
    n = shape[-2]
    flat = lines.reshape(-1, n, 3)
    seg = np.linalg.norm(np.diff(flat, axis=1), axis=-1)  # (m, n-1)
    s = np.concatenate([np.zeros((len(flat), 1)), np.cumsum(seg, axis=1)], axis=1)
    total = s[:, -1]
    ok = total > 0.0
    out = flat.copy()
    if ok.any():
        f = s[ok] / total[ok, None]  # (m, n) nondecreasing in [0, 1]
        pts = flat[ok]
        m = len(f)
        offset = 2.0 * np.arange(m)[:, None]  # one monotone array for all lines
        idx = np.searchsorted((f + offset).ravel(), (t[None, :] + offset).ravel(), side="right") - 1
        local = np.clip(idx.reshape(m, n) - n * np.arange(m)[:, None], 0, n - 2)
        rows = np.arange(m)[:, None]
        f0, f1 = f[rows, local], f[rows, local + 1]
        span = f1 - f0
        w = np.where(span > 0.0, (t[None, :] - f0) / np.where(span > 0.0, span, 1.0), 0.0)
        new = pts[rows, local] + w[..., None] * (pts[rows, local + 1] - pts[rows, local])
        new[:, 0], new[:, -1] = pts[:, 0], pts[:, -1]
        out[ok] = new
    return out.reshape(shape)


def plan(
    blocks: list[HexBlock], block: str, edge: str, t: np.ndarray
) -> list[tuple[HexBlock, str, np.ndarray]]:
    """Every (block, axis, fractions) the redistribution reaches through seams."""
    from shore.mesh.seams import shared_faces

    by_label = {b.label: b for b in blocks}
    owner = {id(f): b for b in blocks for f in b.faces()}
    shared = shared_faces(blocks)
    start = by_label[block]
    todo = [(start, edge_axis(edge), t)]
    done: dict[tuple[str, str], np.ndarray] = {}
    out = []
    while todo:
        b, axis, tt = todo.pop()
        key = (b.label, axis)
        if key in done:
            if not np.allclose(done[key], tt, rtol=0.0, atol=1e-12):
                raise ParameterError(
                    f"redistribute: {b.label} along {axis} would need two distributions "
                    "(the seams loop back with another orientation)."
                )
            continue
        done[key] = tt
        out.append((b, axis, tt))
        for face in b.faces():
            along = FACE_ALONG_AXES[face.name]
            if axis not in along:
                continue  # a face across the axis: its nodes are line ends, kept
            if face.bc is FaceBC.PERIODIC:
                raise ParameterError(
                    f"redistribute: {b.label}.{face.name} is PERIODIC; redistributing "
                    "across periodic seams is not supported."
                )
            if face.bc not in _COINCIDENT or face.partner is None:
                continue
            if (b.label, face.name) in shared or any(
                (owner[id(face.partner)].label, face.partner.name) == s for s in shared
            ):
                raise ParameterError(
                    f"redistribute: {b.label}.{face.name} is part of a face shared with "
                    "several blocks; redistributing across it is not supported."
                )
            partner_axis, flipped = face.axis_map.along[along.index(axis)]
            partner = owner[id(face.partner)]
            todo.append((partner, partner_axis, (1.0 - tt[::-1]) if flipped else tt))
    return out


def redistribute(
    blocks: list[HexBlock], block: str, edge: str, spacing: Spacing1D
) -> list[tuple[str, str]]:
    """Apply *spacing* along *edge* of *block*, in place; the (block, axis) changed.

    New positions are computed from the old ones for every block before any
    is written, so blocks sharing seam nodes as views are not resampled twice.
    """
    from shore.mesh.edge import _edge_slice

    by_label = {b.label: b for b in blocks}
    if block not in by_label:
        raise ParameterError(f"redistribute: no block {block!r}; blocks: {sorted(by_label)}.")
    names = [e.name for e in by_label[block].edges()]
    if edge not in names:
        raise ParameterError(f"redistribute: {edge!r} is not an edge of {block!r}; edges: {names}.")
    t = fractions(spacing, np.asarray(_edge_slice(by_label[block], edge), dtype=np.float64))
    targets = plan(blocks, block, edge, t)
    new = []
    for b, axis, tt in targets:
        a = "ijk".index(axis)
        new.append(np.moveaxis(resample(np.moveaxis(b.nodes, a, -2), tt), -2, a))
    for (b, _, _), nodes in zip(targets, new, strict=True):
        b.nodes[...] = nodes
    return [(b.label, axis) for b, axis, _ in targets]
