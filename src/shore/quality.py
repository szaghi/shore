"""Structured-grid quality metrics (roadmap phase 4).

All functions take one block's node array ``(ni, nj, nk, 3)`` and follow the
marcher's orientation convention: a valid cell has a positive Jacobian
``(e_i x e_j) . e_k`` with ``k`` pointing away from the wall.

- ``scaled_jacobian``: per cell, the minimum over its 8 corners of
  ``det[e_i, e_j, e_k] / (|e_i| |e_j| |e_k|)`` with the edges leaving that
  corner (the Verdict / CUBIT definition). 1 for a rectangular cell, the
  cosine of the shear angle for a sheared one, <= 0 for a folded corner.
- ``k_orthogonality_deg``: per node of layer ``k``, the angle between the
  k-segment to layer ``k + 1`` and the layer's normal ``t_i x t_j``
  (central-difference tangents). 0 is orthogonal; > 90 means the segment
  points back through the layer.
- ``k_angle_deg``: the same angle against given directions (e.g. the body's
  surface normals), unbiased at block-edge nodes.
- ``k_spacing`` / ``k_growth``: k-segment lengths and their ratios.
- ``aspect_ratio``: per cell, longest over shortest of its 12 edges.
- ``summary``: the extremes of all of these for a list of blocks.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import numpy as np


def _corner(nodes: np.ndarray, a: int, b: int, c: int) -> np.ndarray:
    """Corner ``(a, b, c)`` (each 0 or 1) of every cell, (ni-1, nj-1, nk-1, 3)."""
    ni, nj, nk, _ = nodes.shape
    return nodes[a : ni - 1 + a, b : nj - 1 + b, c : nk - 1 + c]


def scaled_jacobian(nodes: np.ndarray) -> np.ndarray:
    """Minimum corner scaled Jacobian per cell, shape (ni-1, nj-1, nk-1)."""
    out = np.full(tuple(s - 1 for s in nodes.shape[:3]), np.inf)
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                p = _corner(nodes, a, b, c)
                # edges leaving this corner, oriented along +i, +j, +k
                ei = (_corner(nodes, 1 - a, b, c) - p) * (1 - 2 * a)
                ej = (_corner(nodes, a, 1 - b, c) - p) * (1 - 2 * b)
                ek = (_corner(nodes, a, b, 1 - c) - p) * (1 - 2 * c)
                det = np.einsum("...j,...j->...", np.cross(ei, ej), ek)
                lengths = (
                    np.linalg.norm(ei, axis=-1)
                    * np.linalg.norm(ej, axis=-1)
                    * np.linalg.norm(ek, axis=-1)
                )
                with np.errstate(divide="ignore", invalid="ignore"):
                    sj = np.where(lengths > 0.0, det / lengths, -1.0)
                out = np.minimum(out, sj)
    return out


def _layer_normals(layer: np.ndarray) -> np.ndarray:
    """Unit normals ``t_i x t_j`` of a layer (ni, nj, 3), central differences."""
    ti = np.gradient(layer, axis=0)
    tj = np.gradient(layer, axis=1)
    n = np.cross(ti, tj)
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    return n / np.where(norm > 0.0, norm, 1.0)


def k_orthogonality_deg(nodes: np.ndarray, k: int = 0) -> np.ndarray:
    """Angle (deg) between the k-segment ``k -> k+1`` and layer ``k``'s normal, (ni, nj)."""
    seg = nodes[:, :, k + 1] - nodes[:, :, k]
    n = _layer_normals(nodes[:, :, k])
    length = np.linalg.norm(seg, axis=-1)
    cos = np.einsum("...j,...j->...", seg, n) / np.where(length > 0.0, length, 1.0)
    return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))


def k_angle_deg(nodes: np.ndarray, directions: np.ndarray, k: int = 0) -> np.ndarray:
    """Angle (deg) between the k-segment ``k -> k+1`` and given *directions* (ni, nj, 3).

    With the body's outward surface normals at layer 0 this is the wall
    orthogonality against the true surface: unlike ``k_orthogonality_deg``
    it is not biased at block-edge nodes, where a one-sided layer normal is
    not the surface normal across the seam.
    """
    seg = nodes[:, :, k + 1] - nodes[:, :, k]
    d = np.asarray(directions, dtype=np.float64)
    d = d / np.where(
        np.linalg.norm(d, axis=-1, keepdims=True) > 0.0,
        np.linalg.norm(d, axis=-1, keepdims=True),
        1.0,
    )
    length = np.linalg.norm(seg, axis=-1)
    cos = np.einsum("...j,...j->...", seg, d) / np.where(length > 0.0, length, 1.0)
    return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))


def k_spacing(nodes: np.ndarray) -> np.ndarray:
    """Length of every k-segment, shape (ni, nj, nk-1)."""
    return np.linalg.norm(np.diff(nodes, axis=2), axis=-1)


def k_growth(nodes: np.ndarray) -> np.ndarray:
    """Ratio of consecutive k-segment lengths, shape (ni, nj, nk-2)."""
    s = k_spacing(nodes)
    return s[:, :, 1:] / np.where(s[:, :, :-1] > 0.0, s[:, :, :-1], np.nan)


def aspect_ratio(nodes: np.ndarray) -> np.ndarray:
    """Longest over shortest of each cell's 12 edges, shape (ni-1, nj-1, nk-1)."""
    edges = []
    for a in (0, 1):
        for b in (0, 1):
            edges.append(_corner(nodes, 1, a, b) - _corner(nodes, 0, a, b))  # along i
            edges.append(_corner(nodes, a, 1, b) - _corner(nodes, a, 0, b))  # along j
            edges.append(_corner(nodes, a, b, 1) - _corner(nodes, a, b, 0))  # along k
    lengths = np.stack([np.linalg.norm(e, axis=-1) for e in edges])
    shortest = lengths.min(axis=0)
    with np.errstate(divide="ignore"):
        return np.where(shortest > 0.0, lengths.max(axis=0) / shortest, np.inf)


def summary(blocks: Iterable[Any]) -> dict[str, float]:
    """Extremes over *blocks* (objects with ``.nodes``, or node arrays)."""
    arrays = [getattr(b, "nodes", b) for b in blocks]
    sj = np.concatenate([scaled_jacobian(a).ravel() for a in arrays])
    wall = np.concatenate([k_orthogonality_deg(a, 0).ravel() for a in arrays])
    ar = np.concatenate([aspect_ratio(a).ravel() for a in arrays])
    growth = [k_growth(a).ravel() for a in arrays if a.shape[2] > 2]
    g = np.concatenate(growth) if growth else np.array([np.nan])
    return {
        "scaled_jacobian_min": float(sj.min()),
        "folded_cells": int((sj <= 0.0).sum()),
        "wall_orthogonality_max_deg": float(wall.max()),
        "wall_orthogonality_mean_deg": float(wall.mean()),
        "aspect_ratio_max": float(ar.max()),
        "k_growth_max": float(np.nanmax(g)),
    }
