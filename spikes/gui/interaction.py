"""Simulated user zoom for testing the rebase trigger (ADR 0006 spike).

The zoom touches only what a mouse interaction touches: the camera pose and
the orbit controller's (float32) target, in scene coordinates. It never calls
``set_origin`` itself, so any rebase must come from the viewport's trigger.

Throwaway spike code: not part of the ``shore`` package.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np
from synthetic import Block
from viewport import PygfxViewport


@dataclass
class WallCell:
    """A boundary-layer cell on a block's ``i_lo`` face."""

    block: str
    centre: np.ndarray  # float64 world
    normal: np.ndarray  # face normal (view direction)
    height: float  # wall-normal size


def wall_cell(block: Block, k: int = 0) -> WallCell:
    """The ``i_lo`` cell at the polar cluster centre (see ``aspect_ratio``)."""
    nodes = block.nodes
    b = (nodes.shape[1] - 1) // 2
    quad = np.array(
        [nodes[0, b, k], nodes[0, b + 1, k], nodes[0, b + 1, k + 1], nodes[0, b, k + 1]]
    )
    normal = np.cross(quad[1] - quad[0], quad[3] - quad[0])
    return WallCell(
        block=block.name,
        centre=quad.mean(axis=0),
        normal=normal / np.linalg.norm(normal),
        height=float(np.linalg.norm(quad[3] - quad[0])),
    )


def final_distance(cell: WallCell, height_px: float = 720.0, cell_px: float = 20.0) -> float:
    """Camera distance at which the cell is ``cell_px`` tall (fov 45 deg)."""
    return cell.height * height_px / (cell_px * 2.0 * np.tan(np.radians(22.5)))


def zoom_to(vp: PygfxViewport, cell: WallCell, d_start: float, n_steps: int) -> Iterator[float]:
    """Geometric zoom onto ``cell``; yields the distance after each step.

    The caller draws after each yield (the trigger runs inside the draw).
    The pivot is set once, like a user "set pivot" click, through the
    controller's float32 target; every later step reads the target back from
    the controller, so a rebase in between is followed transparently.
    """
    vp.controller.target = tuple(vp.to_scene(cell.centre))
    d = cell.normal
    up = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    vp.camera.world.reference_up = up
    for dist in np.geomspace(d_start, final_distance(cell), n_steps):
        tgt = np.asarray(vp.controller.target, dtype=np.float64)
        vp.camera.local.position = tgt - dist * d
        vp.camera.look_at(tgt)
        vp.camera.depth_range = (dist * 1e-3, dist * 1e3)
        yield float(dist)
