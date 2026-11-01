"""Circular-bore H-grid primitive (ADR 0060).

A straight circular tube filled by five hexahedral blocks in an O-H
("butterfly") cross-section:

- ``core``: a square in the middle, its corners at ``core_frac`` of the
  radius, a Cartesian grid;
- ``ring_0`` .. ``ring_3``: four blocks joining the core's sides to the
  outer wall, each a quarter of the circle. Their lines run straight from
  the core's side to the wall, along the 45 deg diagonals at the corners.

No singular axis: every cell is a proper hexahedron, and every block is
right-handed whatever the streamwise axis. Faces:

- ``i_lo`` / ``i_hi`` (streamwise inflow / outflow) of every block: ``FREE``;
- the rings' ``k_hi`` (the outer wall): ``WALL``;
- every face between two blocks: ``SHARED``, with its axis map.

Before ADR 0060 the tube was one block on a polar grid: its ``k = 0`` row
collapsed onto the centreline, which was marked ``WALL``; its first ring
of cells was degenerate and the block was left-handed for the x and z
axes.
"""

from __future__ import annotations

import warnings

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC
from shore.primitive.base import PrimitiveMesh
from shore.volume.spacing import build_i_coords

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
# (cross_a, cross_b) for each streamwise axis: cyclic, so (stream, a, b)
# is a right-handed frame.
_CROSS_AXES = {
    "x": ("y", "z"),
    "y": ("z", "x"),
    "z": ("x", "y"),
}

#: The tube's blocks, in order: the four rings, then the core.
TUBE_BLOCKS: tuple[str, ...] = ("ring_0", "ring_1", "ring_2", "ring_3", "core")


class HGridTubeMesh(PrimitiveMesh):
    """Circular-bore H-grid: a straight tube, O-H cross-section, five blocks.

    Node arrays: ``ring_q`` is ``(n_axial, n_circ // 4 + 1, n_radial, 3)``,
    i streamwise, j along its quarter of the wall, k from the core's side
    (k = 0) to the wall (k = n_radial - 1); ``core`` is ``(n_axial,
    n_circ // 4 + 1, n_circ // 4 + 1, 3)``. ``ring_q`` covers the angles
    -45 + 90 q to 45 + 90 q degrees in the cross-section's plane.

    Examples
    --------

    >>> # 1m-long tube of radius 0.05m with 60 cells round the wall.
    >>> ht = HGridTubeMesh.build(
    ...     length=1.0, radius=0.05,
    ...     n_axial=80, n_circ=60, n_radial=12,
    ...     streamwise_axis="z",
    ...     k_law="tanh", k_beta=3.0,  # cluster the ring's nodes at the wall
    ... )
    >>> [b.label for b in ht.blocks]
    ['ring_0', 'ring_1', 'ring_2', 'ring_3', 'core']
    >>> ht.blocks[0].nodes.shape
    (80, 16, 12, 3)
    """

    label_default = "h_grid_tube"

    def __init__(self, blocks: list[HexBlock]) -> None:
        self.blocks = blocks

    @classmethod
    def build(
        cls,
        length: float,
        radius: float,
        n_axial: int,
        n_circ: int,
        n_radial: int,
        streamwise_axis: str = "x",
        i_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        k_beta: float = 3.0,
        label: str | None = None,
        core_frac: float = 0.5,
    ) -> HGridTubeMesh:
        """Build a circular-bore H-grid.

        Parameters
        ----------
        length :
            Tube length along the streamwise axis.  Must be positive.
        radius :
            Tube outer radius.  Must be positive.
        n_axial :
            Streamwise node count, ``>= 2``.
        n_circ :
            Cells round the wall, a multiple of 4 and ``>= 4``: each ring
            block takes a quarter, and the core has ``n_circ // 4`` cells a
            side.
        n_radial :
            The ring blocks' node count from the core to the wall, ``>= 2``.
        streamwise_axis :
            World axis the tube is aligned with: ``"x"`` (default),
            ``"y"``, or ``"z"``.  The cross-section spans the other two.
        i_law, k_law :
            Spacing laws (``"uniform"`` / ``"tanh"`` / ``"tanh2"``) along
            the stream (``i``) and across the rings, from the core to the
            wall (``k``). ``"tanh"`` clusters at the start of the stream,
            and at the wall across the rings.
        i_beta, k_beta :
            Clustering strengths for tanh / tanh2.
        label :
            Ignored since ADR 0060 (a ``DeprecationWarning`` when given):
            the blocks are ``TUBE_BLOCKS``.
        core_frac :
            Radius of the core's corners as a fraction of the tube's, in
            (0, 1). Default ``0.5``.

        Returns
        -------
        HGridTubeMesh
            Five right-handed blocks (``TUBE_BLOCKS``).

        Raises
        ------
        ParameterError
            On invalid dimensions or node counts.
        """
        if label is not None:
            warnings.warn(
                "HGridTubeMesh.build: label is ignored since ADR 0060; the blocks "
                "are ring_0 .. ring_3 and core.",
                DeprecationWarning,
                stacklevel=2,
            )
        if streamwise_axis not in _AXIS_INDEX:
            raise ParameterError(
                f"HGridTubeMesh: streamwise_axis must be 'x' / 'y' / 'z', got {streamwise_axis!r}"
            )
        if length <= 0:
            raise ParameterError(f"HGridTubeMesh: length must be > 0, got {length}")
        if radius <= 0:
            raise ParameterError(f"HGridTubeMesh: radius must be > 0, got {radius}")
        if n_axial < 2:
            raise ParameterError(f"HGridTubeMesh: n_axial must be >= 2, got {n_axial}")
        if n_circ < 4 or n_circ % 4:
            raise ParameterError(
                f"HGridTubeMesh: n_circ must be a multiple of 4 and >= 4 (one "
                f"quarter per ring block), got {n_circ}"
            )
        if n_radial < 2:
            raise ParameterError(f"HGridTubeMesh: n_radial must be >= 2, got {n_radial}")
        if not 0.0 < core_frac < 1.0:
            raise ParameterError(f"HGridTubeMesh: core_frac must be in (0, 1), got {core_frac}")
        for law in (i_law, k_law):
            if law not in ("uniform", "tanh", "tanh2"):
                raise ParameterError(f"HGridTubeMesh: unknown spacing law {law!r}")

        m = n_circ // 4
        s = build_i_coords(n_axial, 0.0, length, spacing=i_law, beta=i_beta)
        raw = build_i_coords(n_radial, 0.0, 1.0, spacing=k_law, beta=k_beta)
        t_k = 1.0 - raw[::-1] if k_law == "tanh" else raw  # tanh: cluster at the wall

        # The cross-section, in the (a, b) plane of the frame (stream, a, b),
        # from one wall polygon and one core perimeter that every block
        # slices, so shared nodes are the same floats.
        u = np.linspace(
            -core_frac * radius / np.sqrt(2.0), core_frac * radius / np.sqrt(2.0), m + 1
        )
        lo, hi, rev = u[0], u[-1], u[::-1]
        square = np.concatenate(
            [
                np.stack([np.full(m, hi), u[:-1]], axis=-1),  # x = +half, y up
                np.stack([rev[:-1], np.full(m, hi)], axis=-1),  # y = +half, x down
                np.stack([np.full(m, lo), rev[:-1]], axis=-1),  # x = -half, y down
                np.stack([u[:-1], np.full(m, lo)], axis=-1),  # y = -half, x up
            ]
        )
        theta = np.radians(-45.0) + np.arange(4 * m) * (2.0 * np.pi / (4 * m))
        wall = radius * np.stack([np.cos(theta), np.sin(theta)], axis=-1)
        square, wall = np.vstack([square, square[:1]]), np.vstack([wall, wall[:1]])
        rings2d = []
        for q in range(4):
            side, arc = square[q * m : (q + 1) * m + 1], wall[q * m : (q + 1) * m + 1]
            ring = (1.0 - t_k)[None, :, None] * side[:, None] + t_k[None, :, None] * arc[:, None]
            ring[:, 0], ring[:, -1] = side, arc
            rings2d.append(ring)  # (m + 1, n_radial, 2): j along the wall, k outwards
        ua, ub = np.meshgrid(u, u, indexing="ij")
        core2d = np.stack([ua, ub], axis=-1)  # (m + 1, m + 1, 2), axes x and y
        core2d[-1, :] = square[0 : m + 1]
        core2d[:, -1] = square[m : 2 * m + 1][::-1]
        core2d[0, :] = square[2 * m : 3 * m + 1][::-1]
        core2d[:, 0] = square[3 * m : 4 * m + 1]

        a_lbl, b_lbl = _CROSS_AXES[streamwise_axis]
        idx = (_AXIS_INDEX[streamwise_axis], _AXIS_INDEX[a_lbl], _AXIS_INDEX[b_lbl])

        def lift(plane: np.ndarray) -> np.ndarray:
            nodes = np.empty((n_axial, *plane.shape[:2], 3), dtype=np.float64)
            nodes[..., idx[0]] = s[:, None, None]
            nodes[..., idx[1]] = plane[None, ..., 0]
            nodes[..., idx[2]] = plane[None, ..., 1]
            e_i = nodes[1, 0, 0] - nodes[0, 0, 0]
            e_j = nodes[0, 1, 0] - nodes[0, 0, 0]
            e_k = nodes[0, 0, 1] - nodes[0, 0, 0]
            right = float(np.cross(e_i, e_j) @ e_k) > 0.0
            return nodes if right else np.ascontiguousarray(nodes[:, ::-1])

        from shore.mesh.seams import wire_by_corners

        blocks = [
            HexBlock(label=name, nodes=lift(plane), periodic_i=False, periodic_j=False)
            for name, plane in zip(TUBE_BLOCKS, [*rings2d, core2d], strict=True)
        ]
        wire_by_corners(blocks, ("j_lo", "j_hi", "k_lo", "k_hi"))
        for b in blocks:
            b._set_face_bc("i_lo", FaceBC.FREE)  # inflow
            b._set_face_bc("i_hi", FaceBC.FREE)  # outflow
        for b in blocks[:4]:
            b._set_face_bc("k_hi", FaceBC.WALL)  # the tube's wall
        return cls(blocks=blocks)
