"""Rectangular-duct H-grid primitive.

A single hex block filling a Cartesian box, with default BCs that
match the typical internal-flow / duct configuration:

- ``i_lo`` / ``i_hi`` (streamwise inflow / outflow) → ``FREE``.
- ``j_lo`` / ``j_hi`` / ``k_lo`` / ``k_hi`` (the four duct walls) → ``WALL``.

Geometrically the mesh is bit-identical to :class:`shore.primitive.BoxMesh`;
the only differences are the BC defaults and the user-facing
``streamwise_axis`` parameter.  Internally :func:`HGridRectMesh.build`
delegates to ``BoxMesh.build`` for the node array, then permutes the
axes so ``i`` becomes the streamwise direction and overrides the
default face BCs.

For more general duct geometries (variable cross-section, non-trivial
inner shapes), see the planned STL-driven H-grid topology
(``shore mesh --topology hgrid``).
"""

from __future__ import annotations

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC
from shore.mesh.spacing import Spacing1D
from shore.primitive.base import PrimitiveMesh
from shore.primitive.box import BoxMesh

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
# For each streamwise-axis label, the (cross_axis_a, cross_axis_b) labels
# in the order they will become the block's j and k axes.
_CROSS_AXES = {
    "x": ("y", "z"),
    "y": ("x", "z"),
    "z": ("x", "y"),
}


class HGridRectMesh(PrimitiveMesh):
    """Single-block rectangular H-grid (duct without an internal body).

    Node array shape: ``(ni, nj, nk, 3)`` where ``i`` is the streamwise
    direction.  All three axes can carry independent spacing laws
    (uniform / tanh / tanh2); tanh / tanh2 on the ``j`` and ``k`` axes
    typically clusters near the duct walls.

    Examples
    --------

    >>> # 1m-long square duct, 0.1m x 0.1m cross-section, streamwise = z.
    >>> hg = HGridRectMesh.build(
    ...     length=1.0, width=0.1, height=0.1,
    ...     ni=80, nj=20, nk=20,
    ...     streamwise_axis="z",
    ...     j_law="tanh2", j_beta=4.0,
    ...     k_law="tanh2", k_beta=4.0,
    ... )
    >>> hg.blocks[0].nodes.shape
    (80, 20, 20, 3)
    >>> hg.blocks[0].face("k_lo").bc.value
    'wall'
    """

    label_default = "h_grid_rect"

    def __init__(self, blocks: list[HexBlock]) -> None:
        self.blocks = blocks

    @classmethod
    def build(
        cls,
        length: float,
        width: float,
        height: float,
        ni: int,
        nj: int,
        nk: int,
        streamwise_axis: str = "x",
        i_law: str = "uniform",
        j_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        j_beta: float = 3.0,
        k_beta: float = 3.0,
        edge_pins: dict[str, Spacing1D] | None = None,
        label: str = "h_grid_rect",
    ) -> HGridRectMesh:
        """Build a rectangular-duct H-grid.

        Parameters
        ----------
        length :
            Duct length along the streamwise axis.  Must be positive.
        width, height :
            Cross-section dimensions, perpendicular to the streamwise
            axis.  Both must be positive.  ``width`` is along the
            block's ``j`` axis; ``height`` is along ``k``.
        ni, nj, nk :
            Node counts along the (streamwise, j, k) directions of
            the block.  All must be ``>= 2``.
        streamwise_axis :
            World axis the duct is aligned with: ``"x"``, ``"y"``, or
            ``"z"``.  The cross-section spans the other two axes.
            Default ``"x"``.
        i_law, j_law, k_law :
            Spacing laws for the three block axes
            (``"uniform"`` / ``"tanh"`` / ``"tanh2"``).  Defaults
            uniform.  Tanh / tanh2 on ``j`` and ``k`` cluster near the
            duct walls.
        i_beta, j_beta, k_beta :
            Clustering strengths for tanh / tanh2.
        edge_pins :
            Optional sparse per-edge spacing override; see
            :func:`shore.primitive.box.BoxMesh.build`.
        label :
            Block label written into the ``.geo`` filename when used
            in multi-block contexts.

        Returns
        -------
        HGridRectMesh
            Single-block primitive with default face BCs:
            ``i_lo`` / ``i_hi`` = ``FREE`` (inflow / outflow),
            ``j_lo`` / ``j_hi`` / ``k_lo`` / ``k_hi`` = ``WALL`` (duct walls).

        Raises
        ------
        ParameterError
            On invalid dimensions, node counts, or streamwise axis.
        """
        if streamwise_axis not in _AXIS_INDEX:
            raise ParameterError(
                f"HGridRectMesh: streamwise_axis must be 'x' / 'y' / 'z', got {streamwise_axis!r}"
            )
        for name, val in (("length", length), ("width", width), ("height", height)):
            if val <= 0:
                raise ParameterError(f"HGridRectMesh: {name} must be > 0, got {val}")

        # Build the block in *block-local* coordinates: i along streamwise,
        # j along width, k along height.  We use BoxMesh.build with vmin
        # at the origin and vmax at (length, width, height) so the node
        # array is in (ni, nj, nk) layout already.  Then we permute the
        # last (xyz) axis so the world-frame streamwise direction lands
        # on the user's chosen axis.
        box = BoxMesh.build(
            vmin=(0.0, 0.0, 0.0),
            vmax=(length, width, height),
            ni=ni,
            nj=nj,
            nk=nk,
            i_law=i_law,
            j_law=j_law,
            k_law=k_law,
            i_beta=i_beta,
            j_beta=j_beta,
            k_beta=k_beta,
            edge_pins=edge_pins,
            label=label,
        )
        block = box.blocks[0]

        # Permute the xyz components so that block-local i lands on the
        # user-chosen world streamwise axis.
        s_idx = _AXIS_INDEX[streamwise_axis]
        cross_a, cross_b = _CROSS_AXES[streamwise_axis]
        ca_idx = _AXIS_INDEX[cross_a]
        cb_idx = _AXIS_INDEX[cross_b]

        # Block-local layout: nodes[..., 0] = i-axis (= length),
        # nodes[..., 1] = j-axis (= width), nodes[..., 2] = k-axis (= height).
        # Map:  block local (length, width, height) → world (s_idx, ca_idx, cb_idx).
        permuted = np.empty_like(block.nodes)
        permuted[..., s_idx] = block.nodes[..., 0]
        permuted[..., ca_idx] = block.nodes[..., 1]
        permuted[..., cb_idx] = block.nodes[..., 2]
        block.nodes[...] = permuted

        # Override default face BCs: streamwise faces FREE, cross-section
        # walls WALL.  The default BoxMesh BC layout has k_lo = WALL and
        # everything else FREE; for H-grid we need FREE / FREE on the
        # streamwise (i) faces and WALL on the four cross-section faces.
        block._set_face_bc("i_lo", FaceBC.FREE)
        block._set_face_bc("i_hi", FaceBC.FREE)
        block._set_face_bc("j_lo", FaceBC.WALL)
        block._set_face_bc("j_hi", FaceBC.WALL)
        block._set_face_bc("k_lo", FaceBC.WALL)
        block._set_face_bc("k_hi", FaceBC.WALL)

        return cls(blocks=[block])
