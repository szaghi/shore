"""Circular-bore H-grid primitive.

A single hex block filling a straight circular tube, with default BCs
that match the typical internal-flow / duct configuration:

- ``i_lo`` / ``i_hi`` (streamwise inflow / outflow) → ``FREE``.
- ``j_lo`` / ``j_hi`` (circumferential) → SHARED via ``periodic_j``.
- ``k_lo`` (centerline singular axis) → ``WALL`` (degenerate).
- ``k_hi`` (outer wall) → ``WALL``.

Cross-section is a 2D polar grid (lat-lon style): ``j`` is the
circumferential angle and ``k`` is the radial coordinate.  The
``k = 0`` row is the **polar singular axis** at the duct centerline —
all ``n_circ`` angular nodes collapse to a single point.  This is the
analogue of the polar axis row on
:class:`shore.mesh.OGridTopology` and is acceptable because Xall's
sphere-grid code handles polar singularities cleanly.

For non-trivial duct geometries (variable cross-section, non-circular
bores), see the planned STL-driven H-grid topology
(``shore mesh --topology hgrid``).
"""

from __future__ import annotations

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC
from shore.primitive.base import PrimitiveMesh
from shore.volume.spacing import build_i_coords

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
# (cross_a, cross_b) labels for each streamwise-axis choice.
_CROSS_AXES = {
    "x": ("y", "z"),
    "y": ("x", "z"),
    "z": ("x", "y"),
}


class HGridTubeMesh(PrimitiveMesh):
    """Single-block circular-bore H-grid (straight tube, no internal body).

    Node array shape: ``(n_axial, n_circ, n_radial, 3)``.  The block's
    ``i`` axis is streamwise, ``j`` is the circumferential angle
    (periodic by construction), and ``k`` is the radial coordinate
    from the centerline (``k = 0``) to the outer wall (``k = n_radial - 1``).

    Examples
    --------

    >>> # 1m-long tube of radius 0.05m with 60 circumferential cells.
    >>> ht = HGridTubeMesh.build(
    ...     length=1.0, radius=0.05,
    ...     n_axial=80, n_circ=60, n_radial=12,
    ...     streamwise_axis="z",
    ...     k_law="tanh", k_beta=3.0,  # cluster radial nodes near the wall
    ... )
    >>> ht.blocks[0].nodes.shape
    (80, 60, 12, 3)
    >>> ht.blocks[0].periodic_j
    True
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
        label: str = "h_grid_tube",
    ) -> HGridTubeMesh:
        """Build a circular-bore H-grid.

        Parameters
        ----------
        length :
            Tube length along the streamwise axis.  Must be positive.
        radius :
            Tube outer radius.  Must be positive.
        n_axial, n_circ, n_radial :
            Node counts along the streamwise, circumferential, and
            radial directions.  ``n_axial >= 2``, ``n_circ >= 4``,
            ``n_radial >= 2``.
        streamwise_axis :
            World axis the tube is aligned with: ``"x"`` (default),
            ``"y"``, or ``"z"``.  The cross-section spans the other
            two axes.
        i_law, k_law :
            Spacing laws for the streamwise (``i``) and radial (``k``)
            axes (``"uniform"`` / ``"tanh"`` / ``"tanh2"``).
            Circumferential (``j``) is always uniform — clustering on
            a periodic axis is unsupported.
        i_beta, k_beta :
            Clustering strengths for tanh / tanh2.
        label :
            Block label.

        Returns
        -------
        HGridTubeMesh
            Single-block primitive with default face BCs:
            ``i_lo`` / ``i_hi`` = ``FREE`` (inflow / outflow),
            ``j`` direction = SHARED (periodic), ``k_lo`` = ``WALL``
            (centerline singular axis), ``k_hi`` = ``WALL`` (outer wall).

        Raises
        ------
        ParameterError
            On invalid dimensions or node counts.
        """
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
        if n_circ < 4:
            raise ParameterError(
                f"HGridTubeMesh: n_circ must be >= 4 (a quad's worth of "
                f"circumferential cells), got {n_circ}"
            )
        if n_radial < 2:
            raise ParameterError(f"HGridTubeMesh: n_radial must be >= 2, got {n_radial}")

        # ------------------------------------------------------------
        # 1D coordinates per axis (block-local).
        # ------------------------------------------------------------
        # Streamwise: 0..length with the user's law.
        s_coords = build_i_coords(n_axial, 0.0, length, spacing=i_law, beta=i_beta)
        # Radial: 0..radius with the user's law.  k=0 is the singular
        # centerline; k=n_radial-1 is the outer wall.
        r_coords = build_i_coords(n_radial, 0.0, radius, spacing=k_law, beta=k_beta)
        # Circumferential: uniform 0..2*pi (periodic, drop the last node).
        theta = np.linspace(0.0, 2.0 * np.pi, n_circ, endpoint=False)

        # ------------------------------------------------------------
        # Build the (n_axial, n_circ, n_radial, 3) node array in
        # block-local coordinates: i along streamwise, the polar grid
        # (cos*r, sin*r) lives in the (j=theta, k=radius) plane.
        # ------------------------------------------------------------
        s_grid = s_coords[:, None, None]  # (n_axial, 1, 1)
        cos_t = np.cos(theta)[None, :, None]  # (1, n_circ, 1)
        sin_t = np.sin(theta)[None, :, None]  # (1, n_circ, 1)
        r_grid = r_coords[None, None, :]  # (1, 1, n_radial)

        cross_a = (cos_t * r_grid).repeat(n_axial, axis=0)
        cross_b = (sin_t * r_grid).repeat(n_axial, axis=0)
        streamwise = s_grid.repeat(n_circ, axis=1).repeat(n_radial, axis=2)

        # ------------------------------------------------------------
        # Permute axes so block-local (streamwise, cross_a, cross_b)
        # land on the user's chosen world axes.
        # ------------------------------------------------------------
        s_idx = _AXIS_INDEX[streamwise_axis]
        cross_a_label, cross_b_label = _CROSS_AXES[streamwise_axis]
        ca_idx = _AXIS_INDEX[cross_a_label]
        cb_idx = _AXIS_INDEX[cross_b_label]

        nodes = np.empty((n_axial, n_circ, n_radial, 3), dtype=np.float64)
        nodes[..., s_idx] = streamwise
        nodes[..., ca_idx] = cross_a
        nodes[..., cb_idx] = cross_b

        # ------------------------------------------------------------
        # Build the HexBlock with j-periodic and override default BCs.
        # ------------------------------------------------------------
        block = HexBlock(
            label=label,
            nodes=nodes,
            periodic_i=False,
            periodic_j=True,
        )
        block._set_face_bc("i_lo", FaceBC.FREE)  # inflow
        block._set_face_bc("i_hi", FaceBC.FREE)  # outflow
        # j_lo / j_hi: defaults already handle periodicity (no BC needed).
        block._set_face_bc("k_lo", FaceBC.WALL)  # centerline singular axis
        block._set_face_bc("k_hi", FaceBC.WALL)  # outer wall

        return cls(blocks=[block])
