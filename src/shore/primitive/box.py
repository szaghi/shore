"""Cartesian box primitive — a single axis-aligned hexahedral block."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.spacing import Spacing1D
from shore.primitive.base import (
    PrimitiveMesh,
    axis_field_tfi,
    resolve_edge_specs,
    stamp_edge_specs_onto_block,
)


@dataclass(frozen=True)
class BoxExtent:
    """Min/max corners of an axis-aligned box.

    The min corner is the (i=0, j=0, k=0) node; the max corner is the
    (ni-1, nj-1, nk-1) node.  All three components of *vmax* must
    strictly exceed *vmin*.
    """

    vmin: tuple[float, float, float]
    vmax: tuple[float, float, float]

    def __post_init__(self) -> None:
        vmin = np.asarray(self.vmin, dtype=float)
        vmax = np.asarray(self.vmax, dtype=float)
        if vmin.shape != (3,) or vmax.shape != (3,):
            raise ParameterError(
                f"BoxExtent: vmin/vmax must be 3-vectors, got {vmin.shape}, {vmax.shape}"
            )
        if not np.all(vmax > vmin):
            raise ParameterError(
                f"BoxExtent: vmax must exceed vmin in every component "
                f"(got vmin={tuple(vmin)}, vmax={tuple(vmax)})"
            )


def _law_to_spacing1d(law: str, beta: float) -> Spacing1D:
    """Map a string law name to a Spacing1D instance for the axis default."""
    if law not in ("uniform", "tanh", "tanh2"):
        raise ParameterError(
            f"BoxMesh: unknown spacing law {law!r} (expected 'uniform', 'tanh', 'tanh2')"
        )
    return Spacing1D(law=law, beta=beta)


class BoxMesh(PrimitiveMesh):
    """Single-block axis-aligned Cartesian box.

    Node array shape: ``(ni, nj, nk, 3)``.  Each axis has its own
    spacing law (uniform / tanh / tanh2) by default.  Per-edge
    overrides via ``edge_pins`` switch the construction to transfinite
    interpolation between the 4 parallel edges of the affected axis.
    No projection, no marching.

    Examples
    --------
    Tensor-product (legacy) construction with global per-axis laws:

    >>> box = BoxMesh.build(
    ...     vmin=(0.0, 0.0, 0.0),
    ...     vmax=(10.0, 5.0, 3.0),
    ...     ni=40, nj=30, nk=20,
    ... )
    >>> box.blocks[0].nodes.shape
    (40, 30, 20, 3)

    Per-edge override (one i-edge gets tanh clustering, others stay
    uniform — interior nodes blend the four edge distributions):

    >>> from shore.mesh.spacing import Spacing1D
    >>> box = BoxMesh.build(
    ...     vmin=(0.0, 0.0, 0.0), vmax=(10.0, 5.0, 3.0),
    ...     ni=40, nj=30, nk=20,
    ...     edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
    ... )
    """

    label_default = "box"

    def __init__(self, blocks: list[HexBlock]) -> None:
        self.blocks = blocks

    @classmethod
    def build(
        cls,
        vmin: tuple[float, float, float],
        vmax: tuple[float, float, float],
        ni: int,
        nj: int,
        nk: int,
        i_law: str = "uniform",
        j_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        j_beta: float = 3.0,
        k_beta: float = 3.0,
        edge_pins: dict[str, Spacing1D] | None = None,
        label: str = "box",
    ) -> BoxMesh:
        """Build a box primitive.

        Parameters
        ----------
        vmin, vmax :
            3-tuples giving the min and max corners (world coordinates).
            Each component of *vmax* must strictly exceed the
            corresponding component of *vmin*.
        ni, nj, nk :
            Node counts along the three axes.  All must be >= 2.
        i_law, j_law, k_law :
            Default spacing laws for each axis (``"uniform"``,
            ``"tanh"``, ``"tanh2"``).  Used for any of the 4 parallel
            edges of an axis that aren't explicitly overridden in
            *edge_pins*.
        i_beta, j_beta, k_beta :
            Clustering strengths for tanh / tanh2.
        edge_pins :
            Optional sparse override mapping edge name → ``Spacing1D``.
            Edge names are the canonical 12 in :mod:`shore.mesh.edge`
            (e.g. ``"j_lo_k_lo"`` for the i-axis edge at j=0, k=0).
            When any edge is pinned, construction switches to
            transfinite interpolation between the 4 parallel edges of
            that edge's axis.  Edges from other axes that aren't
            mentioned use the global ``*_law`` default.
        label :
            Block label written into the ``.geo`` filename when used in
            multi-block contexts.

        Returns
        -------
        BoxMesh
        """
        extent = BoxExtent(vmin=vmin, vmax=vmax)
        for n, name in ((ni, "ni"), (nj, "nj"), (nk, "nk")):
            if n < 2:
                raise ParameterError(f"BoxMesh: {name} must be >= 2, got {n}")

        i_default = _law_to_spacing1d(i_law, i_beta)
        j_default = _law_to_spacing1d(j_law, j_beta)
        k_default = _law_to_spacing1d(k_law, k_beta)

        # Per-axis edge specs.  When edge_pins is None or empty, every
        # spec falls back to the axis default and TFI reduces to the old
        # tensor-product construction.
        i_specs = resolve_edge_specs("i", edge_pins, i_default)
        j_specs = resolve_edge_specs("j", edge_pins, j_default)
        k_specs = resolve_edge_specs("k", edge_pins, k_default)

        # Build each Cartesian component as a 3D field via TFI.  When all
        # 4 parallel edges of an axis agree, the field collapses to the
        # legacy 1D distribution broadcast across the other two axes.
        x_field = axis_field_tfi(
            ni,
            nj,
            nk,
            extent.vmin[0],
            extent.vmax[0],
            edge_specs=i_specs,
            axis_label="i",
        )  # (ni, nj, nk)
        y_field = axis_field_tfi(
            nj,
            ni,
            nk,
            extent.vmin[1],
            extent.vmax[1],
            edge_specs=j_specs,
            axis_label="j",
        )  # (nj, ni, nk)  — note the axis ordering inside the helper
        z_field = axis_field_tfi(
            nk,
            ni,
            nj,
            extent.vmin[2],
            extent.vmax[2],
            edge_specs=k_specs,
            axis_label="k",
        )  # (nk, ni, nj)

        # Reorder y_field and z_field to (ni, nj, nk) and stack.
        nodes = np.empty((ni, nj, nk, 3), dtype=np.float64)
        nodes[..., 0] = x_field  # x: (ni, nj, nk)
        nodes[..., 1] = np.transpose(y_field, (1, 0, 2))  # y: (nj, ni, nk) → (ni, nj, nk)
        nodes[..., 2] = np.transpose(z_field, (1, 2, 0))  # z: (nk, ni, nj) → (ni, nj, nk)

        # Box has no periodic direction.
        block = HexBlock(label=label, nodes=nodes, periodic_i=False, periodic_j=False)

        # Stamp the user's per-edge specs onto the block so downstream
        # readers can recover the schedule that was actually used.
        stamp_edge_specs_onto_block(block, edge_pins=edge_pins)

        return cls(blocks=[block])
