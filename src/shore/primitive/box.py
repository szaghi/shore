"""Cartesian box primitive — a single axis-aligned hexahedral block."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.primitive.base import PrimitiveMesh, coords_1d


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


class BoxMesh(PrimitiveMesh):
    """Single-block axis-aligned Cartesian box.

    Node array shape: ``(ni, nj, nk, 3)``.  Each axis has its own
    spacing law (uniform / tanh / tanh2).  No projection, no marching;
    construction is a tensor product of three 1D distributions.

    Examples
    --------
    >>> box = BoxMesh.build(
    ...     vmin=(0.0, 0.0, 0.0),
    ...     vmax=(10.0, 5.0, 3.0),
    ...     ni=40, nj=30, nk=20,
    ... )
    >>> box.blocks[0].nodes.shape
    (40, 30, 20, 3)
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
            Spacing laws for each axis (``"uniform"``, ``"tanh"``,
            ``"tanh2"``).
        i_beta, j_beta, k_beta :
            Clustering strengths for tanh / tanh2.
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

        xs = coords_1d(ni, extent.vmin[0], extent.vmax[0], law=i_law, beta=i_beta)
        ys = coords_1d(nj, extent.vmin[1], extent.vmax[1], law=j_law, beta=j_beta)
        zs = coords_1d(nk, extent.vmin[2], extent.vmax[2], law=k_law, beta=k_beta)

        # Assemble (ni, nj, nk, 3).  np.meshgrid with indexing='ij'
        # gives the right axis ordering on first try.
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
        nodes = np.stack([X, Y, Z], axis=-1)

        # Box has no periodic direction.
        block = HexBlock(label=label, nodes=nodes, periodic_i=False, periodic_j=False)
        return cls(blocks=[block])
