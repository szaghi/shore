"""Annular cylinder primitive — a single block with optional theta-periodicity.

Geometry: an annular cylinder bounded by an inner radius, an outer
radius, and two axial planes.  The structured grid is parametrised by
``(r, theta, z)``; the *theta* direction can be either fully periodic
(``2*pi``) or a partial sector ``[theta_start, theta_end]``.

Index convention (i,j,k) → (radial, azimuthal, axial).  This matches
the legacy single-block O-grid: the j direction is the periodic /
sector direction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.primitive.base import PrimitiveMesh, axis_frame, coords_1d, resolve_axis


@dataclass(frozen=True)
class AnnulusExtent:
    """Inner radius, outer radius, axial extent of an annular cylinder."""

    r_in: float
    r_out: float
    z0: float
    z1: float

    def __post_init__(self) -> None:
        if self.r_in < 0.0:
            raise ParameterError(f"AnnulusExtent: r_in must be >= 0, got {self.r_in}")
        if self.r_out <= self.r_in:
            raise ParameterError(
                f"AnnulusExtent: r_out ({self.r_out}) must exceed r_in ({self.r_in})"
            )
        if self.z1 <= self.z0:
            raise ParameterError(f"AnnulusExtent: z1 ({self.z1}) must exceed z0 ({self.z0})")


class AnnularCylinderMesh(PrimitiveMesh):
    """Single-block annular cylinder.

    Index meaning
    -------------
    - ``i`` (axis 0): radial direction, ``r_in`` → ``r_out``.
    - ``j`` (axis 1): azimuthal direction, ``theta_start`` → ``theta_end``.
    - ``k`` (axis 2): axial direction along *axis* unit vector.

    For a full annulus (``theta_end - theta_start == 2*pi``) the j
    direction is periodic and the last azimuthal node coincides with
    the first.  For a partial sector j is non-periodic.

    The cylinder axis defaults to z.  Pass ``axis="x"``, ``axis="y"``,
    or a 3-vector to rotate the geometry into world frame.

    Examples
    --------
    >>> ann = AnnularCylinderMesh.build(
    ...     r_in=0.5, r_out=2.0, z0=0.0, z1=4.0,
    ...     ni=20, nj=64, nk=30,
    ... )
    >>> ann.blocks[0].periodic_j
    True
    """

    label_default = "annulus"

    def __init__(self, blocks: list[HexBlock]) -> None:
        self.blocks = blocks

    @classmethod
    def build(
        cls,
        r_in: float,
        r_out: float,
        z0: float,
        z1: float,
        ni: int,
        nj: int,
        nk: int,
        theta_start: float = 0.0,
        theta_end: float | None = None,
        axis: str | tuple[float, float, float] = "z",
        center: tuple[float, float, float] = (0.0, 0.0, 0.0),
        i_law: str = "uniform",
        j_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        j_beta: float = 3.0,
        k_beta: float = 3.0,
        label: str = "annulus",
    ) -> AnnularCylinderMesh:
        """Build an annular cylinder primitive.

        Parameters
        ----------
        r_in, r_out :
            Inner and outer radii (``r_out > r_in >= 0``).
        z0, z1 :
            Axial extent along *axis* (``z1 > z0``).
        ni, nj, nk :
            Node counts (radial, azimuthal, axial).  All >= 2.
        theta_start, theta_end :
            Azimuthal sector in radians.  When *theta_end* is None,
            the cylinder is a full annulus (``theta_start`` to
            ``theta_start + 2*pi``) and the j direction is periodic.
            When given explicitly, the sector is partial and j is
            non-periodic; ``theta_end > theta_start`` is required.
        axis :
            ``"x"``, ``"y"``, ``"z"`` (case-insensitive), or a
            3-component vector giving an arbitrary axis direction.
        center :
            World-space anchor for the cylinder axis.  The axis line
            passes through *center* with direction *axis*; ``z=0`` of
            the local frame corresponds to *center*.
        i_law, j_law, k_law, *beta :
            Per-axis spacing laws and clustering strengths.

        Returns
        -------
        AnnularCylinderMesh
        """
        extent = AnnulusExtent(r_in=r_in, r_out=r_out, z0=z0, z1=z1)
        for n, name in ((ni, "ni"), (nj, "nj"), (nk, "nk")):
            if n < 2:
                raise ParameterError(f"AnnularCylinderMesh: {name} must be >= 2, got {n}")

        # Resolve theta range and periodicity.
        if theta_end is None:
            theta_end_eff = theta_start + 2.0 * np.pi
            periodic_j = True
        else:
            if theta_end <= theta_start:
                raise ParameterError(
                    f"AnnularCylinderMesh: theta_end ({theta_end}) must exceed "
                    f"theta_start ({theta_start})"
                )
            theta_end_eff = theta_end
            full = abs((theta_end - theta_start) - 2.0 * np.pi) < 1e-12
            periodic_j = full

        # Build the three 1D distributions.
        rs = coords_1d(ni, extent.r_in, extent.r_out, law=i_law, beta=i_beta)
        if periodic_j:
            # For a full annulus the last azimuthal node is the wrap of
            # the first — drop the duplicate by sampling [start, end)
            # via uniform spacing (no nonlinear law at the wrap point).
            if j_law != "uniform":
                raise ParameterError(
                    "AnnularCylinderMesh: j_law must be 'uniform' for a full "
                    "(periodic) annulus; partial sectors accept any law."
                )
            thetas = np.linspace(theta_start, theta_end_eff, nj, endpoint=False)
        else:
            thetas = coords_1d(nj, theta_start, theta_end_eff, law=j_law, beta=j_beta)

        zs = coords_1d(nk, extent.z0, extent.z1, law=k_law, beta=k_beta)

        # Local frame: e1, e2 span the radial plane; e3 is the axis.
        axis_vec = resolve_axis(axis)
        e1, e2, e3 = axis_frame(axis_vec)
        c = np.asarray(center, dtype=float).reshape(1, 1, 1, 3)

        # Assemble the (ni, nj, nk, 3) node array.
        R = rs[:, None, None]  # (ni, 1, 1)
        T = thetas[None, :, None]  # (1, nj, 1)
        Z = zs[None, None, :]  # (1, 1, nk)

        cosT = np.cos(T)
        sinT = np.sin(T)
        # nodes = center + R*cosT*e1 + R*sinT*e2 + Z*e3
        nodes = c + (R * cosT)[..., None] * e1 + (R * sinT)[..., None] * e2 + Z[..., None] * e3

        block = HexBlock(
            label=label,
            nodes=nodes,
            periodic_i=False,
            periodic_j=periodic_j,
        )
        return cls(blocks=[block])
