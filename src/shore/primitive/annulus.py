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
from shore.mesh.spacing import Spacing1D
from shore.primitive.base import (
    PrimitiveMesh,
    axis_field_tfi,
    axis_frame,
    resolve_axis,
    resolve_edge_specs,
    stamp_edge_specs_onto_block,
)


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


_J_EDGE_NAMES = ("i_lo_k_lo", "i_hi_k_lo", "i_lo_k_hi", "i_hi_k_hi")


def _law_to_spacing1d(law: str, beta: float) -> Spacing1D:
    if law not in ("uniform", "tanh", "tanh2"):
        raise ParameterError(
            f"AnnularCylinderMesh: unknown spacing law {law!r} "
            "(expected 'uniform', 'tanh', 'tanh2')"
        )
    return Spacing1D(law=law, beta=beta)


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

    Per-edge spacing
    ----------------
    The ``edge_pins`` argument lets individual edges carry a different
    ``Spacing1D`` than the global ``i_law`` / ``j_law`` / ``k_law``.
    When any edge of an axis is pinned, that axis's parametric
    coordinate becomes a transfinite-interpolation blend of the 4
    parallel edges' distributions.

    Periodic-j constraint: for a full annulus, the j-direction is
    azimuthally periodic and pinning any j-edge (``i_lo_k_lo``,
    ``i_hi_k_lo``, ``i_lo_k_hi``, ``i_hi_k_hi``) is rejected — the
    "first" and "last" cells of a periodic edge are the same cell, so
    per-edge endpoint pinning is meaningless there.  Use a partial
    sector if you need azimuthal clustering control.

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
        edge_pins: dict[str, Spacing1D] | None = None,
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
            Node counts (radial, azimuthal, axial).  All >= 2.  These
            are the array shape of the produced ``(ni, nj, nk, 3)``
            block, not cell counts.  ``ni = 20`` means 20 radial nodes
            and 19 cells.
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
            Per-axis default spacing laws and clustering strengths.
        edge_pins :
            Sparse per-edge overrides; see class docstring.

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

        # Periodic-j cannot be per-edge-pinned: see class docstring.
        if periodic_j:
            if j_law != "uniform":
                raise ParameterError(
                    "AnnularCylinderMesh: j_law must be 'uniform' for a full "
                    "(periodic) annulus; partial sectors accept any law."
                )
            if edge_pins is not None:
                offending = [name for name in edge_pins if name in _J_EDGE_NAMES]
                if offending:
                    raise ParameterError(
                        "AnnularCylinderMesh: edge_pins on j-edges "
                        f"({offending}) is not allowed for a periodic full annulus. "
                        "Use a partial sector (set theta_end) for azimuthal "
                        "per-edge spacing control."
                    )

        # Per-axis defaults.
        i_default = _law_to_spacing1d(i_law, i_beta)
        k_default = _law_to_spacing1d(k_law, k_beta)

        # ---------------------------------------------------------------
        # Radial parametric field (i-axis)
        # ---------------------------------------------------------------
        i_specs = resolve_edge_specs("i", edge_pins, i_default)
        r_field = axis_field_tfi(
            ni,
            nj,
            nk,
            extent.r_in,
            extent.r_out,
            edge_specs=i_specs,
            axis_label="i",
        )  # (ni, nj, nk)

        # ---------------------------------------------------------------
        # Azimuthal parametric field (j-axis)
        # ---------------------------------------------------------------
        if periodic_j:
            # Full circle: no pinning, no TFI.  Drop the duplicate endpoint.
            thetas_1d = np.linspace(theta_start, theta_end_eff, nj, endpoint=False)
            theta_field = np.broadcast_to(thetas_1d[None, :, None], (ni, nj, nk)).copy()
        else:
            j_default = _law_to_spacing1d(j_law, j_beta)
            j_specs = resolve_edge_specs("j", edge_pins, j_default)
            theta_field_local = axis_field_tfi(
                nj,
                ni,
                nk,
                theta_start,
                theta_end_eff,
                edge_specs=j_specs,
                axis_label="j",
            )  # (nj, ni, nk)
            theta_field = np.transpose(theta_field_local, (1, 0, 2))  # (ni, nj, nk)

        # ---------------------------------------------------------------
        # Axial parametric field (k-axis)
        # ---------------------------------------------------------------
        k_specs = resolve_edge_specs("k", edge_pins, k_default)
        z_field_local = axis_field_tfi(
            nk,
            ni,
            nj,
            extent.z0,
            extent.z1,
            edge_specs=k_specs,
            axis_label="k",
        )  # (nk, ni, nj)
        z_field = np.transpose(z_field_local, (1, 2, 0))  # (ni, nj, nk)

        # ---------------------------------------------------------------
        # Local frame and Cartesian assembly.
        # ---------------------------------------------------------------
        axis_vec = resolve_axis(axis)
        e1, e2, e3 = axis_frame(axis_vec)
        c = np.asarray(center, dtype=float).reshape(1, 1, 1, 3)

        cosT = np.cos(theta_field)
        sinT = np.sin(theta_field)
        # nodes = center + r*cosT*e1 + r*sinT*e2 + z*e3
        nodes = (
            c
            + (r_field * cosT)[..., None] * e1
            + (r_field * sinT)[..., None] * e2
            + z_field[..., None] * e3
        )

        block = HexBlock(
            label=label,
            nodes=nodes,
            periodic_i=False,
            periodic_j=periodic_j,
        )

        stamp_edge_specs_onto_block(block, edge_pins=edge_pins)
        return cls(blocks=[block])
