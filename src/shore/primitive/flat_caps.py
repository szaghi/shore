"""Flat-caps cylinder primitive — 5-block conforming alternative to a
hollow Chimera background.

Geometry
--------

At every axial slice the disk of radius ``r_out`` is decomposed into:

- 1 *central square* block, axis-aligned in the (e1, e2) plane and
  centred on the axis.  Its 4 corners sit at radius ``r_inner = inner_frac
  * r_out`` (default ``inner_frac = 0.5``) at azimuths
  ``pi/4 + k*pi/2`` (k = 0..3) — i.e. the "diamond" orientation.
  Side length ``s = sqrt(2) * r_inner``.
- 4 *trapezoidal sub-blocks* around the central square.  Each sub-block
  has:
    - inner edge: one side of the central square (a straight line of
      length ``s``, both endpoints at radius ``r_inner``);
    - outer edge: the matching 90° arc on the outer circle at
      radius ``r_out``;
    - two side edges: straight radial segments connecting a corner of
      the central square at radius ``r_inner`` to the same azimuth on
      the outer circle at radius ``r_out``.

The decomposition is then extruded along *axis* with ``nk`` axial
layers.  This produces 5 blocks total — no hole at the centre, no
overset fringes between blocks (all seams are SHARED structured
boundaries with bit-exact node coincidence).

C1 first-step continuity
------------------------

The radial direction inside each sub-block has its first cell
**pinned** to the central square's first cell (in the corresponding
direction), so that the radial cell-size is C1-continuous across the
central / sub-block seam:

    central_first_step = sqrt(2) * r_inner / (nc - 1)   (for c_law='uniform')
    sub_first_radial   = central_first_step             (always pinned)

The remaining ``ni - 2`` cells of the sub-block radial distribution
absorb the difference, distributed by the chosen ``i_law``.  When the
pin is incompatible with ``i_law``/``ni`` (would require a non-monotone
or trivially-empty interior), SHORE falls back to the unpinned
distribution and emits a ``UserWarning``.

Use case
--------

When you need a cylindrical background mesh but want to avoid Chimera
fringes (high-order accuracy, conservation, GPU memory access), this
primitive replaces the
``annular cylinder + hole`` overset combination with a single
conforming structured assembly.

Block labels
------------

- ``center``           — central square block.
- ``sub_e``, ``sub_n``, ``sub_w``, ``sub_s`` — surrounding
  sub-blocks (east / north / west / south in the (e1, e2) plane).

Each sub-block uses index ``i`` for the radial direction (square edge
at ``i = 0``, circle arc at ``i = ni - 1``) and ``j`` for the
along-edge / along-arc direction.  ``k`` is the axial direction in
every block.

Index counts
------------

- Center block:    ``(nc, nc, nk, 3)``
- Each sub-block:  ``(ni, nc, nk, 3)`` — ``nc`` is shared with the
  central block to guarantee SHARED conformity at every seam.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.spacing import Spacing1D
from shore.primitive.base import PrimitiveMesh, axis_frame, coords_1d, resolve_axis
from shore.volume.spacing import tanh_one_sided, tanh_two_sided


@dataclass(frozen=True)
class FlatCapsExtent:
    """Outer radius and axial extent of a flat-caps cylinder."""

    r_out: float
    z0: float
    z1: float

    def __post_init__(self) -> None:
        if self.r_out <= 0.0:
            raise ParameterError(f"FlatCapsExtent: r_out must be > 0, got {self.r_out}")
        if self.z1 <= self.z0:
            raise ParameterError(f"FlatCapsExtent: z1 ({self.z1}) must exceed z0 ({self.z0})")


def _pinned_distribution(
    n: int,
    start: float,
    stop: float,
    first_step: float,
    law: str,
    beta: float,
) -> np.ndarray:
    """Return *n* coords on ``[start, stop]`` with the first step pinned.

    ``out[1] - out[0] == first_step`` exactly; the remaining ``n - 1``
    coords lie on ``[start + first_step, stop]`` distributed by *law*.

    Falls back to the unpinned ``coords_1d`` distribution and emits a
    ``UserWarning`` when the pin would produce a non-monotone result
    (``first_step <= 0``, ``first_step >= stop - start``, or
    ``first_step`` more than 4x the natural first step of *law*).
    """
    if n < 2:
        raise ValueError(f"_pinned_distribution: n must be >= 2, got {n}")

    total = stop - start
    natural_first = total / (n - 1)
    incompatible = (
        first_step <= 0.0
        or first_step >= total
        or (first_step / max(natural_first, 1e-30) > 4.0)
        or (natural_first / max(first_step, 1e-30) > 4.0)
    )
    if incompatible:
        warnings.warn(
            f"FlatCapsCylinder: cannot pin sub-block first radial step "
            f"({first_step:.4g}) — incompatible with i_law={law!r}, "
            f"ni={n}, radial extent={total:.4g} (natural first step "
            f"~{natural_first:.4g}); falling back to unpinned distribution.  "
            f"C1 match at the central / sub-block seam is approximate.  "
            f"Adjust ni, nc, or inner_frac so the pinned step fits.",
            UserWarning,
            stacklevel=2,
        )
        return coords_1d(n, start, stop, law=law, beta=beta)

    out = np.empty(n, dtype=np.float64)
    out[0] = start
    out[1] = start + first_step
    if n == 2:
        return out

    inner_start = start + first_step
    if law == "uniform":
        out[1:] = np.linspace(inner_start, stop, n - 1)
    elif law == "tanh":
        out[1:] = tanh_one_sided(n - 1, inner_start, stop, beta)
    elif law == "tanh2":
        out[1:] = tanh_two_sided(n - 1, inner_start, stop, beta)
    else:
        raise ValueError(f"_pinned_distribution: unknown law {law!r}")
    return out


class FlatCapsCylinder(PrimitiveMesh):
    """5-block flat-caps cylinder primitive.

    Examples
    --------
    >>> fc = FlatCapsCylinder.build(r_out=2.0, z0=0.0, z1=4.0,
    ...                             ni=10, nc=12, nk=20)
    >>> [b.label for b in fc.blocks]
    ['center', 'sub_e', 'sub_n', 'sub_w', 'sub_s']
    """

    label_default = "flat_caps"

    def __init__(self, blocks: list[HexBlock]) -> None:
        self.blocks = blocks

    @classmethod
    def build(
        cls,
        r_out: float,
        z0: float,
        z1: float,
        ni: int,
        nc: int,
        nk: int,
        inner_frac: float = 0.5,
        axis: str | tuple[float, float, float] = "z",
        center: tuple[float, float, float] = (0.0, 0.0, 0.0),
        i_law: str = "uniform",
        c_law: str = "uniform",
        k_law: str = "uniform",
        i_beta: float = 3.0,
        c_beta: float = 3.0,
        k_beta: float = 3.0,
        edge_pins: dict[str, Spacing1D] | None = None,
    ) -> FlatCapsCylinder:
        """Build a flat-caps cylinder primitive.

        Parameters
        ----------
        r_out :
            Outer radius of the disk.
        z0, z1 :
            Axial extent along *axis* (``z1 > z0``).
        ni :
            Radial node count for each sub-block (square edge → outer
            circle).  Must be >= 2.
        nc :
            Shared node count along each side of the central square.
            Each sub-block uses *nc* nodes along the arc; the central
            block is ``(nc, nc)`` in (e1, e2).  Must be >= 2.
        nk :
            Axial node count.  Must be >= 2.
        inner_frac :
            Ratio ``r_inner / r_out`` of the central square's
            corner-radius to the outer disk radius.  Must lie in
            ``(0, 1)``.  Default 0.5.
        axis :
            Cylinder axis (``"x"``, ``"y"``, ``"z"``, or a 3-vector).
        center :
            World-space anchor for the cylinder axis.
        i_law, c_law, k_law, *beta :
            Spacing laws.  ``i`` controls the radial sub-block
            distribution (with first cell pinned to the central first
            cell); ``c`` controls the central-square / sub-block-arc
            distribution (shared between center and sub-blocks);
            ``k`` controls the axial distribution.
        edge_pins :
            Optional sparse override mapping k-axis edge name →
            ``Spacing1D``.  When any k-edge is pinned the axial
            coordinate is built by transfinite interpolation between
            the 4 parallel k-edges of every block, allowing different
            axial spacing on (e.g.) the inner radius vs the outer
            circumference at the same k-layer.

            Pins on i- or c-axis edges are **rejected**: the radial
            (i) and tangential (c) distributions on this primitive are
            shaped by the geometric construction (per-j first-step
            matching at the central-square seam) and overriding them
            would require rewriting that geometry.

        Returns
        -------
        FlatCapsCylinder
        """
        extent = FlatCapsExtent(r_out=r_out, z0=z0, z1=z1)
        for n, name in ((ni, "ni"), (nc, "nc"), (nk, "nk")):
            if n < 2:
                raise ParameterError(f"FlatCapsCylinder: {name} must be >= 2, got {n}")
        if not (0.0 < inner_frac < 1.0):
            raise ParameterError(
                f"FlatCapsCylinder: inner_frac must be in (0, 1), got {inner_frac}"
            )

        # FlatCapsCylinder pins all three axes geometrically: i (radial)
        # is per-j-pinned to the central square's tangential step,
        # c (tangential) is set by the central square / outer circle
        # geometry, and k is uniform across the 5 blocks because radial
        # SHARED seams must coincide layer by layer.  Per-edge spacing
        # therefore has no degree of freedom on this primitive.
        # ``edge_pins`` is accepted for API symmetry with Box / Annulus
        # but raises if any entry is supplied.
        if edge_pins:
            raise ParameterError(
                "FlatCapsCylinder: edge_pins is not supported on this "
                "primitive.  The radial (i), tangential (c) and axial (k) "
                "distributions are pinned by the geometric construction "
                "and the radial seams across the 5 blocks must remain "
                "coincident at every k-layer.  Use i_law / c_law / k_law "
                "for global control."
            )

        axis_vec = resolve_axis(axis)
        e1, e2, e3 = axis_frame(axis_vec)
        c0 = np.asarray(center, dtype=float)

        # Central-square geometry.
        # Corners sit at radius r_inner at azimuths pi/4 + k*pi/2.
        # The two perpendicular sides therefore have half-length
        # ``half = r_inner * cos(pi/4) = r_inner / sqrt(2)``.
        r_inner = inner_frac * extent.r_out
        half = r_inner * np.cos(np.pi / 4.0)
        s = 2.0 * half  # full side length

        # Shared 1D distributions.
        # ``u_c`` is the spatial coordinate along a side of the central
        # square (in [-half, +half]).  ``u_eta01`` is the same
        # distribution remapped to [0, 1] — used as the parametric
        # coordinate along each sub-block's inner edge / outer arc.
        u_c = coords_1d(nc, -half, half, law=c_law, beta=c_beta)
        u_eta01 = (u_c - (-half)) / s  # exact remap, preserves the law's spacing
        zs = coords_1d(nk, extent.z0, extent.z1, law=k_law, beta=k_beta)

        # Central-square along-side cell sizes.  ``c_first[j]`` is the
        # local tangential step in the central square at column j of a
        # given side — used to pin the sub-block's first radial step
        # *per-j*, so seam cells have radial == tangential size on
        # both sides of the central / sub-block seam.
        c_diff = np.diff(u_c)  # (nc - 1,)
        c_first = np.empty(nc, dtype=np.float64)
        c_first[0] = c_diff[0]
        c_first[-1] = c_diff[-1]
        if nc >= 3:
            c_first[1:-1] = 0.5 * (c_diff[:-1] + c_diff[1:])  # centred

        # ------------------------------------------------------------------
        # Center block: (nc, nc, nk, 3) — Cartesian in (e1, e2).
        # ------------------------------------------------------------------
        U, V = np.meshgrid(u_c, u_c, indexing="ij")  # (nc, nc)
        center_planar = U[..., None] * e1 + V[..., None] * e2  # (nc, nc, 3)
        Z = zs[None, None, :, None]
        center_nodes = center_planar[:, :, None, :] + Z * e3 + c0
        center_block = HexBlock(
            label="center",
            nodes=center_nodes,
            periodic_i=False,
            periodic_j=False,
        )

        # ------------------------------------------------------------------
        # Sub-blocks: 4 trapezoidal slabs around the central square.
        #
        # Each j-column of the sub-block is a *straight radial line* from
        # ``inner[j]`` (on the central square's side) to ``outer[j]``
        # (on the outer circle at the matching azimuth).  The radial
        # distribution along the line is per-j: the first step is
        # pinned to ``c_first[j]`` (the local tangential step of the
        # central square at column j) so seam cells have aspect ratio 1
        # on both sides of the seam.  The remaining ``ni - 2`` cells
        # absorb the difference, distributed by ``i_law``.
        #
        # The choice of azimuth at column j is consistent with the
        # central square: ``theta(j) = atan2(outer[j].y, outer[j].x)``
        # where ``outer[j]`` projects ``inner[j]`` from the origin onto
        # the outer circle.  This guarantees that the radial line at
        # column j is exactly aligned with the line from the origin
        # through ``inner[j]`` (so the "radial direction" deserves its
        # name).
        # ------------------------------------------------------------------
        sub_specs = [
            # (label, corner_a_xy, corner_b_xy)
            ("sub_e", (+half, -half), (+half, +half)),
            ("sub_n", (+half, +half), (-half, +half)),
            ("sub_w", (-half, +half), (-half, -half)),
            ("sub_s", (-half, -half), (+half, -half)),
        ]

        sub_blocks: list[HexBlock] = []
        warned_pin = [False]  # one warning per build, not per-sub-block by per-j

        for label, ca, cb in sub_specs:
            ca_planar = ca[0] * e1 + ca[1] * e2  # (3,)
            cb_planar = cb[0] * e1 + cb[1] * e2  # (3,)

            # Inner edge: straight line ca → cb along the central-square
            # side (eta in [0, 1] under the user's c_law).
            inner = (1.0 - u_eta01)[:, None] * ca_planar + u_eta01[:, None] * cb_planar
            # inner: (nc, 3)

            # Inner-edge radius at column j.  Equals ``r_inner`` only
            # at the corners; at the mid-side it equals ``half`` (the
            # square's apothem).  The outer projection runs along the
            # unit vector ``inner[j] / r_inner_j[j]`` so the j-column
            # of the sub-block is a straight radial line through the
            # origin.
            r_inner_j = np.linalg.norm(inner, axis=-1)  # (nc,)
            unit_radial = inner / r_inner_j[:, None]  # (nc, 3)

            # Per-j radial extent and first-step pin.
            radial_extent = extent.r_out - r_inner_j  # (nc,)

            # Build the (ni, nc) radius-along-radial-line distribution.
            # node[i, j] = inner[j] + r_local[i, j] * unit_radial[j]
            # where r_local[i, j] in [r_inner_j[j], r_out].
            r_local = np.empty((ni, nc), dtype=np.float64)
            for j in range(nc):
                first_step = float(c_first[j])
                extent_j = float(radial_extent[j])
                if extent_j <= 0.0:
                    raise ParameterError(
                        f"FlatCapsCylinder: degenerate radial extent at j={j} "
                        f"of sub-block {label!r} (inner_frac={inner_frac} "
                        f"too large?)"
                    )
                # Use coords_1d directly for j=0 / nc-1 corners (where
                # the outer azimuth coincides with the corner azimuth
                # by construction; pinning works as expected) and for
                # the rest.  The pinned distribution lives in the
                # absolute-radius space [r_inner_j[j], r_out].
                dist_j = _pinned_distribution(
                    n=ni,
                    start=float(r_inner_j[j]),
                    stop=extent.r_out,
                    first_step=first_step,
                    law=i_law,
                    beta=i_beta,
                )
                r_local[:, j] = dist_j
                # _pinned_distribution emits its own warning per call;
                # collapse repeated identical warnings to a single one
                # per build to avoid spamming.
                if first_step >= extent_j or first_step <= 0.0:
                    warned_pin[0] = True

            # Build sub-block planar nodes:
            #   nodes[i, j] = r_local[i, j] * unit_radial[j]
            # which lives at radius r_local[i, j] along the j-th radial
            # line.  This automatically gives ``nodes[0, j] == inner[j]``
            # and ``nodes[-1, j] == outer[j]``.
            nodes_planar = r_local[:, :, None] * unit_radial[None, :, :]
            # (ni, nc, 3)

            # Extrude axially.
            nodes = nodes_planar[:, :, None, :] + Z * e3 + c0
            # (ni, nc, nk, 3)

            sub_blocks.append(
                HexBlock(
                    label=label,
                    nodes=nodes,
                    periodic_i=False,
                    periodic_j=False,
                )
            )

        return cls(blocks=[center_block, *sub_blocks])
