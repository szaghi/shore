"""Topology ABC and concrete topologies (cubed-sphere, O-grid)."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from shore.mesh.block import HexBlock
from shore.mesh.face import AxisMap, FaceBC, FaceName
from shore.mesh.spacing import Spacing1D
from shore.surface.cap import build_cap_k0
from shore.surface.equator import build_equator
from shore.topology.sphere_cap import _adjust_nj


class Topology(ABC):
    """Abstract topology builder.

    Receives an already-projected k=0 surface layer and returns a list of
    initialised blocks with shared-face views wired.  Marching is not
    performed here.
    """

    @abstractmethod
    def build(
        self,
        surface: np.ndarray,
        nk: int,
        **kwargs,
    ) -> list[HexBlock]:
        """Build blocks from the k=0 surface layer.

        Parameters
        ----------
        surface :
            Projected k=0 surface, shape ``(ni, nj, 3)``.
        nk :
            Total number of wall-normal layers (including k=0).

        Returns
        -------
        list[HexBlock]
            Initialised blocks with k=0 filled and shared views wired.
            Volume layers k=1..nk-1 are allocated but not filled.
        """


# Cubed-sphere cap-equator seam tables.  Each entry pairs a sub-block
# index q (0..3) with (cap_face_name, sub_side_axis_map).  Derived from
# ``shore.topology.sphere_cap._overwrite_cap_boundaries``: the cap layer
# k=0 row equals the sub-block's i-extreme row, optionally with j-axis
# reversal.  AxisMaps live on the *sub* side of the seam (sub.i_hi for
# north pole, sub.i_lo for south); the cap-side AxisMap is the inverse.
_CAP_EQUATOR_NORTH: tuple[tuple[FaceName, AxisMap], ...] = (
    # q=0: cap_north[0, :, :]  = sub0[ni-1, :, :]      (no flip)
    ("i_lo", AxisMap.from_pairs(("j", False), ("k", False))),
    # q=1: cap_north[:, -1, :] = sub1[ni-1, :, :]      (cross-axis, no flip)
    ("j_hi", AxisMap.from_pairs(("i", False), ("k", False))),
    # q=2: cap_north[-1, :, :] = sub2[ni-1, ::-1, :]   (j-flip)
    ("i_hi", AxisMap.from_pairs(("j", True), ("k", False))),
    # q=3: cap_north[:, 0, :]  = sub3[ni-1, ::-1, :]   (cross-axis, j-flip)
    ("j_lo", AxisMap.from_pairs(("i", True), ("k", False))),
)
_CAP_EQUATOR_SOUTH: tuple[tuple[FaceName, AxisMap], ...] = (
    # q=0: cap_south[:, 0, :]  = sub0[0, :, :]         (cross-axis, no flip)
    ("j_lo", AxisMap.from_pairs(("i", False), ("k", False))),
    # q=1: cap_south[-1, :, :] = sub1[0, :, :]         (no flip)
    ("i_hi", AxisMap.from_pairs(("j", False), ("k", False))),
    # q=2: cap_south[:, -1, :] = sub2[0, ::-1, :]      (cross-axis, j-flip)
    ("j_hi", AxisMap.from_pairs(("i", True), ("k", False))),
    # q=3: cap_south[0, :, :]  = sub3[0, ::-1, :]      (j-flip)
    ("i_lo", AxisMap.from_pairs(("j", True), ("k", False))),
)


class CubedSphereTopology(Topology):
    """6-block cubed-sphere topology (Sadourny 1972 / Ronchi-Iacono-Paolucci 1996).

    A six-block decomposition of the body surface: 4 equatorial sub-blocks
    plus 2 polar caps.  Each cap k=0 face is built by gnomonic flat-square
    projection (a tangent-plane square radially projected onto the body),
    matching the construction of cubed-sphere atmospheric models.

    Layout: ``[sub0, sub1, sub2, sub3, cap_north, cap_south]``

    Construction order (gnomonic algorithm):

    1. Build the south cap k=0 by gnomonic flat-square projection (4 corners
       at colatitude ``theta_cap``, axis-aligned with sub-block junctions).
    2. Build the north cap k=0 the same way.
    3. Build the equator k=0 by slerp-meridian interpolation between the two
       cap seams, ray-cast onto the STL.
    4. Wire the blocks: sub_q ↔ sub_(q+1) periodic SHARED face seams, and
       sub_q ↔ cap DIRICHLET face links at the seam.

    Shared faces use numpy views into a contiguous equator-ring allocation.
    """

    def build(
        self,
        surface: np.ndarray | None,
        nk: int,
        ni: int | None = None,
        nj: int | None = None,
        spacing_k: Spacing1D | None = None,
        mesh: object | None = None,
        projector: object | None = None,
        theta_cap: float = 0.05 * np.pi,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        **kwargs,
    ) -> list[HexBlock]:
        """Build 6 blocks from the body's STL via gnomonic cap projection.

        Parameters
        ----------
        surface :
            Ignored — preserved as the first positional for backwards-compat
            with the old CubedSphereTopology.build(surface, nk, ...) signature.  The
            new construction projects directly through the cached Projector,
            so a pre-projected surface array is not used.  When ``ni`` or
            ``nj`` are missing, they are derived from ``surface.shape``.
        nk :
            Total number of wall-normal layers (including k=0).
        ni, nj :
            Equator resolution.  ``nj`` must be a multiple of 4.  If omitted,
            both are derived from ``surface.shape``.
        spacing_k :
            Wall-normal spacing law shared by all blocks.  Defaults to
            ``Spacing1D(law='geometric', ds=1e-3)``.
        mesh :
            ``trimesh.Trimesh`` instance.  Legacy alternative to *projector*;
            creates a temporary Projector using ``mesh.centroid``.
        projector :
            Pre-built ``Projector`` instance (shares cached BVH with the
            projection step).  Preferred over *mesh* — avoids a second BVH
            build.  Required for the gnomonic algorithm; if ``None`` and no
            ``mesh`` is supplied, a ValueError is raised.
        theta_cap :
            Polar exclusion angle in radians.  The 4 cap corners sit at
            colatitude ``theta_cap`` on a sphere body; cap edges sweep
            poleward at midpoints (the natural gnomonic-square geometry).
        i_spacing :
            Latitude spacing law of the equator (``'uniform'`` | ``'tanh'`` |
            ``'tanh2'``).  Affects only the equator's meridional distribution.
        i_beta :
            Clustering strength for tanh/tanh2 i-spacing.
        j_spacing :
            Tangent-plane spacing law of the cap edges (``'uniform'`` or
            ``'tanh2'``).  Inherited by the equator seam since the equator's
            seam rows are the cap edges bit-exactly.
        j_beta :
            Clustering strength for tanh2 j-spacing.

        Returns
        -------
        list[HexBlock]
            ``[sub0, sub1, sub2, sub3, cap_north, cap_south]``
        """
        # ------------------------------------------------------------------
        # Resolve ni, nj
        # ------------------------------------------------------------------
        if ni is None or nj is None:
            if surface is None:
                raise ValueError(
                    "CubedSphereTopology.build needs ni and nj "
                    "(or a surface array to derive them from)."
                )
            ni = ni if ni is not None else int(surface.shape[0])
            nj = nj if nj is not None else int(surface.shape[1])
        nj = _adjust_nj(nj)
        nc = nj // 4 + 1

        # ------------------------------------------------------------------
        # Resolve projector
        # ------------------------------------------------------------------
        cap_proj = (
            projector
            if projector is not None
            else (_projector_from_mesh(mesh) if mesh is not None else None)
        )
        if cap_proj is None:
            raise ValueError(
                "CubedSphereTopology.build requires a Projector (pass projector= or mesh=)."
            )

        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

        # ------------------------------------------------------------------
        # Cap k=0 seeds — gnomonic flat-square projection
        # ------------------------------------------------------------------
        cap_south_k0 = build_cap_k0(
            theta_cap=theta_cap,
            nc=nc,
            pole="south",
            projector=cap_proj,
            j_spacing=j_spacing,
            j_beta=j_beta,
        )
        cap_north_k0 = build_cap_k0(
            theta_cap=theta_cap,
            nc=nc,
            pole="north",
            projector=cap_proj,
            j_spacing=j_spacing,
            j_beta=j_beta,
        )

        # ------------------------------------------------------------------
        # Equator k=0 — great-circle meridian slerp from south seam to north
        # seam, ray-cast onto the STL.
        # ------------------------------------------------------------------
        equator_k0 = build_equator(
            cap_south=cap_south_k0,
            cap_north=cap_north_k0,
            ni=ni,
            theta_cap=theta_cap,
            projector=cap_proj,
            i_spacing=i_spacing,
            i_beta=i_beta,
        )  # (ni, nj, 3)

        # ------------------------------------------------------------------
        # Equatorial ring: contiguous (ni, 4*(nc-1)+1, nk, 3) allocation
        # with 4 sub-block views sharing memory at the j-boundaries.
        # ------------------------------------------------------------------
        ring_nj = 4 * (nc - 1) + 1
        ring = np.empty((ni, ring_nj, nk, 3), dtype=np.float64)

        sector = nc - 1
        sub_nodes: list[np.ndarray] = []
        for q in range(4):
            j0 = q * sector
            sub_nodes.append(ring[:, j0 : j0 + nc, :, :])

        # Fill k=0 from the equator surface.  ``equator_k0[:, j_idx, :]`` for
        # j_idx ∈ [q*sector, q*sector+sector] gives sub_q's nc longitude
        # nodes; the last node is shared with sub_(q+1)'s first node, hence
        # the ``+ 1`` in the slice end.
        for q in range(4):
            j0 = q * sector
            j_indices = [*range(j0, j0 + sector), (j0 + sector) % nj]
            sub_nodes[q][:, :, 0, :] = equator_k0[:, j_indices, :]
        ring[:, -1, 0, :] = ring[:, 0, 0, :]

        # ------------------------------------------------------------------
        # Cap blocks
        # ------------------------------------------------------------------
        cap_north_nodes = np.empty((nc, nc, nk, 3), dtype=np.float64)
        cap_south_nodes = np.empty((nc, nc, nk, 3), dtype=np.float64)
        cap_north_nodes[:, :, 0, :] = cap_north_k0
        cap_south_nodes[:, :, 0, :] = cap_south_k0

        # ------------------------------------------------------------------
        # Construct HexBlocks
        # ------------------------------------------------------------------
        sub_blocks: list[HexBlock] = []
        for q in range(4):
            sub_blocks.append(
                HexBlock(
                    label=f"sub{q}",
                    nodes=sub_nodes[q],
                    periodic_i=False,
                    periodic_j=False,
                    spacing_k=_sp_k,
                )
            )

        cap_north_blk = HexBlock(
            label="cap_north",
            nodes=cap_north_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        cap_south_blk = HexBlock(
            label="cap_south",
            nodes=cap_south_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )

        # ------------------------------------------------------------------
        # Wire face BCs
        # ------------------------------------------------------------------
        # Sub-block adjacency: sub_q.j_hi ↔ sub_{q+1 mod 4}.j_lo (SHARED views).
        # Same-axis seam, no flips: identity AxisMap on each side.
        for q in range(4):
            q_next = (q + 1) % 4
            f_hi = sub_blocks[q].face("j_hi")
            f_lo = sub_blocks[q_next].face("j_lo")
            sub_blocks[q]._set_face_bc(
                "j_hi", FaceBC.SHARED, partner=f_lo, axis_map=AxisMap.identity("j_hi")
            )
            sub_blocks[q_next]._set_face_bc(
                "j_lo", FaceBC.SHARED, partner=f_hi, axis_map=AxisMap.identity("j_lo")
            )

        # Cap-equator seams: 8 DIRICHLET face-to-face partnerships.
        #
        # Each sub-block's i_hi (north pole) and i_lo (south pole) is a 2D
        # face that coincides node-by-node with one lateral face of the
        # corresponding cap block, possibly with j-axis reversal (per
        # ``shore.topology.sphere_cap._overwrite_cap_boundaries``).  The
        # cap-side partner is therefore the cap's lateral face — *not*
        # the cap's k_lo (which is the cap's own wall face).
        #
        # Each table entry pairs a sub-block index q with the cap-side
        # face name and the AxisMap on the *sub side* of the seam.  The
        # cap-side AxisMap is derived by ``axis_map.inverse(...)``.
        sub_face: FaceName  # narrows the literal for mypy
        for sub_face, cap_blk, table in (
            ("i_hi", cap_north_blk, _CAP_EQUATOR_NORTH),
            ("i_lo", cap_south_blk, _CAP_EQUATOR_SOUTH),
        ):
            for q, (cap_face, sub_axis_map) in enumerate(table):
                sub_blk = sub_blocks[q]
                cap_blk._set_face_bc(
                    cap_face,
                    FaceBC.DIRICHLET,
                    partner=sub_blk.face(sub_face),
                    axis_map=sub_axis_map.inverse(sub_face, cap_face),
                )
                sub_blk._set_face_bc(
                    sub_face,
                    FaceBC.DIRICHLET,
                    partner=cap_blk.face(cap_face),
                    axis_map=sub_axis_map,
                )

        return [*sub_blocks, cap_north_blk, cap_south_blk]


class OGridTopology(Topology):
    """Single-block lat-lon O-grid topology.

    A single hexahedral block with shape ``(ni, nj, nk, 3)``, j-periodic
    (the longitude direction wraps), and polar singularities at ``i = 0``
    (south pole at colatitude ``theta_cap``) and ``i = ni - 1`` (north
    pole at colatitude ``pi - theta_cap``).  This is the simpler and
    older sphere covering — useful as a Chimera background or in
    downstream tools that expect a single periodic-j block — and it
    coexists with :class:`CubedSphereTopology` as a peer, not a
    component.

    The k=0 surface is built by ray-cast from the body's anchor through
    a structured lat-lon sphere template.  The same ``i_spacing`` /
    ``j_spacing`` flags as the cubed sphere apply, with the same
    semantics: ``i`` controls the meridional distribution, ``j`` the
    longitude.

    Per-edge spacing
    ----------------
    Per-edge **k-spacing** works on O-grid via the same identity-check
    mechanism as :class:`CubedSphereTopology`: replace the four k-edges'
    ``spacing`` field with a shared ``Spacing1D`` instance, and
    :meth:`shore.mesh.Mesh.march` reads it instead of its ``ds``/``growth``
    arguments.  See :func:`shore.mesh.mesh._resolve_k_steps` for the rules.

    Per-edge **i-spacing** has no degree of freedom on O-grid: the j
    direction is periodic, so the four parallel i-edges
    (``j_lo_k_lo`` / ``j_hi_k_lo`` / ``j_lo_k_hi`` / ``j_hi_k_hi``) lie
    at *adjacent* azimuths under the j-wrap rather than at distinct
    transverse corners.  Use the per-axis ``i_spacing`` / ``i_beta``
    arguments instead.

    Per-edge **j-spacing** is meaningless on a periodic axis (the
    "first" and "last" cells of any j-edge are the same cell).
    """

    def build(
        self,
        surface: np.ndarray | None,
        nk: int,
        ni: int | None = None,
        nj: int | None = None,
        spacing_k: Spacing1D | None = None,
        mesh: object | None = None,
        projector: object | None = None,
        theta_cap: float = 0.05 * np.pi,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        **kwargs,
    ) -> list[HexBlock]:
        """Build a single-block O-grid topology from the body's STL.

        Parameters
        ----------
        surface :
            Optional pre-projected k=0 surface, shape ``(ni, nj, 3)``.
            When supplied, *ni* and *nj* are derived from the array
            and no ray-cast is performed.  When ``None``, the surface
            is built by ray-cast through *projector* (or *mesh*).
        nk :
            Total number of wall-normal layers (including k=0).
        ni, nj :
            Latitude (meridional) and longitude (azimuthal) resolution.
            Ignored when *surface* is supplied.
        spacing_k :
            Wall-normal spacing law.  Defaults to
            ``Spacing1D(law='geometric', ds=1e-3)``.
        mesh :
            ``trimesh.Trimesh`` instance.  Legacy alternative to
            *projector*.
        projector :
            Pre-built ``Projector`` instance.  Required when *surface*
            is None and *mesh* is not supplied.
        theta_cap :
            Polar exclusion angle in radians.  ``i = 0`` sits at
            colatitude *theta_cap*, ``i = ni - 1`` at
            ``pi - theta_cap``.
        i_spacing, j_spacing :
            Surface spacing laws.  ``i_spacing`` controls the
            meridional distribution; ``j_spacing`` controls the
            longitude distribution (only ``"uniform"`` and ``"tanh2"``
            are valid for j when the j direction is periodic).
        i_beta, j_beta :
            Clustering strengths for tanh / tanh2.

        Returns
        -------
        list[HexBlock]
            Single-element list containing the O-grid block.
        """
        # Resolve the k=0 surface -------------------------------------
        if surface is None:
            from shore.surface.remesh import sphere_latlon_directions

            if ni is None or nj is None:
                raise ValueError("OGridTopology.build needs ni and nj when no surface is supplied.")
            proj = (
                projector
                if projector is not None
                else (_projector_from_mesh(mesh) if mesh is not None else None)
            )
            if proj is None:
                raise ValueError(
                    "OGridTopology.build requires a Projector "
                    "(pass projector= or mesh=) when surface is None."
                )
            dirs = sphere_latlon_directions(
                ni,
                nj,
                theta_cap=theta_cap,
                i_spacing=i_spacing,
                i_beta=i_beta,
                j_spacing=j_spacing,
                j_beta=j_beta,
            )
            result = proj.project(dirs.reshape(-1, 3))  # type: ignore[union-attr]
            surface = result.locations.reshape(ni, nj, 3)
        else:
            ni = int(surface.shape[0])
            nj = int(surface.shape[1])

        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

        # Allocate (ni, nj, nk, 3) and fill k=0 ------------------------
        nodes = np.empty((ni, nj, nk, 3), dtype=np.float64)
        nodes[:, :, 0, :] = surface

        block = HexBlock(
            label="ogrid",
            nodes=nodes,
            periodic_i=False,
            periodic_j=True,
            spacing_k=_sp_k,
        )

        # No multi-block seams to wire; non-periodic i boundaries are
        # left at their default DIRICHLET BC, j is periodic.
        return [block]


def _projector_from_mesh(mesh: object) -> object:
    """Build a temporary Projector using mesh.centroid (legacy mesh= path)."""
    from shore.surface.projector import Projector

    return Projector(mesh, mesh.centroid)  # type: ignore[attr-defined]
