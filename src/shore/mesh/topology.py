"""Topology ABC and CubedSphereTopology: block wiring from a k=0 surface."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC
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
                "CubedSphereTopology.build requires a Projector "
                "(pass projector= or mesh=)."
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
        for q in range(4):
            q_next = (q + 1) % 4
            f_hi = sub_blocks[q].face("j_hi")
            f_lo = sub_blocks[q_next].face("j_lo")
            sub_blocks[q]._set_face_bc("j_hi", FaceBC.SHARED, partner=f_lo)
            sub_blocks[q_next]._set_face_bc("j_lo", FaceBC.SHARED, partner=f_hi)

        # Cap boundaries: DIRICHLET from equatorial pole rows.
        for blk in sub_blocks:
            blk._set_face_bc("i_hi", FaceBC.DIRICHLET, partner=cap_north_blk.face("k_lo"))
            blk._set_face_bc("i_lo", FaceBC.DIRICHLET, partner=cap_south_blk.face("k_lo"))

        for name in ("i_lo", "i_hi", "j_lo", "j_hi"):
            cap_north_blk._set_face_bc(name, FaceBC.DIRICHLET)
            cap_south_blk._set_face_bc(name, FaceBC.DIRICHLET)

        return [*sub_blocks, cap_north_blk, cap_south_blk]


def _projector_from_mesh(mesh: object) -> object:
    """Build a temporary Projector using mesh.centroid (legacy mesh= path)."""
    from shore.surface.projector import Projector

    return Projector(mesh, mesh.centroid)  # type: ignore[attr-defined]
