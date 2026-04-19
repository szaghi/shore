"""Topology ABC and CapsTopology: block wiring from a k=0 surface."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from shore.mesh.block import HexBlock
from shore.mesh.face import FaceBC
from shore.mesh.spacing import Spacing1D
from shore.topology.sphere_cap import (
    _adjust_nj,
    _build_cap_k0,
    _pull_cap_interior,
    _squareness_rows,
)


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


class CapsTopology(Topology):
    """6-block O-grid with polar caps.

    Layout: ``[sub0, sub1, sub2, sub3, cap_north, cap_south]``

    Shared faces are numpy views — writing to one block's boundary
    automatically updates the adjacent block.

    Parameters
    ----------
    cap_pull :
        Dome factor for cap interior (0 = flat Laplacian).
    cap_squareness :
        Straighten cap boundary edges toward geodesic chord (0 = circular arc).
    """

    def __init__(
        self,
        cap_pull: float = 0.0,
        cap_squareness: float = 0.0,
    ) -> None:
        self.cap_pull       = cap_pull
        self.cap_squareness = cap_squareness

    def build(
        self,
        surface: np.ndarray,
        nk: int,
        nc_cap: int | None = None,
        spacing_k: Spacing1D | None = None,
        mesh: object | None = None,
        **kwargs,
    ) -> list[HexBlock]:
        """Build 6 blocks from the k=0 surface layer.

        Parameters
        ----------
        surface :
            Projected k=0 surface.  When ``nc_cap`` is given, shape is
            ``(ni_total, nj, 3)`` covering the full sphere; otherwise
            ``(ni, nj, 3)`` covering only the equatorial band (legacy path,
            caps built by Laplacian).
        nk :
            Total number of wall-normal layers (including k=0).
        nc_cap :
            Number of rows in each polar cap block as returned by
            ``Body.project_caps`` / ``project_caps_to_surface``.  When
            supplied the cap k=0 interior is ray-cast from the STL (requires
            ``mesh``); otherwise falls back to Laplacian.
        spacing_k :
            Wall-normal spacing law shared by all blocks.  Defaults to
            ``Spacing1D(law='geometric', ds=1e-3)``.
        mesh :
            ``trimesh.Trimesh`` instance.  Required for the ray-cast cap
            interior path (when ``nc_cap`` is supplied).

        Returns
        -------
        list[HexBlock]
            ``[sub0, sub1, sub2, sub3, cap_north, cap_south]``
        """
        ni_total_raw, nj_raw, _ = surface.shape
        nj = _adjust_nj(nj_raw)
        nc_j = nj // 4 + 1   # j-size of every equatorial sub-block and cap edge
        sector = nj // 4

        if nj < nj_raw:
            surface = surface[:, :nj, :]

        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

        # ------------------------------------------------------------------
        # Split full surface into equatorial band
        # When nc_cap is supplied the surface covers the full sphere; the
        # equatorial band is the middle ni rows (shared boundary rows with caps).
        # ------------------------------------------------------------------
        if nc_cap is not None:
            ni = ni_total_raw - 2 * (nc_cap - 1)
            surf_eq = surface[nc_cap - 1 : nc_cap - 1 + ni]   # (ni, nj, 3)
        else:
            surf_eq = surface
            ni = ni_total_raw

        # ------------------------------------------------------------------
        # Equatorial ring: one contiguous allocation, 4 views
        # ring shape: (ni, 4*(nc_j-1)+1, nk, 3)  — full periodic j-ring
        # sub_q gets columns [q*(nc_j-1) : q*(nc_j-1)+nc_j]
        # ------------------------------------------------------------------
        ring_nj = 4 * (nc_j - 1) + 1
        ring = np.empty((ni, ring_nj, nk, 3), dtype=np.float64)

        sub_nodes: list[np.ndarray] = []
        for q in range(4):
            j0 = q * (nc_j - 1)
            sub_nodes.append(ring[:, j0 : j0 + nc_j, :, :])

        # Fill k=0 from equatorial surface
        for q in range(4):
            j0 = q * sector
            j_indices = [*range(j0, j0 + sector), (j0 + sector) % nj]
            sub_nodes[q][:, :, 0, :] = surf_eq[:, j_indices, :]

        ring[:, -1, 0, :] = ring[:, 0, 0, :]

        # ------------------------------------------------------------------
        # Cap k=0 seeds
        # Sub-block pole rows are the cap boundary.  For cap_squareness > 0
        # the boundary rows fed to the cap seed are squareness-modified copies
        # (so the cap k=0 is consistent with what marching will produce at k=1),
        # but the sub_nodes ring is NOT modified — equatorial nodes stay on STL.
        # ------------------------------------------------------------------
        north_eq_rows_k0 = [sub_nodes[q][ni - 1, :, 0, :] for q in range(4)]
        south_eq_rows_k0 = [sub_nodes[q][0, :, 0, :] for q in range(4)]

        if self.cap_squareness != 0.0:
            north_cap_rows_k0 = _squareness_rows(north_eq_rows_k0, self.cap_squareness)
            south_cap_rows_k0 = _squareness_rows(south_eq_rows_k0, self.cap_squareness)
            # Re-stitch shared corners on the cap-seed rows (copies only)
            for q in range(4):
                q_next = (q + 1) % 4
                shared_n = 0.5 * (north_cap_rows_k0[q][-1] + north_cap_rows_k0[q_next][0])
                north_cap_rows_k0[q][-1]       = shared_n
                north_cap_rows_k0[q_next][0]   = shared_n
                shared_s = 0.5 * (south_cap_rows_k0[q][-1] + south_cap_rows_k0[q_next][0])
                south_cap_rows_k0[q][-1]       = shared_s
                south_cap_rows_k0[q_next][0]   = shared_s
        else:
            north_cap_rows_k0 = north_eq_rows_k0
            south_cap_rows_k0 = south_eq_rows_k0

        cap_north_nodes = np.empty((nc_j, nc_j, nk, 3), dtype=np.float64)
        cap_south_nodes = np.empty((nc_j, nc_j, nk, 3), dtype=np.float64)

        if mesh is not None:
            from shore.surface.remesh import project_gnomonic_cap_to_surface
            cap_north_k0 = project_gnomonic_cap_to_surface(
                mesh, north_cap_rows_k0, nc_j, pole="north"
            )
            cap_south_k0 = project_gnomonic_cap_to_surface(
                mesh, south_cap_rows_k0, nc_j, pole="south"
            )
        else:
            cap_north_k0 = _build_cap_k0(north_cap_rows_k0, nc_j, pole="north")
            cap_south_k0 = _build_cap_k0(south_cap_rows_k0, nc_j, pole="south")

            if self.cap_pull != 0.0:
                cap_north_k0 = _pull_cap_interior(
                    cap_north_k0, north_cap_rows_k0, self.cap_pull, pole="north"
                )
                cap_south_k0 = _pull_cap_interior(
                    cap_south_k0, south_cap_rows_k0, self.cap_pull, pole="south"
                )

        cap_north_nodes[:, :, 0, :] = cap_north_k0
        cap_south_nodes[:, :, 0, :] = cap_south_k0

        # ------------------------------------------------------------------
        # Construct HexBlocks
        # ------------------------------------------------------------------
        sub_blocks: list[HexBlock] = []
        for q in range(4):
            blk = HexBlock(
                label=f"sub{q}",
                nodes=sub_nodes[q],
                periodic_i=False,
                periodic_j=False,
                spacing_k=_sp_k,
            )
            sub_blocks.append(blk)

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
        # Sub-block adjacency: sub_q.j_hi ↔ sub_{q+1 mod 4}.j_lo
        # These are already numpy views into the same ring memory — mark as SHARED.
        for q in range(4):
            q_next = (q + 1) % 4
            f_hi   = sub_blocks[q].face("j_hi")
            f_lo   = sub_blocks[q_next].face("j_lo")
            sub_blocks[q]._set_face_bc("j_hi", FaceBC.SHARED, partner=f_lo)
            sub_blocks[q_next]._set_face_bc("j_lo", FaceBC.SHARED, partner=f_hi)

        # Cap boundaries: DIRICHLET from equatorial pole rows
        for blk in sub_blocks:
            blk._set_face_bc("i_hi", FaceBC.DIRICHLET,
                             partner=cap_north_blk.face("k_lo"))
            blk._set_face_bc("i_lo", FaceBC.DIRICHLET,
                             partner=cap_south_blk.face("k_lo"))

        for name in ("i_lo", "i_hi", "j_lo", "j_hi"):
            cap_north_blk._set_face_bc(name, FaceBC.DIRICHLET)
            cap_south_blk._set_face_bc(name, FaceBC.DIRICHLET)

        return [*sub_blocks, cap_north_blk, cap_south_blk]
