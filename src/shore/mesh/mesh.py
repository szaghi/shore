"""Mesh: multi-block orchestrator — projection, marching, I/O."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

from shore.mesh.block import HexBlock
from shore.mesh.body import Body
from shore.mesh.spacing import Spacing1D
from shore.mesh.topology import Topology
from shore.topology.sphere_cap import (
    _check_jacobian,
    _march_cap_layer,
    _march_one_layer,
)
from shore.volume.hyperbolic import _compute_normals
from shore.volume.spacing import build_k_steps

# ---------------------------------------------------------------------------
# Ghost-row helpers for seam-aware normals during marching
# ---------------------------------------------------------------------------
# Each block computes its outward normal from i- and j-tangents.  Without
# ghost rows the boundary tangents fall back to one-sided differences, which
# disagree with the partner block's centred difference at the same seam.  The
# disagreement accumulates over many marching layers as a visible spacing
# jump across block seams (sub-sub j-seams and cap-equator i-seams).
#
# These helpers gather the partner block's "row one inward" data, mapped
# through the index permutation/reversal of ``_overwrite_cap_boundaries``,
# so that ``_compute_normals`` can use centred differences at every seam.


def _sub_ghost_rows(
    q: int,
    sub_k: list[np.ndarray],
    cap_north_k: np.ndarray,
    cap_south_k: np.ndarray,
) -> dict[str, np.ndarray]:
    """Ghost rows for sub_q at the current k-layer.

    j_lo / j_hi come from the adjacent sub-blocks' row one inward from the
    shared seam.  i_lo / i_hi come from the cap blocks; the index mapping
    mirrors ``_overwrite_cap_boundaries``: for each sub_q the cap edge that
    receives sub_q's pole row has a known orientation, and the cap "row one
    inward" is the row offset by 1 in the perpendicular direction.
    """
    q_prev = (q - 1) % 4
    q_next = (q + 1) % 4
    g: dict[str, np.ndarray] = {}
    g["j_lo"] = sub_k[q_prev][:, -2, :]
    g["j_hi"] = sub_k[q_next][:, 1, :]

    # i_hi = cap_north's "row one inward from sub_q's pole edge", indexed
    # through the inverse of _overwrite_cap_boundaries (north pole).
    if q == 0:
        g["i_hi"] = cap_north_k[1, :, :]
    elif q == 1:
        g["i_hi"] = cap_north_k[:, -2, :]
    elif q == 2:
        g["i_hi"] = cap_north_k[-2, ::-1, :]
    else:  # q == 3:  cap_north[:, 0, :] ← sub3[ni-1, ::-1, :]  ⇒  ghost = cap_north[::-1, 1, :]
        g["i_hi"] = cap_north_k[::-1, 1, :]

    # i_lo = cap_south's "row one inward from sub_q's south pole edge",
    # indexed through the inverse of _overwrite_cap_boundaries (south pole).
    if q == 0:
        g["i_lo"] = cap_south_k[:, 1, :]
    elif q == 1:
        g["i_lo"] = cap_south_k[-2, :, :]
    elif q == 2:  # cap_south[:, -1, :] ← sub2[0, ::-1, :]  ⇒  ghost = cap_south[::-1, -2, :]
        g["i_lo"] = cap_south_k[::-1, -2, :]
    else:  # q == 3
        g["i_lo"] = cap_south_k[1, ::-1, :]

    return g


def _cap_ghost_rows(
    pole: str,
    sub_k: list[np.ndarray],
    ni: int,
) -> dict[str, np.ndarray]:
    """Ghost rows for a cap block (north or south) at the current k-layer.

    Mirrors the inverse of ``_overwrite_cap_boundaries``: each cap edge is
    bound to a sub-block pole row, and the ghost row beyond that edge is
    that sub-block's "row one inward from the pole" — i.e. ``[ni-2, :]`` for
    north cap, ``[1, :]`` for south cap, with the same reversal pattern as
    the override.
    """
    g: dict[str, np.ndarray] = {}
    if pole == "north":
        # cap_north[0,  :] ← sub0[ni-1, :]   → ghost beyond = sub0[ni-2, :]
        # cap_north[:, -1] ← sub1[ni-1, :]   → ghost beyond = sub1[ni-2, :]
        # cap_north[-1, :] ← sub2[ni-1, ::-1] → ghost beyond = sub2[ni-2, ::-1]
        # cap_north[:,  0] ← sub3[ni-1, ::-1] → ghost beyond = sub3[ni-2, ::-1]
        g["i_lo"] = sub_k[0][ni - 2, :, :]
        g["j_hi"] = sub_k[1][ni - 2, :, :]
        g["i_hi"] = sub_k[2][ni - 2, ::-1, :]
        g["j_lo"] = sub_k[3][ni - 2, ::-1, :]
    else:  # south
        # cap_south[:,  0] ← sub0[0, :]      → ghost beyond = sub0[1, :]
        # cap_south[-1, :] ← sub1[0, :]      → ghost beyond = sub1[1, :]
        # cap_south[:, -1] ← sub2[0, ::-1]   → ghost beyond = sub2[1, ::-1]
        # cap_south[0,  :] ← sub3[0, ::-1]   → ghost beyond = sub3[1, ::-1]
        g["j_lo"] = sub_k[0][1, :, :]
        g["i_hi"] = sub_k[1][1, :, :]
        g["j_hi"] = sub_k[2][1, ::-1, :]
        g["i_lo"] = sub_k[3][1, ::-1, :]
    return g


class Mesh:
    """Multi-block structured mesh.

    Attributes
    ----------
    body :
        The triangulated body surface.
    blocks :
        Ordered list of blocks as returned by the topology builder.
    topology :
        The topology object that built the blocks.
    """

    def __init__(
        self,
        body: Body,
        blocks: list[HexBlock],
        topology: Topology,
    ) -> None:
        self.body = body
        self.blocks = blocks
        self.topology = topology

    @classmethod
    def from_stl(
        cls,
        path: str | Path,
        topology: Topology,
        ni: int,
        nj: int,
        nk: int,
        theta_cap: float | None = None,
        theta_cap_deg: float | None = None,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        j_cluster_lon: float | None = None,
        spacing_k: Spacing1D | None = None,
    ) -> Mesh:
        """Load an STL, project to surface, wire blocks.  Does not march.

        Parameters
        ----------
        path :
            Path to the STL file.
        topology :
            Wired topology builder (e.g. ``CapsTopology()``).
        ni, nj :
            Surface resolution.
        nk :
            Number of wall-normal layers (including k=0).
        theta_cap / theta_cap_deg :
            Polar exclusion angle.
        spacing_k :
            Wall-normal spacing law.  Defaults to
            ``Spacing1D(law='geometric', ds=1e-3)``.

        Returns
        -------
        Mesh
            Ready to call ``.march()``.
        """
        body = Body.from_stl(path)

        if theta_cap_deg is not None and theta_cap is None:
            _theta = float(np.deg2rad(theta_cap_deg))
        elif theta_cap is not None:
            _theta = float(theta_cap)
        else:
            _theta = float(0.05 * np.pi)

        from shore.mesh.topology import CapsTopology

        if isinstance(topology, CapsTopology):
            # Gnomonic algorithm builds the surface inside topology.build by
            # projecting through the cached BVH; we only need ni, nj here.
            blocks = topology.build(
                surface=None,
                nk=nk,
                ni=ni,
                nj=nj,
                spacing_k=spacing_k,
                projector=body.projector,
                theta_cap=_theta,
                i_spacing=i_spacing,
                i_beta=i_beta,
                j_spacing=j_spacing,
                j_beta=j_beta,
            )
        else:
            surface = body.project(
                ni,
                nj,
                theta_cap=_theta,
                i_spacing=i_spacing,
                i_beta=i_beta,
                j_spacing=j_spacing,
                j_beta=j_beta,
                j_cluster_lon=j_cluster_lon,
            )
            blocks = topology.build(surface, nk, spacing_k=spacing_k)
        return cls(body=body, blocks=blocks, topology=topology)

    # ------------------------------------------------------------------
    # Marching
    # ------------------------------------------------------------------

    def march(
        self,
        ds: float = 1e-3,
        growth: float = 1.15,
        smooth_strength: float = 0.2,
        smooth_iters: int = 2,
        steps: np.ndarray | None = None,
        layer_callback: Callable[[int, float], None] | None = None,
        blend_normals_k: int = 0,
    ) -> None:
        """Extrude all blocks simultaneously by hyperbolic marching.

        Fills ``block.nodes[:, :, 1:, :]`` in place for every block.

        Parameters
        ----------
        ds :
            First-layer thickness (geometric law).
        growth :
            Geometric growth ratio (geometric law).
        smooth_strength :
            Laplacian blend factor per layer.
        smooth_iters :
            Number of Laplacian sweeps per layer.
        steps :
            Pre-computed step sizes ``(nk-1,)``.  Overrides ``ds``/``growth``.
        layer_callback :
            Optional ``(k, jmin_global)`` called after each accepted layer.
        blend_normals_k :
            Number of layers over which to blend STL surface normals into the
            grid-computed normals.  0 (default) disables blending.  At layer k
            the blend weight is ``alpha = 1 - k / blend_normals_k`` so k=0
            uses pure STL normals and the blend reaches zero at
            k=blend_normals_k.  Requires ``self.body.mesh``.

        Raises
        ------
        MarchingError
            If a negative-Jacobian cell is detected.
        NotImplementedError
            If the topology is not ``CapsTopology``.
        """
        from shore.mesh.topology import CapsTopology

        if not isinstance(self.topology, CapsTopology):
            raise NotImplementedError("Mesh.march() currently only supports CapsTopology.")
        self._march_caps(
            ds=ds,
            growth=growth,
            smooth_strength=smooth_strength,
            smooth_iters=smooth_iters,
            steps=steps,
            layer_callback=layer_callback,
            blend_normals_k=blend_normals_k,
        )

    def _march_caps(
        self,
        ds: float,
        growth: float,
        smooth_strength: float,
        smooth_iters: int,
        steps: np.ndarray | None,
        layer_callback: Callable[[int, float], None] | None,
        blend_normals_k: int = 0,
    ) -> None:
        """Simultaneous marching for the 6-block CapsTopology."""
        sub_blocks = self.blocks[:4]
        cap_north = self.blocks[4]
        cap_south = self.blocks[5]

        nk = sub_blocks[0].nodes.shape[2]
        ni = sub_blocks[0].nodes.shape[0]

        _steps = steps if steps is not None else build_k_steps(nk, ds, growth)

        do_blend = blend_normals_k > 0
        if do_blend:
            from shore.volume.hyperbolic import _stl_normals

            trimesh_obj = self.body.mesh

        for k in range(nk - 1):
            step = float(_steps[k])

            # Snapshot of the k-layer for every block — used as the source
            # of ghost rows so that all 6 blocks see a *consistent* set of
            # neighbour data for centred-difference normals at the seams.
            sub_k = [blk.nodes[:, :, k, :] for blk in sub_blocks]
            cn_k = cap_north.nodes[:, :, k, :]
            cs_k = cap_south.nodes[:, :, k, :]

            # --- 1. Advance equatorial sub-blocks ---
            new_subs: list[np.ndarray] = []
            for q in range(len(sub_blocks)):
                layer_q = sub_k[q]
                ni_q, nj_q = layer_q.shape[:2]

                # Ghost rows for sub_q at this layer (pre-march).  i_lo and
                # i_hi come from the cap blocks (rotated/reversed per the
                # _overwrite_cap_boundaries mapping); j_lo and j_hi come from
                # adjacent sub-blocks via the periodic ring.
                ghost = _sub_ghost_rows(q, sub_k, cn_k, cs_k)

                if do_blend and k < blend_normals_k:
                    alpha = 1.0 - k / blend_normals_k
                    grid_n = _compute_normals(
                        layer_q,
                        periodic=(False, False),
                        ghost_rows=ghost,
                    )
                    stl_n = _stl_normals(trimesh_obj, layer_q.reshape(-1, 3)).reshape(ni_q, nj_q, 3)
                    blended = alpha * stl_n + (1.0 - alpha) * grid_n
                    mag = np.linalg.norm(blended, axis=-1, keepdims=True)
                    normals_override: np.ndarray | None = blended / np.where(mag < 1e-15, 1.0, mag)
                else:
                    normals_override = None
                new_q = _march_one_layer(
                    layer_q,
                    step,
                    smooth_strength,
                    smooth_iters,
                    periodic=(False, False),
                    normals_override=normals_override,
                    ghost_rows=ghost,
                )
                # Note: a previous "pole-row STL-normal override" used to run
                # here whenever do_blend was True, on the assumption that the
                # Laplacian smoothing's boundary clamp pulled the cap-equator
                # seam rows inward.  With seam-aware ghost rows in place that
                # clamp produces only sub-FP-noise drift, and the override
                # itself amplifies the C1 cap-equator seam deviation by a
                # factor of ~2 at k=NK-1.  See
                # tests/test_cap_gnomonic.py::TestBlendNormalsK_C1SeamSurvives
                # for the regression that pinned this down.  The override has
                # been removed; the alpha-weighted blend at lines above is
                # the only remaining STL-normal contribution.
                new_subs.append(new_q)

            # Re-stitch shared edges to absorb FP noise from the
            # Laplacian-smoothing boundary clamp.  With seam-aware normals
            # this should be sub-FP-noise; we keep the averaging as a
            # belt-and-braces guard against any residual asymmetry.
            for q in range(4):
                q_next = (q + 1) % 4
                shared = 0.5 * (new_subs[q][:, -1, :] + new_subs[q_next][:, 0, :])
                new_subs[q][:, -1, :] = shared
                new_subs[q_next][:, 0, :] = shared

            # Jacobian check for sub-blocks.
            # Junction corners produce O(-1e-6) noise from degenerate winding
            # at the 4-block junction; exclude them from the check.
            jmin_global = np.inf
            for q, blk in enumerate(sub_blocks):
                jmin = _check_jacobian(
                    blk.nodes[:, :, k, :],
                    new_subs[q],
                    periodic_j=False,
                    label=blk.label,
                    k=k,
                    exclude_corners=True,
                    tol=-1e-6,
                )
                jmin_global = min(jmin_global, jmin)

            # --- 2. Extract pole rows for cap boundary pins ---
            north_eq_rows = [new_subs[q][ni - 1, :, :] for q in range(4)]
            south_eq_rows = [new_subs[q][0, :, :] for q in range(4)]

            # --- 3. Advance cap blocks ---
            cap_north_ghost = _cap_ghost_rows("north", sub_k, ni)
            cap_south_ghost = _cap_ghost_rows("south", sub_k, ni)

            if do_blend and k < blend_normals_k:
                alpha = 1.0 - k / blend_normals_k
                cn_layer = cn_k
                nc_c = cn_layer.shape[0]
                grid_n = _compute_normals(
                    cn_layer,
                    periodic=(False, False),
                    ghost_rows=cap_north_ghost,
                )
                stl_n = _stl_normals(trimesh_obj, cn_layer.reshape(-1, 3)).reshape(nc_c, nc_c, 3)
                blended = alpha * stl_n + (1.0 - alpha) * grid_n
                mag = np.linalg.norm(blended, axis=-1, keepdims=True)
                cap_north_no: np.ndarray | None = blended / np.where(mag < 1e-15, 1.0, mag)

                cs_layer = cs_k
                grid_n = _compute_normals(
                    cs_layer,
                    periodic=(False, False),
                    ghost_rows=cap_south_ghost,
                )
                stl_n = _stl_normals(trimesh_obj, cs_layer.reshape(-1, 3)).reshape(nc_c, nc_c, 3)
                blended = alpha * stl_n + (1.0 - alpha) * grid_n
                mag = np.linalg.norm(blended, axis=-1, keepdims=True)
                cap_south_no: np.ndarray | None = blended / np.where(mag < 1e-15, 1.0, mag)
            else:
                cap_north_no = None
                cap_south_no = None

            new_cap_north = _march_cap_layer(
                cn_k,
                step,
                smooth_strength,
                smooth_iters,
                north_eq_rows,
                pole="north",
                normals_override=cap_north_no,
                ghost_rows=cap_north_ghost,
            )
            new_cap_south = _march_cap_layer(
                cs_k,
                step,
                smooth_strength,
                smooth_iters,
                south_eq_rows,
                pole="south",
                normals_override=cap_south_no,
                ghost_rows=cap_south_ghost,
            )

            jmin_cap_n = _check_jacobian(
                cap_north.nodes[:, :, k, :],
                new_cap_north,
                periodic_j=False,
                label="cap_north",
                k=k,
                exclude_corners=True,
                tol=-1e-6,
            )
            jmin_cap_s = _check_jacobian(
                cap_south.nodes[:, :, k, :],
                new_cap_south,
                periodic_j=False,
                label="cap_south",
                k=k,
                exclude_corners=True,
                tol=-1e-6,
            )
            jmin_global = min(jmin_global, jmin_cap_n, jmin_cap_s)

            # --- 4. Commit ---
            for q, blk in enumerate(sub_blocks):
                blk.nodes[:, :, k + 1, :] = new_subs[q]
            cap_north.nodes[:, :, k + 1, :] = new_cap_north
            cap_south.nodes[:, :, k + 1, :] = new_cap_south

            if layer_callback is not None:
                layer_callback(k + 1, jmin_global)

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    def write_geo(self, stem: str | Path) -> None:
        """Write each block as a separate ``<stem>_<label>.geo`` file."""
        from shore.io.geo import write_geo

        stem = Path(stem)
        for blk in self.blocks:
            write_geo(stem.parent / f"{stem.name}_{blk.label}.geo", blk.nodes)

    def write_vtk(self, stem: str | Path, with_jacobian: bool = True) -> None:
        """Write all blocks as a single PyVista MultiBlock ``.vtm`` file.

        Parameters
        ----------
        stem :
            Output path (without extension).  Writes ``<stem>.vtm``.
        with_jacobian :
            Attach per-cell Jacobian as CellData (default True).
        """
        from shore.io.vtk import write_caps_vtm

        path = Path(stem).with_suffix(".vtm")
        write_caps_vtm(path, [blk.nodes for blk in self.blocks], with_jacobian=with_jacobian)

    # ------------------------------------------------------------------
    # Inspection (in-memory)
    # ------------------------------------------------------------------

    def info(self) -> dict:
        """Return a dict of mesh metadata.

        Returns
        -------
        dict with keys:
            ``n_blocks``, ``labels``, and per-block dicts under each label
            with ``shape``, ``n_points``, ``n_cells``,
            ``xyz_min``, ``xyz_max``.
        """
        result: dict = {"n_blocks": len(self.blocks), "labels": [], "blocks": {}}
        for blk in self.blocks:
            ni, nj, nk, _ = blk.nodes.shape
            pts = blk.nodes.reshape(-1, 3)
            result["labels"].append(blk.label)
            result["blocks"][blk.label] = {
                "shape": (ni, nj, nk),
                "n_points": ni * nj * nk,
                "n_cells": (ni - 1) * (nj - 1) * (nk - 1),
                "xyz_min": pts.min(axis=0).tolist(),
                "xyz_max": pts.max(axis=0).tolist(),
            }
        return result

    def check(self, warn_threshold: float = 1e-6) -> dict:
        """Check cell Jacobians for all blocks layer by layer.

        Parameters
        ----------
        warn_threshold :
            Warn if jmin drops below this value (but stay positive).

        Returns
        -------
        dict with keys ``ok`` (bool) and ``blocks`` (per-block layer data).

        Raises
        ------
        MarchingError
            If any block has a negative Jacobian layer.
        """
        from shore.volume.jacobian import _jacobian_field_open

        result: dict = {"ok": True, "blocks": {}}
        for blk in self.blocks:
            jmins = []
            _ni, _nj, nk, _ = blk.nodes.shape
            # layer-by-layer min Jacobian
            for k in range(nk - 1):
                two = np.stack([blk.nodes[:, :, k, :], blk.nodes[:, :, k + 1, :]], axis=2)
                jac = _jacobian_field_open(two)
                jmins.append(float(jac.min()))
            any_neg = any(j <= 0.0 for j in jmins)
            any_warn = any(0.0 < j < warn_threshold for j in jmins)
            if any_neg:
                result["ok"] = False
            result["blocks"][blk.label] = {
                "jmins": jmins,
                "any_negative": any_neg,
                "any_warn": any_warn,
            }
        return result

    def export(
        self,
        stem: str | Path,
        with_jacobian: bool = True,
    ) -> None:
        """Export the mesh to VTK MultiBlock format (.vtm) for ParaView.

        Equivalent to :meth:`write_vtk`.

        Parameters
        ----------
        stem :
            Output path stem (extension set to ``.vtm`` automatically).
        with_jacobian :
            Attach per-cell Jacobian as CellData.
        """
        self.write_vtk(stem, with_jacobian=with_jacobian)

    def view(
        self,
        *,
        jupyter: bool = False,
        off_screen: bool = False,
        show_edges: bool = True,
        opacity: float = 0.6,
        colorby: str = "Jacobian",
    ) -> object:
        """Open an interactive PyVista window showing all blocks.

        Parameters
        ----------
        jupyter :
            Use HTML backend for inline Jupyter display.
        off_screen :
            Render off-screen (for testing).
        show_edges :
            Draw cell edges.
        opacity :
            Block opacity [0-1].
        colorby :
            Cell data field to colour by (``"Jacobian"`` or ``"none"``).

        Returns
        -------
        pyvista.Plotter
        """
        try:
            import pyvista as pv

            from shore.io.vtk import open_to_pyvista
        except ImportError as exc:
            raise ImportError(
                "pyvista is required for visualisation. Install with: pip install shore[vis]"
            ) from exc

        if jupyter:
            pv.set_jupyter_backend("html")

        pl = pv.Plotter(off_screen=off_screen)
        for blk in self.blocks:
            sg = open_to_pyvista(blk.nodes, with_jacobian=(colorby == "Jacobian"))
            pl.add_mesh(
                sg,
                scalars=colorby if colorby != "none" else None,
                cmap="RdYlGn",
                show_edges=show_edges,
                opacity=opacity,
            )

        pl.add_axes()
        pl.add_bounding_box()

        if not jupyter and not off_screen:
            pl.show()

        return pl
