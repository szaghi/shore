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
    _squareness_rows,
)
from shore.volume.hyperbolic import _compute_normals
from shore.volume.spacing import build_k_steps


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
        self.body     = body
        self.blocks   = blocks
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
        body    = Body.from_stl(path)
        surface = body.project(
            ni, nj,
            theta_cap=theta_cap,
            theta_cap_deg=theta_cap_deg,
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
            raise NotImplementedError(
                "Mesh.march() currently only supports CapsTopology."
            )
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
        cap_north  = self.blocks[4]
        cap_south  = self.blocks[5]

        nk = sub_blocks[0].nodes.shape[2]
        ni = sub_blocks[0].nodes.shape[0]

        _steps = steps if steps is not None else build_k_steps(nk, ds, growth)

        cap_squareness = self.topology.cap_squareness

        do_blend = blend_normals_k > 0
        if do_blend:
            from shore.volume.hyperbolic import _stl_normals
            trimesh_obj = self.body.mesh

        for k in range(nk - 1):
            step = float(_steps[k])

            # --- 1. Advance equatorial sub-blocks ---
            new_subs: list[np.ndarray] = []
            for blk in sub_blocks:
                layer_q = blk.nodes[:, :, k, :]
                ni_q, nj_q = layer_q.shape[:2]
                if do_blend and k < blend_normals_k:
                    alpha = 1.0 - k / blend_normals_k
                    grid_n = _compute_normals(layer_q, periodic=(False, False))
                    stl_n  = _stl_normals(
                        trimesh_obj, layer_q.reshape(-1, 3)
                    ).reshape(ni_q, nj_q, 3)
                    blended = alpha * stl_n + (1.0 - alpha) * grid_n
                    mag = np.linalg.norm(blended, axis=-1, keepdims=True)
                    normals_override: np.ndarray | None = blended / np.where(mag < 1e-15, 1.0, mag)
                else:
                    normals_override = None
                new_q = _march_one_layer(
                    layer_q, step, smooth_strength, smooth_iters,
                    periodic=(False, False),
                    normals_override=normals_override,
                )
                if do_blend:
                    # Laplacian smoothing corrupts the pole rows regardless of
                    # the blend window: the one-sided boundary clamp pulls i=0
                    # and i=-1 inward.  Recompute those rows post-smoothing
                    # using pure STL normals so they are always perpendicular
                    # to the body surface.
                    pole_pts = np.concatenate([layer_q[0], layer_q[-1]], axis=0)
                    stl_poles = _stl_normals(trimesh_obj, pole_pts)
                    new_q[0]  = layer_q[0]  + step * stl_poles[:nj_q]
                    new_q[-1] = layer_q[-1] + step * stl_poles[nj_q:]
                new_subs.append(new_q)

            # Re-stitch shared edges (SHARED views, but averaging removes
            # any floating-point divergence from the smoothing passes)
            for q in range(4):
                q_next = (q + 1) % 4
                shared = 0.5 * (new_subs[q][:, -1, :] + new_subs[q_next][:, 0, :])
                new_subs[q][:, -1, :]     = shared
                new_subs[q_next][:, 0, :] = shared

            # Jacobian check for sub-blocks.
            # Junction corners produce O(-1e-6) noise from degenerate winding
            # at the 4-block junction; exclude them from the check.
            jmin_global = np.inf
            for q, blk in enumerate(sub_blocks):
                jmin = _check_jacobian(
                    blk.nodes[:, :, k, :], new_subs[q],
                    periodic_j=False, label=blk.label, k=k,
                    exclude_corners=True, tol=-1e-6,
                )
                jmin_global = min(jmin_global, jmin)

            # --- 2. Extract pole rows for cap boundary pins ---
            north_eq_rows = [new_subs[q][ni - 1, :, :] for q in range(4)]
            south_eq_rows = [new_subs[q][0, :, :]      for q in range(4)]

            if cap_squareness != 0.0:
                # Squareness-modified rows used only as cap boundary pins.
                # Equatorial sub-blocks keep their unmodified positions to avoid
                # a geometric kink at the seam.
                north_eq_rows = _squareness_rows(north_eq_rows, cap_squareness)
                south_eq_rows = _squareness_rows(south_eq_rows, cap_squareness)
                for q in range(4):
                    q_next = (q + 1) % 4
                    shared_n = 0.5 * (north_eq_rows[q][-1] + north_eq_rows[q_next][0])
                    north_eq_rows[q][-1]       = shared_n
                    north_eq_rows[q_next][0]   = shared_n
                    shared_s = 0.5 * (south_eq_rows[q][-1] + south_eq_rows[q_next][0])
                    south_eq_rows[q][-1]       = shared_s
                    south_eq_rows[q_next][0]   = shared_s

            # --- 3. Advance cap blocks ---
            if do_blend and k < blend_normals_k:
                alpha = 1.0 - k / blend_normals_k
                cn_layer = cap_north.nodes[:, :, k, :]
                nc_c = cn_layer.shape[0]
                grid_n = _compute_normals(cn_layer, periodic=(False, False))
                stl_n  = _stl_normals(trimesh_obj, cn_layer.reshape(-1, 3)).reshape(nc_c, nc_c, 3)
                blended = alpha * stl_n + (1.0 - alpha) * grid_n
                mag = np.linalg.norm(blended, axis=-1, keepdims=True)
                cap_north_no: np.ndarray | None = blended / np.where(mag < 1e-15, 1.0, mag)

                cs_layer = cap_south.nodes[:, :, k, :]
                grid_n = _compute_normals(cs_layer, periodic=(False, False))
                stl_n  = _stl_normals(trimesh_obj, cs_layer.reshape(-1, 3)).reshape(nc_c, nc_c, 3)
                blended = alpha * stl_n + (1.0 - alpha) * grid_n
                mag = np.linalg.norm(blended, axis=-1, keepdims=True)
                cap_south_no: np.ndarray | None = blended / np.where(mag < 1e-15, 1.0, mag)
            else:
                cap_north_no = None
                cap_south_no = None

            new_cap_north = _march_cap_layer(
                cap_north.nodes[:, :, k, :], step,
                smooth_strength, smooth_iters,
                north_eq_rows, pole="north",
                normals_override=cap_north_no,
            )
            new_cap_south = _march_cap_layer(
                cap_south.nodes[:, :, k, :], step,
                smooth_strength, smooth_iters,
                south_eq_rows, pole="south",
                normals_override=cap_south_no,
            )

            jmin_cap_n = _check_jacobian(
                cap_north.nodes[:, :, k, :], new_cap_north,
                periodic_j=False, label="cap_north", k=k,
                exclude_corners=True, tol=-1e-6,
            )
            jmin_cap_s = _check_jacobian(
                cap_south.nodes[:, :, k, :], new_cap_south,
                periodic_j=False, label="cap_south", k=k,
                exclude_corners=True, tol=-1e-6,
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
        write_caps_vtm(path, [blk.nodes for blk in self.blocks],
                       with_jacobian=with_jacobian)

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
                two = np.stack(
                    [blk.nodes[:, :, k, :], blk.nodes[:, :, k + 1, :]], axis=2
                )
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
                "pyvista is required for visualisation. "
                "Install with: pip install shore[vis]"
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
