"""Mesh: multi-block orchestrator — projection, marching, I/O."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np

from shore.exceptions import ParameterError
from shore.geometry import Geometry
from shore.mesh.association import (
    BODY_GEOMETRY,
    Association,
    AssociationResolver,
    associate_wall,
    snap_block,
)
from shore.mesh.block import HexBlock
from shore.mesh.body import Body
from shore.mesh.face import FaceBC
from shore.mesh.spacing import Spacing1D
from shore.mesh.topology import Topology
from shore.topology.sphere_cap import (
    _check_jacobian,
    _march_cap_layer,
    _march_one_layer,
)
from shore.volume.spacing import build_k_steps

# ---------------------------------------------------------------------------
# Per-edge k-spacing resolver
# ---------------------------------------------------------------------------
# All blocks in a single Mesh share one k-step schedule because adjacent
# blocks have k-coupled SHARED / DIRICHLET seams that must remain coincident
# layer-by-layer.  The resolver below picks the schedule from per-edge data
# when available, falling back to the (ds, growth) defaults.
#
# Precedence (highest first):
#   1. Caller-supplied ``steps`` argument (overrides everything).
#   2. Block 0's k-edges, if they all carry a consistent Spacing1D.  Other
#      blocks' k-edges are checked for agreement; disagreement → warn and
#      fall back.
#   3. Block 0's ``spacing_k`` attribute (the topology-level default).
#   4. The ``(ds, growth)`` geometric defaults.

_K_EDGE_NAMES = ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi")

BLEND_REMOVED = (
    "blend_normals_k was removed (ADR 0034; deprecated by ADR 0019): blending body normals "
    "into the marching direction showed no benefit over the grid normals and folded "
    "C-grids. Delete the option; a project that used it loads only without replay."
)
BLEND_ZERO = "blend_normals_k was removed (ADR 0034) and 0 is all it accepts; delete the option."


def check_blend_normals_k(blend_normals_k: int | None) -> None:
    """Refuse the removed ``blend_normals_k`` (ADR 0034); warn on an explicit 0."""
    if blend_normals_k is None:
        return
    if blend_normals_k != 0:
        raise ParameterError(BLEND_REMOVED)
    warnings.warn(BLEND_ZERO, DeprecationWarning, stacklevel=3)


def _block_user_k_spacing(block: HexBlock) -> Spacing1D | None:
    """Return *block*'s k-spacing if user-overridden on every k-edge, else None.

    Override detection is by identity: ``HexBlock.__init__`` assigns the same
    ``Spacing1D`` instance (``block._default_spacing_k``) to all 4 k-edges, so
    any edge whose ``spacing`` is no longer that instance has been replaced
    by the caller post-construction.

    Policy:
        - All 4 k-edges still pointing at the default → no override (return None).
        - All 4 k-edges replaced with the same instance → return that instance.
        - Mixed (some replaced, some still default; or replaced with different
          instances) → emit a warning and return None (fall back to default).
    """
    default = block._default_spacing_k
    edges = [block.edge(name).spacing for name in _K_EDGE_NAMES]

    overridden = [e is not default for e in edges]
    if not any(overridden):
        return None  # all edges still on the default — no user override

    if not all(overridden):
        warnings.warn(
            f"Block '{block.label}': k-spacing overridden on some k-edges but "
            "not others; falling back to block default.  Set all 4 k-edges "
            "(i_lo_j_lo, i_hi_j_lo, i_lo_j_hi, i_hi_j_hi) to the same "
            "Spacing1D instance for a per-block override.",
            UserWarning,
            stacklevel=3,
        )
        return None

    first = edges[0]
    if any(e is not first and not _spacings_equal(e, first) for e in edges[1:]):
        warnings.warn(
            f"Block '{block.label}': the 4 k-edges carry different Spacing1D "
            "values; falling back to block default.  All 4 k-edges must agree "
            "for a per-block override.",
            UserWarning,
            stacklevel=3,
        )
        return None

    return first


def _resolve_k_steps(
    blocks: list[HexBlock],
    nk: int,
    ds: float,
    growth: float,
    steps: np.ndarray | None,
) -> np.ndarray:
    """Resolve the wall-normal k-step schedule for a list of blocks.

    Precedence (highest first):
        1. Caller-supplied ``steps`` argument.
        2. Block 0's k-edges — if all 4 carry a user-overridden ``Spacing1D``
           that agrees, use it.
        3. Block 0's ``spacing_k`` default if it has been customised away from
           the constructor fallback.
        4. The ``(ds, growth)`` geometric defaults.

    All blocks in the mesh share one schedule because k-coupled seams require
    coincident node positions layer-by-layer.  When other blocks disagree
    with block 0's choice, a ``UserWarning`` is emitted but block 0 still wins.
    """
    if steps is not None:
        return steps

    block0 = blocks[0]
    canonical = _block_user_k_spacing(block0)

    # Cross-block consistency: every block's k-edges should agree with block 0.
    for blk in blocks[1:]:
        their = _block_user_k_spacing(blk)
        if their is not None and not _spacings_equal(their, canonical):
            warnings.warn(
                f"Block '{blk.label}' k-edges request a different spacing than "
                f"block '{block0.label}' ({their} vs {canonical}).  K-coupled "
                f"seams require one schedule; using block '{block0.label}'.",
                UserWarning,
                stacklevel=3,
            )

    if canonical is None:
        # No per-edge override: fall through to the legacy (ds, growth) path.
        return build_k_steps(nk, ds, growth)
    return canonical.build_steps(nk)


def _spacings_equal(a: Spacing1D | None, b: Spacing1D | None) -> bool:
    """Compare two Spacing1D instances by *value* (dataclass equality)."""
    if a is None or b is None:
        return a is b
    return a == b


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
        body: Body | None,
        blocks: list[HexBlock],
        topology: Topology,
    ) -> None:
        self.body = body
        self.blocks = blocks
        self.topology = topology
        if body is not None:
            # every topology's body-fitted wall is its WALL k_lo face (C-grid
            # wake cuts and channel faces are SHARED): associate it with the body
            body_assoc = Association.body()
            for block in blocks:
                if block.face("k_lo").bc is FaceBC.WALL:
                    associate_wall(block, body_assoc)

    @property
    def geometries(self) -> dict[str, Geometry]:
        """Geometry map for associations (``"body"`` when the mesh has one)."""
        return {} if self.body is None else {BODY_GEOMETRY: self.body.geometry}

    def snap(self, exact: bool = False) -> None:
        """Project every associated node onto its geometry (ADR 0003).

        Faces, then edges, then vertices (vertex -> edge -> face priority).
        ``exact=True`` uses the exact CAD geometry (tier 2) instead of its
        tessellation; for STL bodies the two are the same.
        """
        resolver = AssociationResolver(self.geometries, exact=exact)
        for block in self.blocks:
            snap_block(block, resolver)

    @classmethod
    def from_array(
        cls,
        surface: np.ndarray,
        topology: Topology,
        nk: int,
        spacing_k: Spacing1D | None = None,
        **kwargs,
    ) -> Mesh:
        """Build a Mesh from a pre-built surface array (no STL).

        Used by topologies whose input is a numeric array (e.g. a
        cross-section polyline) rather than a triangulated body.  The
        resulting :class:`Mesh` has ``body = None``; methods that
        depend on ``self.body`` (e.g. :meth:`snap`) will fail loudly if
        called.

        Parameters
        ----------
        surface :
            Topology-specific surface array.  For
            :class:`CGridTipCap`: a ``(ni, 3)`` cross-section polyline.
        topology :
            The topology builder.
        nk :
            Topology-specific node count.  For :class:`CGridTipCap`,
            the cap's spanwise (k-axis) node count.
        spacing_k :
            Optional wall-normal / spanwise spacing law forwarded to
            ``topology.build``.
        **kwargs :
            Forwarded verbatim to ``topology.build``.

        Returns
        -------
        Mesh
            With ``body=None`` and the topology-built blocks.
        """
        blocks = topology.build(
            surface=surface,
            nk=nk,
            spacing_k=spacing_k,
            **kwargs,
        )
        return cls(body=None, blocks=blocks, topology=topology)

    @classmethod
    def from_stl(
        cls, path: str | Path, topology: Topology, *args: object, **kwargs: object
    ) -> Mesh:
        """Alias of ``from_file`` (which also reads STEP/IGES/BREP)."""
        return cls.from_file(path, topology, *args, **kwargs)  # type: ignore[arg-type]

    @classmethod
    def from_file(
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
        """Load a body, project to surface, wire blocks.  Does not march.

        Parameters
        ----------
        path :
            STL, or a CAD file (STEP/IGES/BREP, needs ``shore-mesh[cad]``).
            For CAD the k=0 ray casts are snapped to the exact surface.
        topology :
            Wired topology builder (e.g. ``CubedSphereTopology()``).
        ni, nj :
            Surface node counts (i and j directions, the array shape).
        nk :
            Wall-normal node count, including the wall row at k=0.  With
            ``nk = 30`` the marcher extrudes 29 cells from the wall.
            (SHORE convention: ``ni``, ``nj``, ``nk`` are always node
            counts; the only file that stores cell counts is the binary
            ``.grd`` header, where the conversion is automatic.)
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
        body = Body.from_file(path)

        if theta_cap_deg is not None and theta_cap is None:
            _theta = float(np.deg2rad(theta_cap_deg))
        elif theta_cap is not None:
            _theta = float(theta_cap)
        else:
            _theta = float(0.05 * np.pi)

        from shore.mesh.topology import (
            CGridTopology,
            CHTopology,
            CubedSphereTopology,
            OGridTopology,
            OHTopology,
        )

        if isinstance(topology, (CGridTopology, CHTopology)):
            kind = type(topology).__name__
            raise ValueError(
                f"Mesh.from_file does not support {kind} — its "
                "parameters (ni_body, ni_wake, n_stations, wake_length, "
                "span_axis, freestream_world) don't fit the (ni, nj, nk) "
                "signature.  Build the block list directly:\n\n"
                "    import trimesh\n"
                f"    from shore.mesh import Body, Mesh, {kind}\n"
                "    body = Body.from_file(path)\n"
                f"    topo = {kind}(...)\n"
                "    blocks = topo.build(\n"
                "        surface=None, nk=nk, mesh=body.mesh,\n"
                "        ni_body=..., ni_wake=..., n_stations=...,\n"
                "        wake_length=...,\n"
                "    )\n"
                "    mesh = Mesh(body=body, blocks=blocks, topology=topo)\n"
            )

        if isinstance(topology, (CubedSphereTopology, OGridTopology, OHTopology)):
            # Both topologies project through the cached BVH inside
            # ``topology.build`` — single ray-cast pass.
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
        blend_normals_k: int | None = None,
        smooth_mode: str = "limited",
        march_scheme: str = "algebraic",
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
            **Removed** (ADR 0034): only ``None`` (the default) or 0 are
            accepted; 0 warns. It blended body normals into the first layers'
            marching directions, with no benefit (ADR 0019).
        smooth_mode :
            ``'limited'`` (default): the Laplacian smoothing motion is
            limited per node where it would damage the layer (it may pull a
            node back by at most a quarter of the step, and move it sideways
            by at most a quarter of the step); elsewhere it is unchanged, so grids that marched
            cleanly before are identical. ``'legacy'``: never limit (kept to
            reproduce old grids); on curved walls the layer shrinks inward by
            about ``smooth_strength * h**2 * curvature`` and folds layers
            thinner than that.
        march_scheme :
            ``'algebraic'`` (default): extrude along the local normal, then
            Laplacian smoothing (``smooth_*``). ``'hyperbolic'``: the
            Steger-Chaussee scheme (``volume.hyperbolic_scheme``): k-lines
            orthogonal to each layer, prescribed cell volumes, implicit
            coupling and dissipation that grow with step / cell size; the
            ``smooth_*`` arguments are not used.

        Raises
        ------
        MarchingError
            If a negative-Jacobian cell is detected.
        NotImplementedError
            If the topology has no registered marcher.
        ParameterError
            ``blend_normals_k`` other than 0 (removed, ADR 0034).

        Warns
        -----
        DeprecationWarning
            ``blend_normals_k=0`` given explicitly.
        """
        if march_scheme not in ("algebraic", "hyperbolic"):
            raise ParameterError(
                f"march_scheme must be 'algebraic' or 'hyperbolic', got {march_scheme!r}."
            )
        check_blend_normals_k(blend_normals_k)
        from shore.mesh.topology import (
            CGridTipCap,
            CGridTopology,
            CHTopology,
            CubedSphereTopology,
            OGridTopology,
            OHTopology,
        )

        if isinstance(self.topology, CGridTipCap):
            # The tip cap is fully volumetric at build time — nothing to march.
            return

        if isinstance(self.topology, OHTopology):
            # March the wrap (first 6 blocks) using the cubed-sphere marcher;
            # then build channel-block interiors against the bounding box.
            wrap_blocks = self.blocks[:6]
            self._march_caps_blocks(
                wrap_blocks,
                ds=ds,
                growth=growth,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                steps=steps,
                layer_callback=layer_callback,
            )
            self.topology.fill_channel_blocks(self.blocks)
        elif isinstance(self.topology, CHTopology):
            # March the C-grid wrap (first 3 blocks; 4th is wake_base which
            # is fully determined at build time and not marched).  Then fill
            # the channel-block interiors against the chord/normal frame.
            self._march_cgrid_three(
                ds=ds,
                growth=growth,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                steps=steps,
                layer_callback=layer_callback,
            )
            self.topology.fill_channel_blocks(self.blocks)
        elif isinstance(self.topology, CubedSphereTopology):
            self._march_caps(
                ds=ds,
                growth=growth,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                steps=steps,
                layer_callback=layer_callback,
            )
        elif isinstance(self.topology, OGridTopology):
            self._march_ogrid(
                ds=ds,
                growth=growth,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                steps=steps,
                layer_callback=layer_callback,
            )
        elif isinstance(self.topology, CGridTopology):
            self._march_cgrid(
                ds=ds,
                growth=growth,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                steps=steps,
                layer_callback=layer_callback,
            )
        else:
            raise NotImplementedError(
                f"Mesh.march(): no marcher registered for {type(self.topology).__name__}."
            )

    def _march_caps(
        self,
        ds: float,
        growth: float,
        smooth_strength: float,
        smooth_iters: int,
        steps: np.ndarray | None,
        layer_callback: Callable[[int, float], None] | None,
        smooth_mode: str = "limited",
        march_scheme: str = "algebraic",
    ) -> None:
        """Simultaneous marching for the 6-block CubedSphereTopology."""
        self._march_caps_blocks(
            self.blocks,
            ds=ds,
            growth=growth,
            smooth_strength=smooth_strength,
            smooth_mode=smooth_mode,
            march_scheme=march_scheme,
            smooth_iters=smooth_iters,
            steps=steps,
            layer_callback=layer_callback,
        )

    def _march_caps_blocks(
        self,
        wrap_blocks: list[HexBlock],
        ds: float,
        growth: float,
        smooth_strength: float,
        smooth_iters: int,
        steps: np.ndarray | None,
        layer_callback: Callable[[int, float], None] | None,
        smooth_mode: str = "limited",
        march_scheme: str = "algebraic",
    ) -> None:
        """Simultaneous marching for a 6-block cubed-sphere wrap.

        Operates on the supplied block list rather than ``self.blocks``,
        so it can drive the wrap portion of a hybrid topology like
        :class:`shore.mesh.OHTopology` (which has 12 blocks total — only
        the first 6 get marched).
        """
        sub_blocks = wrap_blocks[:4]
        cap_north = wrap_blocks[4]
        cap_south = wrap_blocks[5]

        nk = sub_blocks[0].nodes.shape[2]
        ni = sub_blocks[0].nodes.shape[0]

        _steps = _resolve_k_steps(wrap_blocks, nk, ds, growth, steps)

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

                # Ghost rows for sub_q at this layer (pre-march).  i_lo and
                # i_hi come from the cap blocks (rotated/reversed per the
                # _overwrite_cap_boundaries mapping); j_lo and j_hi come from
                # adjacent sub-blocks via the periodic ring.
                ghost = _sub_ghost_rows(q, sub_k, cn_k, cs_k)
                new_q = _march_one_layer(
                    layer_q,
                    step,
                    smooth_strength,
                    smooth_iters,
                    periodic=(False, False),
                    ghost_rows=ghost,
                    smooth_mode=smooth_mode,
                    march_scheme=march_scheme,
                )
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

            new_cap_north = _march_cap_layer(
                cn_k,
                step,
                smooth_strength,
                smooth_iters,
                north_eq_rows,
                pole="north",
                ghost_rows=cap_north_ghost,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
            )
            new_cap_south = _march_cap_layer(
                cs_k,
                step,
                smooth_strength,
                smooth_iters,
                south_eq_rows,
                pole="south",
                ghost_rows=cap_south_ghost,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
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

    def _march_cgrid_three(
        self,
        ds: float,
        growth: float,
        smooth_strength: float,
        smooth_iters: int,
        steps: np.ndarray | None,
        layer_callback: Callable[[int, float], None] | None,
        smooth_mode: str = "limited",
        march_scheme: str = "algebraic",
    ) -> None:
        """Hyperbolic march for the 3-block C-grid (upper / body / lower).

        Per layer: snapshot all three blocks; supply ghost rows so
        each block's i-tangent at the upper-TE / lower-TE seams is a
        centred difference against the partner block; advance each
        block; restitch the SHARED TE columns to absorb FP noise from
        the Laplacian boundary clamp.

        The wake-cut SHARED seam (``upper.k_lo ↔ lower.k_lo``) is
        purely structural — it tells Xall the two faces share data,
        but the marcher itself doesn't need cross-block coordination
        at k=0 (it's the seed surface, never modified) or at k>0
        (the SHARED face only exists at k=0, no seam to cross above).
        Per-block Jacobian guards skip k=0 on the wake blocks (upper
        and lower) because the wake-cut k=0 cells have collapsing
        volumes for sharp TE — see the C-grid topology docs.

        For the 4-block blunt-TE layout, the wake-base block W
        (``self.blocks[3]``) is **not marched** — its geometry is
        fully determined at build time by linear interpolation
        across the TE thickness gap.  Only blocks U / B / L are
        advanced layer by layer.
        """
        upper, body, lower = self.blocks[:3]
        nk = body.nodes.shape[2]
        _steps = _resolve_k_steps(self.blocks, nk, ds, growth, steps)

        periodic = (False, False)

        for k in range(nk - 1):
            step = float(_steps[k])

            # Snapshots of every block's k-layer — the consistent
            # neighbour data for centred-difference normals at the
            # TE seams.
            u_k = upper.nodes[:, :, k, :]
            b_k = body.nodes[:, :, k, :]
            l_k = lower.nodes[:, :, k, :]

            # ----- ghost rows at the TE seams ------------------------
            # Each block's i-tangent at its own seam needs the
            # one-row-inside-the-partner data.
            upper_ghost = {"i_hi": b_k[1, :, :]}  # body row 1 (after upper-TE)
            body_ghost = {
                "i_lo": u_k[-2, :, :],  # upper row -2 (before upper-TE)
                "i_hi": l_k[1, :, :],  # lower row 1 (after lower-TE)
            }
            lower_ghost = {"i_lo": b_k[-2, :, :]}  # body row -2 (before lower-TE)

            # ----- advance each block --------------------------------
            new_u = _march_one_layer(
                u_k,
                step,
                smooth_strength,
                smooth_iters,
                periodic=periodic,
                ghost_rows=upper_ghost,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
            )
            new_b = _march_one_layer(
                b_k,
                step,
                smooth_strength,
                smooth_iters,
                periodic=periodic,
                ghost_rows=body_ghost,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
            )
            new_l = _march_one_layer(
                l_k,
                step,
                smooth_strength,
                smooth_iters,
                periodic=periodic,
                ghost_rows=lower_ghost,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
            )

            # ----- re-stitch SHARED TE columns ----------------------
            # FP noise from the Laplacian boundary clamp can let the
            # upper-TE column drift between U and B (likewise B↔L);
            # average to keep the seams bit-coincident.
            te_upper = 0.5 * (new_u[-1, :, :] + new_b[0, :, :])
            new_u[-1, :, :] = te_upper
            new_b[0, :, :] = te_upper
            te_lower = 0.5 * (new_b[-1, :, :] + new_l[0, :, :])
            new_b[-1, :, :] = te_lower
            new_l[0, :, :] = te_lower

            # ----- per-block Jacobian guards -------------------------
            # Skip k=0 on upper / lower because the wake-cut at k=0
            # has zero volume (sharp TE) — that's geometric, not a bug.
            jmin_u = _check_jacobian(
                layer=u_k,
                new_layer=new_u,
                periodic_j=False,
                label="cgrid.upper",
                k=k,
                exclude_corners=False,
                tol=-1e-6,
            )
            jmin_b = _check_jacobian(
                layer=b_k,
                new_layer=new_b,
                periodic_j=False,
                label="cgrid.body",
                k=k,
                exclude_corners=False,
                tol=-1e-6,
            )
            jmin_l = _check_jacobian(
                layer=l_k,
                new_layer=new_l,
                periodic_j=False,
                label="cgrid.lower",
                k=k,
                exclude_corners=False,
                tol=-1e-6,
            )
            jmin = float(min(jmin_u, jmin_b, jmin_l))

            # Commit.
            upper.nodes[:, :, k + 1, :] = new_u
            body.nodes[:, :, k + 1, :] = new_b
            lower.nodes[:, :, k + 1, :] = new_l

            if layer_callback is not None:
                layer_callback(k + 1, jmin)

    def _march_cgrid(
        self,
        ds: float,
        growth: float,
        smooth_strength: float,
        smooth_iters: int,
        steps: np.ndarray | None,
        layer_callback: Callable[[int, float], None] | None,
        smooth_mode: str = "limited",
        march_scheme: str = "algebraic",
    ) -> None:
        """Hyperbolic march for ``CGridTopology`` (1, 3, or 4 blocks; or
        4 + N cap blocks under the conformal-cap path).

        With cap_tips != 'none', the first 4 blocks are the C-grid
        (upper, body, lower, wake_base) and any further blocks are
        cap blocks built fully volumetric at construction time —
        they do not march.  ``_march_cgrid_three`` already restricts
        itself to the first 3 blocks so only the C-grid gets marched.
        """
        # The 1-block legacy path has a single block labelled 'cgrid'.
        # Anything else (3, 4, or 4 + N cap blocks) goes through the
        # 3/4-block marcher which only touches the first 3 entries.
        if len(self.blocks) > 1:
            self._march_cgrid_three(
                ds=ds,
                growth=growth,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                steps=steps,
                layer_callback=layer_callback,
            )
            return
        # blocks == 1: legacy single-block path.
        block = self.blocks[0]
        nk = block.nodes.shape[2]
        _steps = _resolve_k_steps(self.blocks, nk, ds, growth, steps)

        periodic = (False, False)

        for k in range(nk - 1):
            step = float(_steps[k])
            layer = block.nodes[:, :, k, :]

            new_layer = _march_one_layer(
                layer,
                step=step,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                periodic=periodic,
            )

            jmin = _check_jacobian(
                layer=layer,
                new_layer=new_layer,
                periodic_j=False,
                label="cgrid",
                k=k,
                exclude_corners=False,
                tol=-1e-6,
            )

            block.nodes[:, :, k + 1, :] = new_layer

            if layer_callback is not None:
                layer_callback(k + 1, jmin)

    def _march_ogrid(
        self,
        ds: float,
        growth: float,
        smooth_strength: float,
        smooth_iters: int,
        steps: np.ndarray | None,
        layer_callback: Callable[[int, float], None] | None,
        smooth_mode: str = "limited",
        march_scheme: str = "algebraic",
    ) -> None:
        """Single-block hyperbolic march for ``OGridTopology``.

        The j-direction is periodic; the i-boundaries (south / north
        polar singular axis rows) are clamped during Laplacian
        smoothing.
        """
        block = self.blocks[0]
        nk = block.nodes.shape[2]
        _steps = _resolve_k_steps(self.blocks, nk, ds, growth, steps)

        periodic = (False, True)

        for k in range(nk - 1):
            step = float(_steps[k])
            layer = block.nodes[:, :, k, :]

            new_layer = _march_one_layer(
                layer,
                step=step,
                smooth_strength=smooth_strength,
                smooth_mode=smooth_mode,
                march_scheme=march_scheme,
                smooth_iters=smooth_iters,
                periodic=periodic,
            )

            jmin = _check_jacobian(
                layer=layer,
                new_layer=new_layer,
                periodic_j=True,
                label="ogrid",
                k=k,
                exclude_corners=False,
                tol=-1e-6,
            )

            block.nodes[:, :, k + 1, :] = new_layer

            if layer_callback is not None:
                layer_callback(k + 1, jmin)

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
        write_caps_vtm(
            path,
            [blk.nodes for blk in self.blocks],
            with_jacobian=with_jacobian,
            labels=[blk.label for blk in self.blocks],
        )

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
