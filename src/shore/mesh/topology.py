"""Topology ABC and concrete topologies (cubed-sphere, O-grid, C-grid, OH-/CH-hybrid)."""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import AxisMap, FaceBC, FaceName
from shore.mesh.spacing import Spacing1D
from shore.surface.cap import build_cap_k0
from shore.surface.cgrid import build_cgrid_k0
from shore.surface.equator import build_equator
from shore.topology.cgrid_cap import TIP_STATION_MARGIN
from shore.topology.ch_channel import outer_quads_for_cgrid, resolve_frame
from shore.topology.oh_channel import build_channel_block
from shore.topology.sphere_cap import _adjust_nj

#: The cubed sphere's default cap angle (ADR 0054): 30 deg marches the cap
#: corners without folds and joins cap and equator with C1 spacing on the
#: bodies and sizes measured; 0.05 pi folded them (ADR 0052).
CUBED_SPHERE_THETA_CAP = float(np.pi / 6)
#: The cap angle before ADR 0054, still the O-grid's polar opening.
LEGACY_THETA_CAP = float(0.05 * np.pi)


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


def _build_cubed_sphere_blocks(
    *,
    ni: int,
    nj: int,
    nk: int,
    total_k: int,
    cap_proj: object,
    spacing_k: Spacing1D | None,
    theta_cap: float,
    i_spacing: str,
    i_beta: float,
    j_spacing: str,
    j_beta: float,
    phi_offset: float,
) -> list[HexBlock]:
    """Allocate + wire the 6 cubed-sphere blocks.

    ``total_k`` controls the depth of the underlying numpy buffers.  When
    ``total_k > nk`` the extra slots at the end can be carved out as
    *channel* views by callers (e.g. :class:`OHTopology`) so that the
    wrap-channel seam is view-aliased rather than copied.  ``CubedSphereTopology``
    passes ``total_k == nk`` and gets the original behaviour.
    """
    nc = nj // 4 + 1
    _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

    # Cap k=0 seeds — gnomonic flat-square projection.
    cap_south_k0 = build_cap_k0(
        theta_cap=theta_cap,
        nc=nc,
        pole="south",
        projector=cap_proj,
        j_spacing=j_spacing,
        j_beta=j_beta,
        phi_offset=phi_offset,
    )
    cap_north_k0 = build_cap_k0(
        theta_cap=theta_cap,
        nc=nc,
        pole="north",
        projector=cap_proj,
        j_spacing=j_spacing,
        j_beta=j_beta,
        phi_offset=phi_offset,
    )

    # Equator k=0 — great-circle slerp + ray-cast.
    equator_k0 = build_equator(
        cap_south=cap_south_k0,
        cap_north=cap_north_k0,
        ni=ni,
        theta_cap=theta_cap,
        projector=cap_proj,
        i_spacing=i_spacing,
        i_beta=i_beta,
    )

    # Equatorial ring buffer of depth total_k.  Sub-block sees [..., :nk, :].
    ring_nj = 4 * (nc - 1) + 1
    ring = np.zeros((ni, ring_nj, total_k, 3), dtype=np.float64)
    sector = nc - 1
    sub_nodes_full: list[np.ndarray] = []
    for q in range(4):
        j0 = q * sector
        sub_nodes_full.append(ring[:, j0 : j0 + nc, :, :])

    # k=0 fill from equator surface (only writes layer 0 of the deeper buffer).
    for q in range(4):
        j0 = q * sector
        j_indices = [*range(j0, j0 + sector), (j0 + sector) % nj]
        sub_nodes_full[q][:, :, 0, :] = equator_k0[:, j_indices, :]
    ring[:, -1, 0, :] = ring[:, 0, 0, :]

    # Cap blocks of depth total_k.
    cap_north_full = np.zeros((nc, nc, total_k, 3), dtype=np.float64)
    cap_south_full = np.zeros((nc, nc, total_k, 3), dtype=np.float64)
    cap_north_full[:, :, 0, :] = cap_north_k0
    cap_south_full[:, :, 0, :] = cap_south_k0

    # The wrap-block views: depth nk, addressing only the first nk slots.
    sub_wrap_views = [arr[:, :, :nk, :] for arr in sub_nodes_full]
    cap_north_wrap = cap_north_full[:, :, :nk, :]
    cap_south_wrap = cap_south_full[:, :, :nk, :]

    sub_blocks: list[HexBlock] = []
    for q in range(4):
        sub_blocks.append(
            HexBlock(
                label=f"sub{q}",
                nodes=sub_wrap_views[q],
                periodic_i=False,
                periodic_j=False,
                spacing_k=_sp_k,
            )
        )
    cap_north_blk = HexBlock(
        label="cap_north",
        nodes=cap_north_wrap,
        periodic_i=False,
        periodic_j=False,
        spacing_k=_sp_k,
    )
    cap_south_blk = HexBlock(
        label="cap_south",
        nodes=cap_south_wrap,
        periodic_i=False,
        periodic_j=False,
        spacing_k=_sp_k,
    )

    # Sub-block adjacency: sub_q.j_hi ↔ sub_{q+1 mod 4}.j_lo (SHARED).
    for q in range(4):
        q_next = (q + 1) % 4
        sub_blocks[q]._set_face_bc(
            "j_hi",
            FaceBC.SHARED,
            partner=sub_blocks[q_next].face("j_lo"),
            axis_map=AxisMap.identity("j_hi"),
        )
        sub_blocks[q_next]._set_face_bc(
            "j_lo",
            FaceBC.SHARED,
            partner=sub_blocks[q].face("j_hi"),
            axis_map=AxisMap.identity("j_lo"),
        )

    # Cap-equator seams: 8 DIRICHLET pairings.
    sub_face: FaceName
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

    blocks = [*sub_blocks, cap_north_blk, cap_south_blk]

    # Stash the underlying buffers on the blocks so OHTopology can reach
    # them to carve channel views.  Only used when total_k > nk.
    if total_k > nk:
        for q in range(4):
            blocks[q]._oh_full_nodes = sub_nodes_full[q]  # type: ignore[attr-defined]
        cap_north_blk._oh_full_nodes = cap_north_full  # type: ignore[attr-defined]
        cap_south_blk._oh_full_nodes = cap_south_full  # type: ignore[attr-defined]

    return blocks


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
        theta_cap: float = CUBED_SPHERE_THETA_CAP,
        i_spacing: str = "uniform",
        i_beta: float = 3.0,
        j_spacing: str = "uniform",
        j_beta: float = 3.0,
        phi_offset: float = 0.0,
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

        return _build_cubed_sphere_blocks(
            ni=ni,
            nj=nj,
            nk=nk,
            total_k=nk,
            cap_proj=cap_proj,
            spacing_k=spacing_k,
            theta_cap=theta_cap,
            i_spacing=i_spacing,
            i_beta=i_beta,
            j_spacing=j_spacing,
            j_beta=j_beta,
            phi_offset=phi_offset,
        )


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
        nodes = np.zeros((ni, nj, nk, 3), dtype=np.float64)
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


class CGridTopology(Topology):
    """C-grid topology for sharp-TE bodies (airfoils, hydrofoils).

    Two block layouts are supported, selected via the ``blocks``
    argument to :meth:`build`:

    **3-block (default, Xall-correct).** Three hex blocks ``upper``,
    ``body``, ``lower``:

    - ``upper`` ``(ni_wake, n_stations, nk)`` — sweeps i from the
      upper-wake far-field tip toward the upper-TE.
    - ``body`` ``(ni_body, n_stations, nk)`` — sweeps i from the
      upper-TE around the LE to the lower-TE.
    - ``lower`` ``(ni_wake, n_stations, nk)`` — sweeps i from the
      lower-TE toward the lower-wake far-field tip.

    Three SHARED seams:

    - ``upper.i_hi  ↔ body.i_lo``  — upper-TE column (same axis,
      identity AxisMap on along-seam ``(j, k)``).
    - ``body.i_hi   ↔ lower.i_lo`` — lower-TE column (same).
    - ``upper.k_lo  ↔ lower.k_lo`` — wake cut.  Same-axis (k↔k) but
      ``i`` reverses on the partner side: upper's ``i = 0`` is the
      far-field tip while lower's ``i = 0`` is the TE.  AxisMap on
      ``upper.k_lo``: ``from_pairs(("i", True), ("j", False))``.

    For sharp TE the wake-cut nodes coincide bit-exactly in space and
    SHARED is a true memory-aliased view; for blunt TE the upper and
    lower wake k=0 strips are parallel-but-separated and the 3-block
    topology is rejected (use the 4-block variant instead — see the
    ``blocks=4`` arm of :meth:`build`).

    **1-block (legacy, retained for opt-in compatibility).** A single
    hex block of shape ``(ni_total, n_stations, nk)`` whose i-axis
    sweeps ``upper-wake-far-tip → upper-TE → suction → LE → pressure
    → lower-TE → lower-wake-far-tip`` as one continuous parametric
    line.  The wake cut is encoded as an i-direction overlap: the
    upper wake column at ``i = 0..ni_wake - 1`` and the lower wake
    column at ``i = ni_wake + ni_body - 1..ni_total - 1`` are
    geometrically coincident at k=0 but topologically distinct in i.
    Xall cannot apply BCs cleanly on the wake-cut faces in this
    layout (the j-counter "flips direction" along the overlapping
    edge); use ``blocks=3`` instead for new work.

    Conventions, both layouts:

    - ``span_axis`` defaults to ``"z"``.
    - Wake far-field tips (``upper.i_lo``, ``lower.i_hi`` for the
      3-block layout; ``i_lo`` and ``i_hi`` for the 1-block layout)
      are ``FREE`` by default — typically Chimera fringes (BC code
      ``40``) or extrapolation (``-9``) in a downstream ``cc.par``.
    - ``j_lo`` and ``j_hi`` on every block are wing-tip ``FREE``.
    - ``k_lo`` is ``WALL`` on the body block; on upper / lower it
      is ``FREE`` (the wake-cut faces become SHARED in 3-block
      layout, so the FREE default only ever surfaces in the
      legacy 1-block path).
    - ``k_hi`` is ``FREE`` (outer / Chimera fringe).

    See Also
    --------
    :func:`shore.surface.cgrid.build_cgrid_k0` :
        Builds the k=0 layer fed to this topology's :meth:`build`.
    """

    def build(
        self,
        surface: np.ndarray | None,
        nk: int,
        ni_body: int | None = None,
        ni_wake: int | None = None,
        n_stations: int | None = None,
        wake_length: float | None = None,
        span_axis: str = "z",
        freestream_world: np.ndarray | None = None,
        te_dihedral_deg: float = 30.0,
        wake_growth: float = 1.15,
        spacing_k: Spacing1D | None = None,
        mesh: object | None = None,
        blocks: int = 3,
        nk_w: int = 3,
        cap_tips: str = "none",
        cap_te_thickness: float | None = None,
        cap_te_thickness_frac: float = 0.002,
        cap_le_frac: float = 0.05,
        cap_te_frac: float | None = None,
        cap_tip_inset: float | None = None,
        cap_n_cap_j: int = 12,
        cap_j_spacing: str = "tanh",
        cap_j_beta: float = 4.0,
        cap_nk: int | None = None,
        cap_core_frac: float = 0.5,
        cap_flare_deg: float = 30.0,
        **kwargs,
    ) -> list[HexBlock]:
        """Build a C-grid block list from a body STL.

        Parameters
        ----------
        surface :
            Optional pre-built ``(ni_total, n_stations, 3)`` k=0 layer
            from :func:`shore.surface.cgrid.build_cgrid_k0`.  When
            supplied, *ni_body* / *ni_wake* / *n_stations* /
            *wake_length* are ignored.  When ``None``, the k=0 layer
            is built from *mesh* using the supplied parameters.

            For ``blocks=3`` the surface is split internally into per-
            block (upper / body / lower) k=0 layers via
            :meth:`shore.surface.cgrid.CGridK0Result.per_block_layers`;
            pre-built surfaces work equivalently as long as the
            ``ni_total = 2 * (ni_wake - 1) + ni_body`` invariant holds.
        nk :
            Total number of wall-normal layers (including k=0).
        ni_body :
            Per-station body wrap node count.  Required when *surface*
            is None.
        ni_wake :
            Per-branch wake node count (shared TE node included).
            Required when *surface* is None.
        n_stations :
            Number of span-axis stations.  Required when *surface* is
            None.  Becomes each block's ``nj``.
        wake_length :
            Length of each wake branch from the TE.  Required when
            *surface* is None.
        span_axis :
            ``"x"`` / ``"y"`` / ``"z"``.  Default ``"z"``.
        freestream_world :
            Length-3 vector giving the freestream direction in world
            coordinates.  ``None`` = default to chord-axis +1.
        te_dihedral_deg, wake_growth :
            Forwarded to :func:`build_cgrid_k0`.
        spacing_k :
            Wall-normal spacing law.  Defaults to
            ``Spacing1D(law='geometric', ds=1e-3)``.
        mesh :
            ``trimesh.Trimesh`` instance.  Required when *surface* is
            None.
        blocks :
            ``3`` (default, Xall-correct sharp-TE), ``4`` (Xall-correct
            blunt-TE: adds a wake-base block W bridging the TE
            thickness gap), or ``1`` (legacy single-block).  See class
            docstring for the trade-offs.
        nk_w :
            Wake-base block W's resolution across the TE thickness
            gap.  Used only for ``blocks=4``.  Default ``3`` is
            adequate for the thin gaps (≤ 1% chord) typical of
            blunt-TE airfoils.
        cap_tips :
            ``"none"`` (default), ``"jlo"``, ``"jhi"``, or ``"both"``:
            close those wing tips with cap blocks (ADR 0059). Each cap
            is an overset component: five blocks (``core`` and the ring
            ``upper`` / ``le`` / ``lower`` / ``te``) whose wall (``k_lo``)
            covers the tip face, extruded outwards; see
            :mod:`shore.topology.cgrid_cap`. The wing's first and last
            stations are placed on its tips (``TIP_STATION_MARGIN``), and
            the caps' edges are the tip stations' wall nodes, bit for bit.

            A sharp TE is rounded over the whole span by
            ``cap_te_thickness`` (or ``cap_te_thickness_frac * chord``), so
            the TE block has a base to span; ``blocks=3`` is then promoted
            to ``blocks=4`` (a stderr note). ``blocks=1`` with caps is
            rejected. A blunt TE is used as it is.
        cap_te_thickness :
            Absolute TE thickness for rounding a sharp TE. Mutually
            exclusive with ``cap_te_thickness_frac``.
        cap_te_thickness_frac :
            Same, as a fraction of the chord. Default ``0.002`` (0.2 %).
        cap_le_frac :
            Chord fraction of the LE cuts, which bound the ``le`` ring
            block. Default ``0.05``.
        cap_te_frac :
            Ignored since ADR 0059, with a ``DeprecationWarning`` when caps
            are built: the ``te`` block spans the TE base.
        cap_tip_inset :
            Distance of the caps' wall from the tip station, outwards
            along its normal. Default ``0``: the station is on the tip.
        cap_n_cap_j :
            Ring blocks' node count from the core to the section's wall.
            Default ``12``.
        cap_j_spacing, cap_j_beta :
            Spacing law of that direction (``"uniform"`` / ``"tanh"`` /
            ``"tanh2"``; ``"tanh"`` clusters at the wall) and its
            strength. Default ``"tanh"``, beta ``4.0``.
        cap_nk :
            Caps' node count along k (outwards from the tip face), at the
            wrap's wall-normal schedule *spacing_k*. Default ``nk``.
        cap_core_frac :
            The core's half-thickness as a fraction of the section's, in
            (0, 1). Default ``0.5``.
        cap_flare_deg :
            Outward lean of the ring's k-lines at the section's edge, in
            [0, 90) degrees, growing linearly from the core. Default
            ``30``.

        Returns
        -------
        list[HexBlock]
            For ``blocks=3``: ``[upper, body, lower]``.
            For ``blocks=4``: ``[upper, body, lower, wake_base]``.
            For ``blocks=1``: ``[cgrid]``.

            When ``cap_tips != "none"``: the 4-block C-grid followed
            by 5 cap blocks per requested tip: ``cap_<tip>_core``,
            ``cap_<tip>_upper``, ``cap_<tip>_le``, ``cap_<tip>_lower``,
            ``cap_<tip>_te``.
        """
        if blocks not in (1, 3, 4):
            raise ValueError(
                f"CGridTopology.build: `blocks` must be 1, 3, or 4; got {blocks}.  "
                "Use blocks=3 (default, sharp-TE), blocks=4 (blunt-TE), "
                "or blocks=1 (legacy)."
            )

        # ------------------------------------------------------------
        # Build / accept the k=0 layer + locate the TE columns.
        # ------------------------------------------------------------
        if surface is None:
            if mesh is None:
                raise ValueError(
                    "CGridTopology.build needs either a pre-built surface "
                    "or a trimesh.Trimesh `mesh=` argument."
                )
            for name, val in (
                ("ni_body", ni_body),
                ("ni_wake", ni_wake),
                ("n_stations", n_stations),
                ("wake_length", wake_length),
            ):
                if val is None:
                    raise ValueError(f"CGridTopology.build needs `{name}=` when surface is None.")
            k0_result = build_cgrid_k0(
                mesh,  # type: ignore[arg-type]
                ni_body=int(ni_body),  # type: ignore[arg-type]
                ni_wake=int(ni_wake),  # type: ignore[arg-type]
                n_stations=int(n_stations),  # type: ignore[arg-type]
                wake_length=float(wake_length),  # type: ignore[arg-type]
                span_axis=span_axis,
                freestream_world=freestream_world,
                te_dihedral_deg=te_dihedral_deg,
                wake_growth=wake_growth,
                end_margin=TIP_STATION_MARGIN if cap_tips != "none" else 0.02,
            )
            surface = k0_result.nodes  # (ni_total, n_stations, 3)
            _ni_body = k0_result.ni_body
            _ni_wake = k0_result.ni_wake
            _is_blunt = k0_result.te_result.is_blunt
        else:
            if ni_body is None or ni_wake is None:
                raise ValueError(
                    "CGridTopology.build with a pre-built `surface` requires "
                    "explicit `ni_body=` and `ni_wake=` so the wake / body "
                    "split point is unambiguous."
                )
            _ni_body = int(ni_body)
            _ni_wake = int(ni_wake)
            _is_blunt = False  # caller-supplied surface; trust them
            if surface.shape[0] != 2 * (_ni_wake - 1) + _ni_body:
                raise ValueError(
                    f"CGridTopology.build: surface i-extent {surface.shape[0]} "
                    f"does not match 2*(ni_wake-1) + ni_body = "
                    f"{2 * (_ni_wake - 1) + _ni_body}."
                )

        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

        # ------------------------------------------------------------
        # Conformal-cap path.  When cap_tips != "none" we may need to
        # auto-promote sharp 3-block to blunt 4-block + round the body
        # TE so the cap blocks share the airfoil-wall nodes bit-exactly
        # with the body wrap.
        # ------------------------------------------------------------
        if cap_tips not in ("none", "jlo", "jhi", "both"):
            raise ValueError(
                f"CGridTopology.build: cap_tips must be 'none' / 'jlo' / "
                f"'jhi' / 'both', got {cap_tips!r}"
            )
        if cap_te_frac is not None and cap_tips != "none":
            warnings.warn(
                "CGridTopology.build: cap_te_frac is ignored since ADR 0059 "
                "(the cap's TE block spans the TE base).",
                DeprecationWarning,
                stacklevel=2,
            )
        if cap_tips != "none":
            return self._build_with_caps(
                surface=surface,
                nk=nk,
                spacing_k=_sp_k,
                ni_body=_ni_body,
                ni_wake=_ni_wake,
                blocks=blocks,
                nk_w=int(nk_w),
                is_blunt=_is_blunt,
                span_axis=span_axis,
                cap_tips=cap_tips,
                cap_te_thickness=cap_te_thickness,
                cap_te_thickness_frac=cap_te_thickness_frac,
                cap_le_frac=cap_le_frac,
                cap_tip_inset=cap_tip_inset,
                cap_n_cap_j=cap_n_cap_j,
                cap_j_spacing=cap_j_spacing,
                cap_j_beta=cap_j_beta,
                cap_nk=cap_nk,
                cap_core_frac=cap_core_frac,
                cap_flare_deg=cap_flare_deg,
            )

        # ------------------------------------------------------------
        # No-caps dispatch by block count.
        # ------------------------------------------------------------
        if blocks == 1:
            return self._build_single_block(surface, nk, _sp_k)
        if blocks == 3:
            if _is_blunt:
                raise ValueError(
                    "CGridTopology.build: 3-block layout requires a sharp TE, "
                    "but the body's TE was detected as blunt.  Either close the "
                    "TE (e.g. shore profile naca --closed-te) or pass blocks=4 "
                    "to get the wake-base block W that meshes the TE gap."
                )
            return self._build_three_blocks(surface, nk, _sp_k, _ni_body, _ni_wake)
        # blocks == 4
        if not _is_blunt:
            raise ValueError(
                "CGridTopology.build: 4-block layout requires a blunt TE "
                "(non-zero TE thickness), but the body's TE was detected as "
                "sharp.  The wake-base block W has no thickness to bridge — "
                "use blocks=3 (default) for sharp-TE bodies."
            )
        return self._build_four_blocks(surface, nk, _sp_k, _ni_body, _ni_wake, nk_w=int(nk_w))

    def _build_with_caps(
        self,
        *,
        surface: np.ndarray,
        nk: int,
        spacing_k: Spacing1D,
        ni_body: int,
        ni_wake: int,
        blocks: int,
        nk_w: int,
        is_blunt: bool,
        span_axis: str,
        cap_tips: str,
        cap_te_thickness: float | None,
        cap_te_thickness_frac: float,
        cap_le_frac: float,
        cap_tip_inset: float | None,
        cap_n_cap_j: int,
        cap_j_spacing: str,
        cap_j_beta: float,
        cap_nk: int | None,
        cap_core_frac: float,
        cap_flare_deg: float,
    ) -> list[HexBlock]:
        """Build the conformal-cap C-grid: 4-block C-grid + cap blocks per tip.

        See :meth:`build` for the user-facing semantics of every parameter.
        """
        import sys

        from shore.surface.cgrid import round_body_te_in_surface
        from shore.topology.cgrid_cap import (
            CAP_BLOCKS,
            build_cap_blocks,
            check_cap_blocks,
            find_symmetric_le_cuts,
            wire_cap_blocks,
        )

        # ----- Reject 1-block + caps -----
        if blocks == 1:
            raise ValueError(
                "CGridTopology.build: cap_tips != 'none' is incompatible with "
                "blocks=1 (the legacy single-block layout cannot represent the "
                "conformal cap-body interface).  Use blocks=3 (auto-promoted "
                "to 4) or blocks=4 directly."
            )

        # ----- Resolve effective TE thickness -----
        # Derive a chord length from the body region of the surface.
        i_te_upper = ni_wake - 1
        i_te_lower = i_te_upper + ni_body - 1
        body_xs_mid = surface[i_te_upper : i_te_lower + 1, surface.shape[1] // 2, :]
        span_axis_world = "xyz".index(span_axis)
        extents = body_xs_mid.max(axis=0) - body_xs_mid.min(axis=0)
        extents[span_axis_world] = -np.inf
        chord_len = float(extents.max())
        if cap_te_thickness is not None:
            te = float(cap_te_thickness)
        else:
            te = float(cap_te_thickness_frac) * chord_len
        if te <= 0:
            raise ParameterError(
                f"CGridTopology.build: cap-tips TE thickness resolved to "
                f"non-positive value ({te:.6g})"
            )

        # ----- Body rounding (sharp body only) -----
        if not is_blunt:
            # Auto-promote 3-block → 4-block for sharp + caps.
            if blocks == 3:
                sys.stderr.write(
                    f"CGridTopology: cap_tips={cap_tips!r} on a sharp-TE body — "
                    f"auto-promoting blocks=3 to blocks=4 and rounding the "
                    f"body's TE by te_thickness={te:.6g} so the cap and body "
                    f"meshes match conformally.\n"
                )
                blocks = 4
            round_body_te_in_surface(
                surface=surface,
                ni_wake=ni_wake,
                ni_body=ni_body,
                te_thickness=te,
                span_axis=span_axis,
            )

        # ----- Build the 4-block C-grid (always 4-block under cap_tips) -----
        cgrid_blocks = self._build_four_blocks(
            surface=surface,
            nk=nk,
            spacing_k=spacing_k,
            ni_body=ni_body,
            ni_wake=ni_wake,
            nk_w=nk_w,
        )
        body = cgrid_blocks[1]

        # ----- Build the cap blocks per requested tip (ADR 0059) -----
        tips_to_build = ("jlo", "jhi") if cap_tips == "both" else (cap_tips,)
        cap_blocks: list[HexBlock] = []
        n_stations = body.nodes.shape[1]
        nk_cap = int(cap_nk) if cap_nk is not None else int(nk)
        distances = np.concatenate([[0.0], np.cumsum(spacing_k.build_steps(nk_cap))])
        tip_inset = 0.0 if cap_tip_inset is None else float(cap_tip_inset)

        for tip in tips_to_build:
            j_extreme = 0 if tip == "jlo" else n_stations - 1
            j_inboard = 1 if tip == "jlo" else n_stations - 2
            # The tip station's wall nodes, TE rounded above: the cap's edge.
            cross_section = body.nodes[:, j_extreme, 0, :].copy()
            outward = cross_section.mean(axis=0) - body.nodes[:, j_inboard, 0, :].mean(axis=0)
            if float(np.linalg.norm(outward)) < 1e-12:
                raise ParameterError(
                    f"CGridTopology.build: at tip {tip!r}, the body's adjacent "
                    "spanwise stations have coincident centroids; cannot infer "
                    "the outward span direction."
                )
            le_up, le_lo = find_symmetric_le_cuts(
                cross_section, le_frac=cap_le_frac, span_axis_world=span_axis_world
            )
            cap_arrays = build_cap_blocks(
                cross_section,
                outward=outward,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                k_distances=distances,
                n_cap_j=cap_n_cap_j,
                j_spacing=cap_j_spacing,
                j_beta=cap_j_beta,
                core_frac=cap_core_frac,
                flare_deg=cap_flare_deg,
                tip_inset=tip_inset,
            )
            tip_blocks = [
                HexBlock(
                    label=f"cap_{tip}_{name}",
                    nodes=cap_arrays[name],
                    periodic_i=False,
                    periodic_j=False,
                    spacing_k=spacing_k,
                )
                for name in CAP_BLOCKS
            ]
            wire_cap_blocks(tip_blocks)
            cap_blocks.extend(tip_blocks)

        check_cap_blocks(cap_blocks)
        return cgrid_blocks + cap_blocks

    # ------------------------------------------------------------------
    # Internal block builders
    # ------------------------------------------------------------------

    def _build_single_block(
        self,
        surface: np.ndarray,
        nk: int,
        spacing_k: Spacing1D,
    ) -> list[HexBlock]:
        """Legacy 1-block C-grid (single hex block, wake cut as i-overlap)."""
        ni_total, nj, _ = surface.shape
        nodes = np.zeros((ni_total, nj, nk, 3), dtype=np.float64)
        nodes[:, :, 0, :] = surface
        block = HexBlock(
            label="cgrid",
            nodes=nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )
        return [block]

    def _build_three_blocks(
        self,
        surface: np.ndarray,
        nk: int,
        spacing_k: Spacing1D,
        ni_body: int,
        ni_wake: int,
    ) -> list[HexBlock]:
        """3-block C-grid (upper wake / body / lower wake) with SHARED seams."""
        nj = surface.shape[1]
        # The full (ni_total, nj, nk, 3) array.  We allocate one
        # contiguous buffer and carve three views from it so that
        # node coincidence at the upper-TE and lower-TE columns
        # comes for free via numpy view aliasing.
        ni_total = 2 * (ni_wake - 1) + ni_body
        nodes_all = np.zeros((ni_total, nj, nk, 3), dtype=np.float64)
        nodes_all[:, :, 0, :] = surface

        i_te_upper = ni_wake - 1
        i_te_lower = i_te_upper + ni_body - 1
        upper_nodes = nodes_all[: i_te_upper + 1, :, :, :]
        body_nodes = nodes_all[i_te_upper : i_te_lower + 1, :, :, :]
        lower_nodes = nodes_all[i_te_lower:, :, :, :]

        upper = HexBlock(
            label="upper",
            nodes=upper_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )
        body = HexBlock(
            label="body",
            nodes=body_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )
        lower = HexBlock(
            label="lower",
            nodes=lower_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )

        # ------------------------------------------------------------
        # Wire SHARED seams.
        # ------------------------------------------------------------
        # 1. upper.i_hi  ↔ body.i_lo  (upper-TE column, identity).
        upper._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=body.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        body._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=upper.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )

        # 2. body.i_hi   ↔ lower.i_lo (lower-TE column, identity).
        body._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=body.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )

        # 3. upper.k_lo  ↔ lower.k_lo (wake cut).
        # Same axis (k↔k); along-seam axes (i, j); i flipped, j identity.
        wake_cut_axis_map = AxisMap.from_pairs(("i", True), ("j", False))
        upper._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=lower.face("k_lo"),
            axis_map=wake_cut_axis_map,
        )
        lower._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=upper.face("k_lo"),
            axis_map=wake_cut_axis_map.inverse("k_lo", "k_lo"),
        )

        return [upper, body, lower]

    def _build_four_blocks(
        self,
        surface: np.ndarray,
        nk: int,
        spacing_k: Spacing1D,
        ni_body: int,
        ni_wake: int,
        *,
        nk_w: int,
    ) -> list[HexBlock]:
        """4-block C-grid for blunt TE: 3-block layout + wake-base block W.

        The wake-base block W bridges the TE thickness gap between
        the upper and lower wake k=0 strips.  W's k axis runs from
        ``y = -te_thickness/2`` (= L's k=0 line) to ``y = +te_thickness/2``
        (= U's k=0 line).  W's i axis matches U's (far-tip → TE):

        - ``W.i_lo`` (far-tip, *FREE*).
        - ``W.i_hi`` (TE base wall, *WALL*).  The bluff-base of the
          airfoil is meshed by this face — the structural fix the
          single-block C-grid couldn't represent.
        - ``W.j_lo`` / ``W.j_hi`` (wing tips, *FREE*).
        - ``W.k_lo`` (lower face, ``SHARED`` with L.k_lo).  Same-axis
          (k↔k); along-axes (i, j); i flipped because L's i runs the
          opposite direction (TE → far-tip).
        - ``W.k_hi`` (upper face, ``SHARED`` with U.k_lo).  Same-axis;
          identity AxisMap (W's i and U's i match).

        W is **not marched** — its geometry is fully determined at
        build time by linear interpolation in k between the lower-
        wake and upper-wake k=0 strips.  The marcher's per-block
        loop must skip W (see :meth:`Mesh._march_cgrid_three_or_four`).
        """
        if nk_w < 2:
            raise ValueError(f"_build_four_blocks: nk_w must be >= 2, got {nk_w}")

        # ------------------------------------------------------------
        # Build the 3 marching blocks (upper, body, lower) the usual
        # way — same SHARED seams at the upper-TE and lower-TE
        # columns.  The wake-cut SHARED face wiring differs: in the
        # 4-block layout, U.k_lo partners W.k_hi (not L.k_lo) and
        # L.k_lo partners W.k_lo (not U.k_lo).
        # ------------------------------------------------------------
        nj = surface.shape[1]
        ni_total = 2 * (ni_wake - 1) + ni_body
        nodes_all = np.zeros((ni_total, nj, nk, 3), dtype=np.float64)
        nodes_all[:, :, 0, :] = surface

        i_te_upper = ni_wake - 1
        i_te_lower = i_te_upper + ni_body - 1
        upper_nodes = nodes_all[: i_te_upper + 1, :, :, :]
        body_nodes = nodes_all[i_te_upper : i_te_lower + 1, :, :, :]
        lower_nodes = nodes_all[i_te_lower:, :, :, :]

        upper = HexBlock(
            label="upper",
            nodes=upper_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )
        body = HexBlock(
            label="body",
            nodes=body_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )
        lower = HexBlock(
            label="lower",
            nodes=lower_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )

        # ------------------------------------------------------------
        # Build the wake-base block W.
        # ------------------------------------------------------------
        # W's i-axis matches U's (far-tip → TE).  W's k-axis spans
        # the TE thickness from the lower-wake k=0 strip up to the
        # upper-wake k=0 strip.
        upper_k0 = upper_nodes[:, :, 0, :]  # (ni_wake, nj, 3) at y = +te/2
        lower_k0 = lower_nodes[:, :, 0, :]  # (ni_wake, nj, 3) at y = -te/2
        # Reverse lower_k0 along i so W's i orientation matches upper.
        lower_k0_rev = lower_k0[::-1, :, :]

        # Linear interpolation across the gap.
        ts = np.linspace(0.0, 1.0, nk_w)  # k=0 -> lower; k=last -> upper
        wake_base_nodes = np.zeros((ni_wake, nj, nk_w, 3), dtype=np.float64)
        for k_idx, t in enumerate(ts):
            wake_base_nodes[:, :, k_idx, :] = (1.0 - t) * lower_k0_rev + t * upper_k0

        wake_base = HexBlock(
            label="wake_base",
            nodes=wake_base_nodes,
            periodic_i=False,
            periodic_j=False,
            spacing_k=spacing_k,
        )

        # ------------------------------------------------------------
        # Wire SHARED seams (4-block topology).
        # ------------------------------------------------------------
        # 1-2. Upper-TE and lower-TE columns — identical to 3-block.
        upper._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=body.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        body._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=upper.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        body._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=body.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )

        # 3. W.k_hi ↔ upper.k_lo (top of W bridges to upper-wake k=0).
        # Same-axis (k↔k), identity AxisMap on (i, j) — W's i and
        # upper's i both run far-tip → TE.
        wake_base._set_face_bc(
            "k_hi",
            FaceBC.SHARED,
            partner=upper.face("k_lo"),
            axis_map=AxisMap.identity("k_hi"),
        )
        upper._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=wake_base.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )

        # 4. W.k_lo ↔ lower.k_lo (bottom of W bridges to lower-wake k=0).
        # Same-axis (k↔k); along-axes (i, j); i flipped because lower's
        # i runs TE → far-tip while W's runs far-tip → TE.
        i_flip_axis_map = AxisMap.from_pairs(("i", True), ("j", False))
        wake_base._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=lower.face("k_lo"),
            axis_map=i_flip_axis_map,
        )
        lower._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=wake_base.face("k_lo"),
            axis_map=i_flip_axis_map.inverse("k_lo", "k_lo"),
        )

        # ------------------------------------------------------------
        # W's i_hi face is the airfoil's TE-base WALL.  The default
        # FaceBC at construction is FREE; promote it to WALL so
        # downstream cc.par writers see the bluff-base.
        # ------------------------------------------------------------
        wake_base._set_face_bc("i_hi", FaceBC.WALL)

        return [upper, body, lower, wake_base]


class CGridTipCap(Topology):
    """A wing-tip cap from one tip station of a C-grid, as its own component.

    Builds the cap of :mod:`shore.topology.cgrid_cap` (ADR 0059) from a
    cross-section polyline, the C-grid body's wall nodes at its tip station:
    five blocks ``cap_core``, ``cap_upper``, ``cap_le``, ``cap_lower``,
    ``cap_te``. Their wall (``k_lo``) covers the section, moved outwards by
    ``tip_inset`` onto the tip face; k leaves it along the section's normal.
    The faces between cap blocks are SHARED; the others are FREE, so the cap
    meets the C-grid wrap as an overset component (``shore cc-par``).

    The cap is volumetric at build time: :meth:`shore.mesh.Mesh.march` on
    a CGridTipCap mesh is a no-op.

    Parameters
    ----------
    le_cut_upper_idx, le_cut_lower_idx :
        i-indices of the LE cuts in the cross-section, symmetric about its
        middle (``find_symmetric_le_cuts``). They bound the ``le`` block.
    span_dir :
        Length-3 unit vector pointing **outward** along the span at this
        tip: it picks the side of the section's plane the cap grows on.
    tip_inset :
        Distance of the cap's wall from the cross-section, outwards along
        its normal, >= 0: the distance from the station to the tip face.
    n_cap_j :
        Ring blocks' node count from the core to the section's wall, >= 2.
    cap_j_spacing, cap_j_beta :
        That direction's spacing law (``"uniform"`` / ``"tanh"`` /
        ``"tanh2"``; ``"tanh"`` clusters at the wall) and strength.
    is_blunt, wake_base_cross_section, te_cut_upper_idx, te_cut_lower_idx :
        Ignored since ADR 0059 (a ``DeprecationWarning`` when given): the
        ``te`` block spans the TE base, and nothing closes the wake base.
    auto_round_te, te_thickness, te_thickness_frac :
        A sharp TE is rounded, in the cap only, to ``te_thickness`` (or
        ``te_thickness_frac`` of the chord, default 0.002) so the ``te``
        block has a base; with ``auto_round_te=False`` it is refused.
    core_frac :
        The core's half-thickness as a fraction of the section's, in (0, 1).
    flare_deg :
        Outward lean of the ring's k-lines at the section's edge, in [0, 90).
    k_distances :
        Distances of the k layers from the wall, ``k_distances[0] == 0``,
        one per node. Default: *spacing_k*'s steps for ``nk`` nodes.
    """

    def __init__(
        self,
        le_cut_upper_idx: int,
        le_cut_lower_idx: int,
        span_dir: tuple[float, float, float],
        tip_inset: float = 0.0,
        n_cap_j: int = 12,
        cap_j_spacing: str = "tanh",
        cap_j_beta: float = 4.0,
        is_blunt: bool = False,
        wake_base_cross_section: np.ndarray | None = None,
        te_cut_upper_idx: int | None = None,
        te_cut_lower_idx: int | None = None,
        auto_round_te: bool = True,
        te_thickness: float | None = None,
        te_thickness_frac: float | None = 0.002,
        core_frac: float = 0.5,
        flare_deg: float = 30.0,
        k_distances: np.ndarray | None = None,
    ) -> None:
        span_arr = np.asarray(span_dir, dtype=np.float64)
        if span_arr.shape != (3,):
            raise ParameterError(
                f"CGridTipCap: span_dir must have shape (3,); got {span_arr.shape}"
            )
        span_norm = float(np.linalg.norm(span_arr))
        if not np.isclose(span_norm, 1.0, atol=1e-9):
            raise ParameterError(
                f"CGridTipCap: span_dir must be a unit vector; got "
                f"|span_dir|={span_norm:.6g}.  Normalise before passing."
            )
        if tip_inset < 0:
            raise ParameterError(f"CGridTipCap: tip_inset must be >= 0, got {tip_inset}")
        if n_cap_j < 2:
            raise ParameterError(f"CGridTipCap: n_cap_j must be >= 2, got {n_cap_j}")
        if cap_j_spacing not in ("uniform", "tanh", "tanh2"):
            raise ParameterError(
                f"CGridTipCap: cap_j_spacing must be 'uniform' / 'tanh' / "
                f"'tanh2', got {cap_j_spacing!r}"
            )
        if te_thickness is not None and te_thickness_frac is not None:
            raise ParameterError(
                "CGridTipCap: pass either te_thickness (absolute) or "
                "te_thickness_frac (fraction of chord), not both."
            )
        if te_thickness is not None and te_thickness <= 0:
            raise ParameterError(f"CGridTipCap: te_thickness must be > 0, got {te_thickness}")
        if te_thickness_frac is not None and te_thickness_frac <= 0:
            raise ParameterError(
                f"CGridTipCap: te_thickness_frac must be > 0, got {te_thickness_frac}"
            )
        if (
            is_blunt
            or wake_base_cross_section is not None
            or te_cut_upper_idx is not None
            or te_cut_lower_idx is not None
        ):
            warnings.warn(
                "CGridTipCap: is_blunt, wake_base_cross_section and te_cut_*_idx are "
                "ignored since ADR 0059 (the te block spans the TE base).",
                DeprecationWarning,
                stacklevel=2,
            )
        self.le_cut_upper_idx = int(le_cut_upper_idx)
        self.le_cut_lower_idx = int(le_cut_lower_idx)
        self.span_dir = span_arr
        self.tip_inset = float(tip_inset)
        self.n_cap_j = int(n_cap_j)
        self.cap_j_spacing = cap_j_spacing
        self.cap_j_beta = float(cap_j_beta)
        self.auto_round_te = bool(auto_round_te)
        self.te_thickness = None if te_thickness is None else float(te_thickness)
        self.te_thickness_frac = None if te_thickness_frac is None else float(te_thickness_frac)
        self.core_frac = float(core_frac)
        self.flare_deg = float(flare_deg)
        self.k_distances = None if k_distances is None else np.asarray(k_distances, float)

    def build(
        self,
        surface: np.ndarray,
        nk: int,
        spacing_k: Spacing1D | None = None,
        **kwargs,
    ) -> list[HexBlock]:
        """Build the cap blocks from a cross-section polyline.

        Parameters
        ----------
        surface :
            ``(ni_body, 3)`` cross-section polyline at the wing tip, from
            the upper TE round the LE to the lower TE.
        nk :
            The cap's node count along k (outwards). Must be >= 2, and equal
            to ``len(k_distances)`` when those are given.
        spacing_k :
            Wall-normal spacing law of the layers when ``k_distances`` is
            not given. Default ``Spacing1D(law='geometric', ds=1e-3)``.

        Returns
        -------
        list[HexBlock]
            ``[cap_core, cap_upper, cap_le, cap_lower, cap_te]``.
        """
        import sys

        from shore.topology.cgrid_cap import (
            CAP_BLOCKS,
            auto_round_sharp_te,
            build_cap_blocks,
            check_cap_blocks,
            detect_te_thickness,
            wire_cap_blocks,
        )

        if nk < 2:
            raise ParameterError(f"CGridTipCap.build: nk must be >= 2, got {nk}")
        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)
        if self.k_distances is not None:
            if self.k_distances.shape != (nk,):
                raise ParameterError(
                    f"CGridTipCap.build: {self.k_distances.size} k_distances for nk={nk}"
                )
            distances = self.k_distances
        else:
            distances = np.concatenate([[0.0], np.cumsum(_sp_k.build_steps(int(nk)))])

        # A sharp TE is rounded in the cap only, so the te block has a base.
        span_axis_world = int(np.argmax(np.abs(self.span_dir)))
        cs_for_cap = np.asarray(surface, dtype=np.float64)
        if detect_te_thickness(cs_for_cap, span_axis_world) < 1e-9:
            if not self.auto_round_te:
                raise ParameterError(
                    "CGridTipCap.build: cross-section has a sharp TE "
                    "(i=0 and i=ni-1 coincide) but auto_round_te=False; "
                    "the cap's te block needs a TE base.  Either enable "
                    "auto_round_te (default) or pre-round the TE in the "
                    "source body."
                )
            if self.te_thickness is not None:
                te = self.te_thickness
            else:
                extents = cs_for_cap.max(axis=0) - cs_for_cap.min(axis=0)
                extents[span_axis_world] = -np.inf
                frac = self.te_thickness_frac if self.te_thickness_frac is not None else 0.002
                te = frac * float(extents.max())
            cs_for_cap = auto_round_sharp_te(
                cs_for_cap, te_thickness=te, span_axis_world=span_axis_world
            )
            sys.stderr.write(
                f"CGridTipCap: sharp TE detected; auto-rounding cap with "
                f"te_thickness={te:.6g} (the body wrap is unchanged).\n"
            )

        arrays = build_cap_blocks(
            cs_for_cap,
            outward=self.span_dir,
            le_cut_upper_idx=self.le_cut_upper_idx,
            le_cut_lower_idx=self.le_cut_lower_idx,
            k_distances=distances,
            n_cap_j=self.n_cap_j,
            j_spacing=self.cap_j_spacing,
            j_beta=self.cap_j_beta,
            core_frac=self.core_frac,
            flare_deg=self.flare_deg,
            tip_inset=self.tip_inset,
        )
        blocks = [
            HexBlock(
                label=f"cap_{name}",
                nodes=arrays[name],
                periodic_i=False,
                periodic_j=False,
                spacing_k=_sp_k,
            )
            for name in CAP_BLOCKS
        ]
        wire_cap_blocks(blocks)
        check_cap_blocks(blocks)
        return blocks


class OHTopology(Topology):
    """Style-2 conforming hybrid: cubed-sphere wrap + 6-block H-channel.

    A 12-block decomposition: 6 wrap blocks (cubed-sphere, marched outward
    from the STL surface to ``nk_wrap - 1``) and 6 channel blocks (frustums
    bridging the wrap's outer envelope to the user-supplied bounding box).
    Wrap and channel are SHARED-by-view at the seam — they share one
    contiguous numpy buffer per pair, so ``wrap.k_hi`` and ``chan.k_lo`` are
    bit-identical.

    The wrap is built with ``phi_offset = pi/4`` so the four equatorial
    sub-block *centres* (rather than corners) face the world cardinal axes:
    sub0→+y, sub1→-x, sub2→-y, sub3→+x.  Combined with the two cap blocks
    (cap_south→+z, cap_north→-z), each of the 6 wrap blocks faces exactly
    one of the box's 6 faces, allowing a 1:1 wrap-block-to-channel-block
    mapping.

    Block layout returned by :meth:`build`:

        [sub0, sub1, sub2, sub3, cap_north, cap_south,
         chan_sub0, chan_sub1, chan_sub2, chan_sub3, chan_north, chan_south]

    Channel construction is deferred until after marching: the wrap's k_hi
    envelope is unknown at build time, so the channel's k=0 row (which IS
    the wrap's k_hi via view aliasing) gets filled when the wrap finishes
    marching.  :meth:`shore.mesh.Mesh.march` calls
    :meth:`fill_channel_blocks` after the wrap march completes; users who
    drive marching manually must call :meth:`fill_channel_blocks` before
    using the channel block geometry.

    Default face BCs:

    - Wrap k_lo = WALL (STL surface).
    - Wrap k_hi = SHARED (with channel k_lo).
    - Channel k_lo = SHARED (with wrap k_hi).
    - Channel k_hi = FREE (far-field box face).
    - Wrap and channel lateral seams (j_hi↔j_lo, i_hi↔cap, i_lo↔cap)
      mirror the cubed-sphere topology on both sides.
    """

    def __init__(
        self,
        box_vmin: tuple[float, float, float],
        box_vmax: tuple[float, float, float],
        nk_channel: int,
        k_law_channel: str = "uniform",
        k_beta_channel: float = 3.0,
        c1_seam: bool = True,
    ) -> None:
        if nk_channel < 2:
            raise ParameterError(f"OHTopology: nk_channel must be >= 2, got {nk_channel}")
        if k_law_channel not in ("uniform", "tanh", "tanh2"):
            raise ParameterError(
                f"OHTopology: k_law_channel must be 'uniform' / 'tanh' / "
                f"'tanh2', got {k_law_channel!r}"
            )
        vmin = np.asarray(box_vmin, dtype=np.float64)
        vmax = np.asarray(box_vmax, dtype=np.float64)
        if vmin.shape != (3,) or vmax.shape != (3,):
            raise ParameterError(
                f"OHTopology: box_vmin and box_vmax must be length-3, got "
                f"{vmin.shape} and {vmax.shape}"
            )
        if not np.all(vmax > vmin):
            raise ParameterError(
                f"OHTopology: box_vmin must be strictly less than box_vmax "
                f"on every axis; got box_vmin={tuple(vmin)} box_vmax={tuple(vmax)}"
            )
        self.box_vmin = vmin
        self.box_vmax = vmax
        self.nk_channel = int(nk_channel)
        self.k_law_channel = k_law_channel
        self.k_beta_channel = float(k_beta_channel)
        self.c1_seam = bool(c1_seam)

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
        """Build the 12-block OH topology.

        ``nk`` is the wrap depth (number of layers including k=0).  The
        channel depth comes from the OHTopology constructor.  Returns
        ``[sub0..sub3, cap_north, cap_south, chan_sub0..chan_sub3,
        chan_north, chan_south]``.
        """
        if ni is None or nj is None:
            if surface is None:
                raise ParameterError("OHTopology.build needs ni and nj (or a surface array).")
            ni = ni if ni is not None else int(surface.shape[0])
            nj = nj if nj is not None else int(surface.shape[1])
        nj = _adjust_nj(nj)

        cap_proj = (
            projector
            if projector is not None
            else (_projector_from_mesh(mesh) if mesh is not None else None)
        )
        if cap_proj is None:
            raise ParameterError(
                "OHTopology.build requires a Projector (pass projector= or mesh=)."
            )

        nk_wrap = nk
        nk_chan = self.nk_channel
        total_k = nk_wrap + nk_chan - 1  # shared seam slot at index nk_wrap-1

        # OH pins theta_cap = pi/4 so each sub-block envelope maps cleanly
        # onto exactly one lateral cube face under central projection.  At
        # any other theta_cap the cap-equator seam falls inside a cube
        # face, which would force the channel outer faces to be L-shaped
        # across cube edges (giving up structured-block planarity).  The
        # caller's theta_cap is ignored here — cube geometry dictates it.
        theta_cap_oh = 0.25 * np.pi
        wrap_blocks = _build_cubed_sphere_blocks(
            ni=ni,
            nj=nj,
            nk=nk_wrap,
            total_k=total_k,
            cap_proj=cap_proj,
            spacing_k=spacing_k,
            theta_cap=theta_cap_oh,
            i_spacing=i_spacing,
            i_beta=i_beta,
            j_spacing=j_spacing,
            j_beta=j_beta,
            phi_offset=0.25 * np.pi,
        )

        # Carve channel views out of the same buffers.
        channel_blocks: list[HexBlock] = []
        for q, wrap_blk in enumerate(wrap_blocks):
            full_nodes = wrap_blk._oh_full_nodes  # type: ignore[attr-defined]
            # Channel view covers [nk_wrap - 1 : total_k] on the k axis,
            # which has length nk_chan and includes the shared seam at k=0.
            chan_view = full_nodes[:, :, nk_wrap - 1 :, :]
            assert chan_view.shape[2] == nk_chan
            label = f"chan_sub{q}" if q < 4 else "chan_north" if q == 4 else "chan_south"
            chan_blk = HexBlock(
                label=label,
                nodes=chan_view,
                periodic_i=False,
                periodic_j=False,
                spacing_k=spacing_k or Spacing1D(law="geometric", ds=1e-3),
            )
            channel_blocks.append(chan_blk)

        # Wire wrap.k_hi <-> channel.k_lo (SHARED, view-aliased, identity).
        for w_blk, c_blk in zip(wrap_blocks, channel_blocks, strict=True):
            w_blk._set_face_bc(
                "k_hi",
                FaceBC.SHARED,
                partner=c_blk.face("k_lo"),
                axis_map=AxisMap.identity("k_hi"),
            )
            c_blk._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=w_blk.face("k_hi"),
                axis_map=AxisMap.identity("k_lo"),
            )

        # Wire channel-channel lateral seams to mirror the wrap topology.
        # (sub_q.j_hi ↔ sub_(q+1).j_lo for q in 0..3 → same on chan side)
        chan_subs = channel_blocks[:4]
        for q in range(4):
            q_next = (q + 1) % 4
            chan_subs[q]._set_face_bc(
                "j_hi",
                FaceBC.SHARED,
                partner=chan_subs[q_next].face("j_lo"),
                axis_map=AxisMap.identity("j_hi"),
            )
            chan_subs[q_next]._set_face_bc(
                "j_lo",
                FaceBC.SHARED,
                partner=chan_subs[q].face("j_hi"),
                axis_map=AxisMap.identity("j_lo"),
            )

        # Wire chan_sub.i_lo/i_hi ↔ chan_cap.<face> (DIRICHLET, mirror wrap).
        chan_north_blk = channel_blocks[4]
        chan_south_blk = channel_blocks[5]
        sub_face_lit: FaceName
        for sub_face_lit, chan_cap_blk, table in (
            ("i_hi", chan_north_blk, _CAP_EQUATOR_NORTH),
            ("i_lo", chan_south_blk, _CAP_EQUATOR_SOUTH),
        ):
            for q, (cap_face, sub_axis_map) in enumerate(table):
                csub = chan_subs[q]
                chan_cap_blk._set_face_bc(
                    cap_face,
                    FaceBC.DIRICHLET,
                    partner=csub.face(sub_face_lit),
                    axis_map=sub_axis_map.inverse(sub_face_lit, cap_face),
                )
                csub._set_face_bc(
                    sub_face_lit,
                    FaceBC.DIRICHLET,
                    partner=chan_cap_blk.face(cap_face),
                    axis_map=sub_axis_map,
                )

        # k_hi on every channel block defaults to FREE (HexBlock constructor).
        # No further wiring needed.

        return [*wrap_blocks, *channel_blocks]

    def fill_channel_blocks(self, blocks: list[HexBlock]) -> None:
        """Build channel block interiors after the wrap has been marched.

        Walks the 6 wrap-channel pairs, computes each channel's outer face
        on the user-supplied bounding box, and fills the channel's
        ``k=1..nk_channel-1`` interior by ruled-surface frustum interpolation.

        The wrap-channel seam (channel k=0 = wrap k_hi) is already populated
        because the two faces share memory.  Only k=1..nk_channel-1 needs
        filling.

        When ``self.c1_seam`` is true (default), the channel's first step
        at every (i, j) is pinned to the wrap's last marching step at the
        same (i, j) — bit-exact C1 spacing match across the seam.  Pass
        ``c1_seam=False`` to the constructor to disable and use the raw
        spacing law instead (cell-size jump at the seam).
        """
        wrap_blocks = blocks[:6]
        channel_blocks = blocks[6:12]

        # Reference centre for the radial "which face does this block point
        # at?" decision and the corner sidedness snap.  This is the global
        # wrap-envelope centroid (mean over all six k_hi faces) — a robust,
        # purely geometric proxy for the body centre.  It is deliberately
        # NOT the box centre: the body need not sit at the box centre (e.g.
        # a sphere placed upstream in a long downstream-biased wake box), and
        # using the box centre there makes every block snap toward the
        # nearest box face instead of the face it actually points at.
        body_centre = np.stack(
            [w.nodes[:, :, -1, :].reshape(-1, 3).mean(axis=0) for w in wrap_blocks],
            axis=0,
        ).mean(axis=0)

        for w_blk, c_blk in zip(wrap_blocks, channel_blocks, strict=True):
            inner = w_blk.nodes[:, :, -1, :]  # (ni, nj, 3) — wrap envelope
            outer_quad = self._outer_quad_for_block(inner, body_centre)
            if self.c1_seam:
                # Per-(i, j) last marching step on the wrap.  Smoothing means
                # this isn't perfectly uniform, so we pin per-node.
                last_step = np.linalg.norm(
                    w_blk.nodes[:, :, -1, :] - w_blk.nodes[:, :, -2, :],
                    axis=-1,
                )  # (ni, nj)
                c1 = last_step
            else:
                c1 = None
            frustum = build_channel_block(
                inner_face=inner,
                outer_quad=outer_quad,
                nk_channel=self.nk_channel,
                k_law=self.k_law_channel,
                k_beta=self.k_beta_channel,
                c1_first_step=c1,
            )
            # Write only the k>=1 interior; k=0 is the shared seam.
            c_blk.nodes[:, :, 1:, :] = frustum[:, :, 1:, :]

    def _outer_quad_for_block(self, inner_face: np.ndarray, body_centre: np.ndarray) -> np.ndarray:
        """Find the box face this wrap block points at; return its 4 corners.

        ``body_centre`` is the global wrap-envelope centroid (see
        :meth:`fill_channel_blocks`).  Everything here is judged relative to
        it, never to the box centre, so the body may sit anywhere inside the
        box (including far off-centre, e.g. a body placed upstream in a long
        downstream wake box).

        Algorithm: pick the cube face from the block's **outward direction**
        — the vector from ``body_centre`` to this block's envelope centroid.
        The face is the one whose outward normal that direction most aligns
        with (largest signed component), NOT the geometrically nearest box
        plane.  For an off-centre body the nearest plane is often the wrong
        one (every block's nearest face is whichever box wall the body was
        pushed toward), so a distance-based pick collapses multiple blocks
        onto a single face.  Then snap each of the 4 envelope corners onto
        the cube by:

        1. Setting the face-normal component to the picked plane value.
        2. For each in-plane component, taking the sign of the corner's
           offset from ``body_centre`` and pushing it to the corresponding
           ``vmin`` / ``vmax``.

        The result: each corner lands at a cube **vertex**, the 4 corners
        of a single cube face.  Adjacent wrap blocks pick adjacent faces
        (sub0→+y, sub1→-x, ...); their shared envelope seams sit at the
        same cube *edge* under the snap, so the channel-channel lateral
        seams are bit-exact across the SHARED-by-view memory aliases.

        The (a=0,b=0)..(a=0,b=1) corner ordering matches the wrap's (i, j)
        winding so that the wrap-channel SHARED seam is identity-mapped
        and the channel-channel lateral seams are identity-mapped too.
        """
        c00 = inner_face[0, 0, :]
        c10 = inner_face[-1, 0, :]
        c11 = inner_face[-1, -1, :]
        c01 = inner_face[0, -1, :]

        full_centroid = inner_face.reshape(-1, 3).mean(axis=0)
        if not np.all((full_centroid >= self.box_vmin) & (full_centroid <= self.box_vmax)):
            raise ParameterError(
                f"OHTopology: wrap envelope centroid {full_centroid} is "
                f"outside the bounding box [{self.box_vmin}, {self.box_vmax}].  "
                f"The box must contain the marched wrap envelope; either "
                f"enlarge the box or reduce nk / wrap step sizes."
            )

        # Outward direction of this block, relative to the body centre.  The
        # six cubed-sphere blocks point at the six cardinal directions; the
        # dominant signed component selects the box face.
        outward = full_centroid - body_centre
        signed = np.array(
            [
                outward[0],  # +x face
                -outward[0],  # -x face
                outward[1],  # +y face
                -outward[1],  # -y face
                outward[2],  # +z face
                -outward[2],  # -z face
            ]
        )
        face_idx = int(np.argmax(signed))
        axis_idx = face_idx // 2
        plane_val = (self.box_vmax if face_idx % 2 == 0 else self.box_vmin)[axis_idx]

        out_corners = np.stack([c00, c10, c11, c01], axis=0).copy()
        out_corners[:, axis_idx] = plane_val
        for ax in range(3):
            if ax == axis_idx:
                continue
            # Push each corner's in-plane component to vmin or vmax based on
            # its sign relative to the body centre.  This snaps every corner
            # to a cube vertex of the picked face.
            sign = np.sign(out_corners[:, ax] - body_centre[ax])
            # Treat exactly-on-axis corners (sign == 0) as positive; in
            # practice they don't occur for cubed-sphere envelopes but the
            # snap must be deterministic.
            sign = np.where(sign == 0, 1.0, sign)
            out_corners[:, ax] = np.where(sign > 0, self.box_vmax[ax], self.box_vmin[ax])
        return out_corners


class CHTopology(Topology):
    """Hybrid C-grid wrap + H-channel far-field for sharp- or blunt-TE bodies.

    Mirrors :class:`OHTopology` but with a C-grid wrap (3 or 4 blocks)
    instead of a cubed-sphere wrap.  The user-supplied bounding "box"
    is a rectangle in the chord/normal plane (extruded across the wrap's
    span extents), defined by ``chord_min``, ``chord_max``, ``normal_min``,
    ``normal_max`` in world coords.  ``freestream_world`` and
    ``span_axis`` select which world axes carry chord and span.

    **Block layout** returned by :meth:`build`:

    CH-3 (sharp TE) — 12 blocks::

        [upper, body, lower,                              # wrap
         chan_upper, chan_body_top, chan_body_front,
         chan_body_bot, chan_lower,                       # k_hi-side channels
         chan_outflow_top_upper, chan_outflow_wake_upper,
         chan_outflow_wake_lower, chan_outflow_top_lower] # outflow (4 sub-blocks)

    CH-4 (blunt TE) — 14 blocks::

        [upper, body, lower, wake_base,                   # wrap
         chan_upper, chan_body_top, chan_body_front,
         chan_body_bot, chan_lower,                       # k_hi-side channels
         chan_outflow_top_upper, chan_outflow_wake_upper,
         chan_outflow_wake_lower, chan_outflow_top_lower,
         chan_outflow_wake_base]                          # outflow (5 sub-blocks)

    The outflow region is partitioned per side into a *top* sub-block
    (conformal SHARED to ``chan_upper.i_lo`` / ``chan_lower.i_hi`` —
    same ``nk_channel``) and a *wake* sub-block (conformal SHARED to
    the wrap's ``upper.i_lo`` / ``lower.i_hi`` — same ``nk_wrap``).
    Adjacent sub-blocks meet at envelope-far-tip seams (SHARED) and
    at the wake-cut seam (SHARED across upper and lower wake halves;
    via ``chan_outflow_wake_base`` for blunt TE).

    The 5 ``k_hi``-side channels (``chan_upper``, ``chan_body_top``,
    ``chan_body_front``, ``chan_body_bot``, ``chan_lower``) are
    SHARED-by-view with their wrap partners — wrap.k_hi and chan.k_lo
    bit-identical via numpy view aliasing on a single ``nodes_all``
    allocation per wrap segment.  The body wrap's k_hi is split into
    3 i-segments at user-tunable indices ``i_split_body_upper`` and
    ``i_split_body_lower`` (defaults: ``ni_body // 3`` and
    ``2 * ni_body // 3``), each sub-channel sharing only its own
    chunk of body.k_hi via numpy view aliasing.  Adjacent
    chan_body_* sub-channels share their TE-direction seam by
    overlapping on the split column (identity AxisMap).

    The 2 (or 3) outflow channels are NOT view-aliased: their k=0
    inner face has shape ``(nk_wrap, nj, 3)`` derived from the wrap's
    ``i_lo`` / ``i_hi`` faces (a cross-axis remap from the wrap
    convention), so :meth:`fill_channel_blocks` copies the face data
    explicitly each time.  The SHARED face-graph entries are still
    emitted so cc.par writers see the connection.

    **Default face BCs**:

    - Wrap k_lo = WALL on body, SHARED on upper/lower (wake-cut /
      wake-base).
    - Wrap k_hi = SHARED with the corresponding chan.k_lo.
    - Wrap i_lo / i_hi at the wake far-tips = SHARED with chan_outflow_*.
    - Channel k_lo = SHARED.  Channel k_hi = FREE (frame edge).
    - Channel j_lo / j_hi = FREE (wing tips, Phase 2 will optionally
      add spanwise tip channels).

    Channels do NOT march; :meth:`shore.mesh.Mesh.march` runs only
    the C-grid wrap marcher and then calls :meth:`fill_channel_blocks`.

    Parameters
    ----------
    freestream_world :
        Length-3 world-frame freestream vector.  ``None`` defaults to
        the +chord_axis (chord_axis = first non-span world axis).
    span_axis :
        ``"x"``, ``"y"``, or ``"z"``.  Default ``"z"``.
    chord_min, chord_max, normal_min, normal_max :
        Frame extents in world coordinates.  Must contain the marched
        wrap envelope (validated in :meth:`fill_channel_blocks`).
    nk_channel :
        Total channel-block layers (incl. shared-seam slot at k=0).
        Must be >= 2.
    k_law_channel :
        Spacing law along channel k: ``'uniform'`` / ``'tanh'`` /
        ``'tanh2'``.  Default ``'uniform'``.
    k_beta_channel :
        Clustering strength for tanh / tanh2.
    c1_seam :
        When True (default), the channel's first step at every (i, j)
        is pinned to the wrap's last marching step at the same node —
        bit-exact C1 spacing match across the wrap-channel seam.  Only
        applies to the 5 ``k_hi``-side channels (the outflow channels
        bridge wrap *side* faces and have no wrap "last step" along
        their own k-axis to inherit from).
    i_split_body_upper, i_split_body_lower :
        Body-block i-indices at which ``body.k_hi`` is split into
        the 3 chan_body_* sub-channels.  ``None`` (default) =
        ``ni_body // 3`` and ``2 * ni_body // 3``.
    """

    def __init__(
        self,
        *,
        chord_min: float,
        chord_max: float,
        normal_min: float,
        normal_max: float,
        nk_channel: int,
        freestream_world: np.ndarray | None = None,
        span_axis: str = "z",
        k_law_channel: str = "uniform",
        k_beta_channel: float = 3.0,
        c1_seam: bool = True,
        i_split_body_upper: int | None = None,
        i_split_body_lower: int | None = None,
        ni_outflow: int = 8,
    ) -> None:
        if nk_channel < 2:
            raise ParameterError(f"CHTopology: nk_channel must be >= 2, got {nk_channel}")
        if k_law_channel not in ("uniform", "tanh", "tanh2"):
            raise ParameterError(
                f"CHTopology: k_law_channel must be 'uniform' / 'tanh' / "
                f"'tanh2', got {k_law_channel!r}"
            )
        # resolve_frame validates span_axis, freestream_world, and the
        # chord/normal extents.
        self.frame = resolve_frame(
            freestream_world=freestream_world,
            span_axis=span_axis,
            chord_min=chord_min,
            chord_max=chord_max,
            normal_min=normal_min,
            normal_max=normal_max,
        )
        self.freestream_world = (
            None if freestream_world is None else np.asarray(freestream_world, dtype=np.float64)
        )
        self.span_axis = span_axis
        self.nk_channel = int(nk_channel)
        self.k_law_channel = k_law_channel
        self.k_beta_channel = float(k_beta_channel)
        self.c1_seam = bool(c1_seam)
        self.i_split_body_upper = i_split_body_upper
        self.i_split_body_lower = i_split_body_lower
        if ni_outflow < 2:
            raise ParameterError(f"CHTopology: ni_outflow must be >= 2, got {ni_outflow}")
        self.ni_outflow = int(ni_outflow)
        # Set in build() — fill_channel_blocks needs them later.
        self._is_blunt: bool = False
        self._ni_wake: int = 0
        self._ni_body: int = 0
        self._nk_wrap: int = 0

    @property
    def is_blunt(self) -> bool:
        """True when build() detected a blunt TE (CH-4)."""
        return self._is_blunt

    def build(
        self,
        surface: np.ndarray | None,
        nk: int,
        *,
        ni_body: int | None = None,
        ni_wake: int | None = None,
        n_stations: int | None = None,
        wake_length: float | None = None,
        wake_growth: float = 1.15,
        te_dihedral_deg: float = 30.0,
        spacing_k: Spacing1D | None = None,
        mesh: object | None = None,
        nk_w: int = 3,
        **kwargs,
    ) -> list[HexBlock]:
        """Build the CH topology block list.

        ``nk`` is the wrap depth.  The channel depth is taken from
        ``CHTopology.nk_channel`` set at construction.
        """
        # --------------------------------------------------------------
        # 1. Build / accept the C-grid k=0 layer.  Detect TE bluntness.
        # --------------------------------------------------------------
        if surface is None:
            if mesh is None:
                raise ParameterError(
                    "CHTopology.build needs either a pre-built `surface` or a "
                    "trimesh.Trimesh `mesh=` argument."
                )
            for name, val in (
                ("ni_body", ni_body),
                ("ni_wake", ni_wake),
                ("n_stations", n_stations),
                ("wake_length", wake_length),
            ):
                if val is None:
                    raise ParameterError(f"CHTopology.build needs `{name}=` when surface is None.")
            k0 = build_cgrid_k0(
                mesh,  # type: ignore[arg-type]
                ni_body=int(ni_body),  # type: ignore[arg-type]
                ni_wake=int(ni_wake),  # type: ignore[arg-type]
                n_stations=int(n_stations),  # type: ignore[arg-type]
                wake_length=float(wake_length),  # type: ignore[arg-type]
                span_axis=self.span_axis,
                freestream_world=self.freestream_world,
                te_dihedral_deg=te_dihedral_deg,
                wake_growth=wake_growth,
            )
            surface_arr = k0.nodes
            _ni_body = k0.ni_body
            _ni_wake = k0.ni_wake
            _is_blunt = k0.te_result.is_blunt
        else:
            if ni_body is None or ni_wake is None:
                raise ParameterError(
                    "CHTopology.build with a pre-built `surface` requires explicit "
                    "`ni_body=` and `ni_wake=`."
                )
            _ni_body = int(ni_body)
            _ni_wake = int(ni_wake)
            _is_blunt = False
            surface_arr = surface
            if surface_arr.shape[0] != 2 * (_ni_wake - 1) + _ni_body:
                raise ParameterError(
                    f"CHTopology.build: surface i-extent {surface_arr.shape[0]} "
                    f"does not match 2*(ni_wake-1) + ni_body = "
                    f"{2 * (_ni_wake - 1) + _ni_body}."
                )

        self._is_blunt = _is_blunt
        self._ni_wake = _ni_wake
        self._ni_body = _ni_body
        self._nk_wrap = int(nk)

        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

        # Resolve body-split indices.
        if self.i_split_body_upper is None:
            i_split_upper = _ni_body // 3
        else:
            i_split_upper = int(self.i_split_body_upper)
        if self.i_split_body_lower is None:
            i_split_lower = (2 * _ni_body) // 3
        else:
            i_split_lower = int(self.i_split_body_lower)
        if not (0 < i_split_upper < i_split_lower < _ni_body - 1):
            raise ParameterError(
                f"CHTopology.build: invalid body splits "
                f"(0 < i_split_body_upper={i_split_upper} < "
                f"i_split_body_lower={i_split_lower} < ni_body-1={_ni_body - 1}) "
                f"violated.  Adjust the split indices or increase ni_body."
            )
        self._i_split_upper = i_split_upper
        self._i_split_lower = i_split_lower

        # --------------------------------------------------------------
        # 2. Allocate the wrap+channel buffer with depth total_k.
        # --------------------------------------------------------------
        # Channel depth includes the SHARED seam slot at k=0 (= wrap's
        # k=nk_wrap-1).  Total k-buffer length = nk_wrap + nk_channel - 1.
        nk_wrap = int(nk)
        nk_chan = self.nk_channel
        total_k = nk_wrap + nk_chan - 1

        nj = surface_arr.shape[1]
        ni_total = 2 * (_ni_wake - 1) + _ni_body
        nodes_all = np.zeros((ni_total, nj, total_k, 3), dtype=np.float64)
        nodes_all[:, :, 0, :] = surface_arr

        i_te_upper = _ni_wake - 1
        i_te_lower = i_te_upper + _ni_body - 1

        # --- Wrap views (depth nk_wrap, k=0..nk_wrap-1) ---------------
        upper_wrap = nodes_all[: i_te_upper + 1, :, :nk_wrap, :]
        body_wrap = nodes_all[i_te_upper : i_te_lower + 1, :, :nk_wrap, :]
        lower_wrap = nodes_all[i_te_lower:, :, :nk_wrap, :]

        upper = HexBlock(
            label="upper",
            nodes=upper_wrap,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        body = HexBlock(
            label="body",
            nodes=body_wrap,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        lower = HexBlock(
            label="lower",
            nodes=lower_wrap,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )

        # --- Wrap SHARED seams (mirrors CGridTopology) ----------------
        upper._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=body.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        body._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=upper.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        body._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=body.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )

        wrap_blocks: list[HexBlock] = [upper, body, lower]

        # --- Wake-base block (CH-4 only) ------------------------------
        if _is_blunt:
            if nk_w < 2:
                raise ParameterError(f"CHTopology.build: nk_w must be >= 2, got {nk_w}")
            upper_k0 = upper_wrap[:, :, 0, :]
            lower_k0 = lower_wrap[:, :, 0, :]
            lower_k0_rev = lower_k0[::-1, :, :]
            ts = np.linspace(0.0, 1.0, nk_w)
            wake_base_nodes = np.zeros((_ni_wake, nj, nk_w, 3), dtype=np.float64)
            for k_idx, t in enumerate(ts):
                wake_base_nodes[:, :, k_idx, :] = (1.0 - t) * lower_k0_rev + t * upper_k0
            wake_base = HexBlock(
                label="wake_base",
                nodes=wake_base_nodes,
                periodic_i=False,
                periodic_j=False,
                spacing_k=_sp_k,
            )
            wake_base._set_face_bc(
                "k_hi",
                FaceBC.SHARED,
                partner=upper.face("k_lo"),
                axis_map=AxisMap.identity("k_hi"),
            )
            upper._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=wake_base.face("k_hi"),
                axis_map=AxisMap.identity("k_lo"),
            )
            i_flip = AxisMap.from_pairs(("i", True), ("j", False))
            wake_base._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=lower.face("k_lo"),
                axis_map=i_flip,
            )
            lower._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=wake_base.face("k_lo"),
                axis_map=i_flip.inverse("k_lo", "k_lo"),
            )
            wake_base._set_face_bc("i_hi", FaceBC.WALL)
            wrap_blocks.append(wake_base)
        else:
            # Sharp-TE wake cut: upper.k_lo ↔ lower.k_lo direct seam.
            wake_cut_axis_map = AxisMap.from_pairs(("i", True), ("j", False))
            upper._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=lower.face("k_lo"),
                axis_map=wake_cut_axis_map,
            )
            lower._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=upper.face("k_lo"),
                axis_map=wake_cut_axis_map.inverse("k_lo", "k_lo"),
            )

        # --------------------------------------------------------------
        # 3. Carve the 5 k_hi-side channel views from nodes_all.
        # --------------------------------------------------------------
        # Each channel's k-axis runs from the wrap's k=nk_wrap-1 (the
        # shared seam) up through the deeper buffer slots.
        chan_upper_view = nodes_all[: i_te_upper + 1, :, nk_wrap - 1 :, :]
        # Body-block slice in nodes_all: [i_te_upper, i_te_lower].  Split
        # at i_split_upper (relative to body) -> nodes_all index
        # i_te_upper + i_split_upper.  Channel views overlap by 1 row at
        # the split columns to make those seams SHARED-by-view.
        i_body_split_upper_global = i_te_upper + i_split_upper
        i_body_split_lower_global = i_te_upper + i_split_lower
        chan_body_top_view = nodes_all[
            i_te_upper : i_body_split_upper_global + 1, :, nk_wrap - 1 :, :
        ]
        chan_body_front_view = nodes_all[
            i_body_split_upper_global : i_body_split_lower_global + 1, :, nk_wrap - 1 :, :
        ]
        chan_body_bot_view = nodes_all[
            i_body_split_lower_global : i_te_lower + 1, :, nk_wrap - 1 :, :
        ]
        chan_lower_view = nodes_all[i_te_lower:, :, nk_wrap - 1 :, :]

        chan_upper = HexBlock(
            label="chan_upper",
            nodes=chan_upper_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        chan_body_top = HexBlock(
            label="chan_body_top",
            nodes=chan_body_top_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        chan_body_front = HexBlock(
            label="chan_body_front",
            nodes=chan_body_front_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        chan_body_bot = HexBlock(
            label="chan_body_bot",
            nodes=chan_body_bot_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        chan_lower = HexBlock(
            label="chan_lower",
            nodes=chan_lower_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )

        khi_channels = [chan_upper, chan_body_top, chan_body_front, chan_body_bot, chan_lower]

        # --- Wire wrap.k_hi ↔ chan.k_lo (SHARED, view-aliased) --------
        # Each k_hi-side channel partners with a slice of one wrap block.
        # The wrap-side face is the partial wrap.k_hi covering the
        # channel's i-extent; we still emit the SHARED face on the
        # whole wrap block (the face graph treats the whole face as one
        # connection) — for the body wrap, that means body.k_hi has 3
        # SHARED partners.  Face wiring only stores ONE partner per
        # face on a HexBlock, so we use the body's i-range to identify
        # the dominant partner; downstream consumers (cc.par writer)
        # need the .adjacency.json sidecar to reconstruct the full
        # body.k_hi → 3 channels mapping.  For Phase 1 we wire each
        # channel.k_lo back to the relevant wrap's k_hi face and let
        # the wrap's k_hi face partner with the central body-front
        # channel (the one whose chord interval covers the LE — the
        # most representative).  This is a documented Phase-1
        # limitation; Phase 2 will add multi-partner face support.
        upper._set_face_bc(
            "k_hi",
            FaceBC.SHARED,
            partner=chan_upper.face("k_lo"),
            axis_map=AxisMap.identity("k_hi"),
        )
        chan_upper._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=upper.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )
        body._set_face_bc(
            "k_hi",
            FaceBC.SHARED,
            partner=chan_body_front.face("k_lo"),
            axis_map=AxisMap.identity("k_hi"),
        )
        chan_body_top._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=body.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )
        chan_body_front._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=body.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )
        chan_body_bot._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=body.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )
        lower._set_face_bc(
            "k_hi",
            FaceBC.SHARED,
            partner=chan_lower.face("k_lo"),
            axis_map=AxisMap.identity("k_hi"),
        )
        chan_lower._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=lower.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )

        # --- Channel-channel lateral SHARED seams (i↔i identity) ------
        # chan_upper.i_hi  ↔ chan_body_top.i_lo  (upper-TE column)
        chan_upper._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=chan_body_top.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        chan_body_top._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_upper.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        # chan_body_top.i_hi ↔ chan_body_front.i_lo (upper body split)
        chan_body_top._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=chan_body_front.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        chan_body_front._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_body_top.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        # chan_body_front.i_hi ↔ chan_body_bot.i_lo (lower body split)
        chan_body_front._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=chan_body_bot.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        chan_body_bot._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_body_front.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        # chan_body_bot.i_hi ↔ chan_lower.i_lo (lower-TE column)
        chan_body_bot._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=chan_lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        chan_lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_body_bot.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )

        # --------------------------------------------------------------
        # 4. Build the 4 (or 5) outflow channels (Phase-1.2 redesign).
        # --------------------------------------------------------------
        # Each side of the outflow region is partitioned into TWO
        # sub-channels along the normal direction:
        #
        #   chan_outflow_top_upper:  normal in [normal_env_upper, normal_max]
        #                            conformal SHARED to chan_upper.i_lo
        #                            (channel.k count = nk_channel)
        #
        #   chan_outflow_wake_upper: normal in [+te/2, normal_env_upper]
        #                            conformal SHARED to upper.i_lo
        #                            (channel.k count = nk_wrap)
        #
        # Mirror layout on the lower side.  For blunt TE, a wake_base
        # outflow channel sits between wake_lower and wake_upper at
        # normal in [-te/2, +te/2] (channel.k count = nk_w).  All blocks
        # share their lateral seams via numpy view aliasing on a single
        # ``outflow_buf`` allocation per side.  Channel axes for every
        # outflow block:
        #
        #     i (channel)  -> chord  (wrap-far-tip -> rear frame)
        #     j (channel)  -> span
        #     k (channel)  -> normal (low normal -> high normal)
        #
        # Conformal SHARED seams to the wrap and to the k_hi-side
        # channels carry identity AxisMaps (same-axis i_lo <-> i_lo
        # or i_lo <-> i_hi depending on the partner's i orientation).
        ni_outflow = self.ni_outflow

        # Convention: every outflow block has channel-k pointing OUTWARD
        # from the wake-cut.  Upper-side blocks: k=0 at wake-cut (or at
        # +te/2 for blunt), k=last at top frame.  Lower-side blocks:
        # k=0 at wake-cut (or at -te/2 for blunt), k=last at bottom frame.
        # This matches the wrap blocks' k orientation (wrap-upper.k=0 on
        # wake-cut, wrap-upper.k=last on the marched envelope) so the
        # conformal SHARED seams to upper.i_lo / lower.i_hi can use
        # identity AxisMaps without k-flips.
        #
        # Buffer layout per side: stack [wake_half | top_half] along
        # channel-k.  wake_half has depth nk_wrap; top_half has depth
        # nk_chan; they share their seam row at the marched envelope
        # normal (1-row overlap).
        #
        #   upper_buf shape (ni_outflow, nj, nk_wrap + nk_chan - 1, 3)
        #     wake_upper view: [..., :nk_wrap, :]
        #     top_upper  view: [..., nk_wrap - 1 :, :]
        #
        # Lower mirrors that.  Wake-cut SHARED between wake_upper.k_lo
        # and wake_lower.k_lo (sharp) is geometrically coincident but
        # lives in different buffers — the SHARED face graph still
        # captures the relation; the geometry is enforced by writing
        # bit-identical wake-cut data to both in fill_channel_blocks.
        nk_outflow_side = nk_wrap + nk_chan - 1
        upper_buf = np.zeros((ni_outflow, nj, nk_outflow_side, 3), dtype=np.float64)
        lower_buf = np.zeros((ni_outflow, nj, nk_outflow_side, 3), dtype=np.float64)

        wake_up_view = upper_buf[:, :, :nk_wrap, :]
        top_up_view = upper_buf[:, :, nk_wrap - 1 :, :]
        wake_lo_view = lower_buf[:, :, :nk_wrap, :]
        top_lo_view = lower_buf[:, :, nk_wrap - 1 :, :]

        chan_outflow_wake_upper = HexBlock(
            label="chan_outflow_wake_upper",
            nodes=wake_up_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        chan_outflow_top_upper = HexBlock(
            label="chan_outflow_top_upper",
            nodes=top_up_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        chan_outflow_wake_lower = HexBlock(
            label="chan_outflow_wake_lower",
            nodes=wake_lo_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        chan_outflow_top_lower = HexBlock(
            label="chan_outflow_top_lower",
            nodes=top_lo_view,
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )

        # ----- Internal envelope-seam SHARED (view-aliased) ------------------
        # wake_upper.k_hi <-> top_upper.k_lo at envelope-top normal.
        # wake_lower.k_hi <-> top_lower.k_lo at envelope-bot normal.
        chan_outflow_wake_upper._set_face_bc(
            "k_hi",
            FaceBC.SHARED,
            partner=chan_outflow_top_upper.face("k_lo"),
            axis_map=AxisMap.identity("k_hi"),
        )
        chan_outflow_top_upper._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=chan_outflow_wake_upper.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )
        chan_outflow_wake_lower._set_face_bc(
            "k_hi",
            FaceBC.SHARED,
            partner=chan_outflow_top_lower.face("k_lo"),
            axis_map=AxisMap.identity("k_hi"),
        )
        chan_outflow_top_lower._set_face_bc(
            "k_lo",
            FaceBC.SHARED,
            partner=chan_outflow_wake_lower.face("k_hi"),
            axis_map=AxisMap.identity("k_lo"),
        )

        # ----- Wake-cut SHARED -----------------------------------------------
        # Sharp: wake_upper.k_lo (at normal=0) <-> wake_lower.k_lo (at normal=0).
        #   Same-axis k<->k.  Both blocks have channel-k pointing outward
        #   (away from the wake-cut), so their k=0 faces are coincident
        #   but their k axes point in OPPOSITE directions.  AxisMap on
        #   wake_upper.k_lo: along-axes (i, j) of wake_upper map to
        #   along-axes (i, j) of wake_lower.  Channel-i for both runs
        #   chord_far -> chord_max in the same direction; channel-j same.
        #   Identity on (i, j).
        # Blunt: wake_upper.k_lo (at +te/2) <-> wake_base.k_hi.
        #        wake_lower.k_lo (at -te/2) <-> wake_base.k_lo.
        if _is_blunt:
            chan_outflow_wb_nodes = np.zeros((ni_outflow, nj, nk_w, 3), dtype=np.float64)
            chan_outflow_wb = HexBlock(
                label="chan_outflow_wake_base",
                nodes=chan_outflow_wb_nodes,
                periodic_i=False,
                periodic_j=False,
                spacing_k=_sp_k,
            )
            # wake_lower.k_lo (-te/2 from below) <-> wake_base.k_lo (-te/2 from above)
            chan_outflow_wake_lower._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=chan_outflow_wb.face("k_lo"),
                axis_map=AxisMap.from_pairs(
                    ("i", False), ("j", True)
                ),  # wake_base's j reversed (ADR 0057)
            )
            chan_outflow_wb._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=chan_outflow_wake_lower.face("k_lo"),
                axis_map=AxisMap.from_pairs(("i", False), ("j", True)),
            )
            # wake_upper.k_lo (+te/2 from above) <-> wake_base.k_hi
            chan_outflow_wake_upper._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=chan_outflow_wb.face("k_hi"),
                axis_map=AxisMap.identity("k_lo"),
            )
            chan_outflow_wb._set_face_bc(
                "k_hi",
                FaceBC.SHARED,
                partner=chan_outflow_wake_upper.face("k_lo"),
                axis_map=AxisMap.identity("k_hi"),
            )
        else:  # j reversed on the upper side (ADR 0057)
            chan_outflow_wake_upper._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=chan_outflow_wake_lower.face("k_lo"),
                axis_map=AxisMap.from_pairs(("i", False), ("j", True)),
            )
            chan_outflow_wake_lower._set_face_bc(
                "k_lo",
                FaceBC.SHARED,
                partner=chan_outflow_wake_upper.face("k_lo"),
                axis_map=AxisMap.from_pairs(("i", False), ("j", True)),
            )

        # ----- Conformal SHARED seams to the wrap and k_hi-side channels -----
        # Top-half seams to the k_hi-side channels (chan_upper / chan_lower).
        # chan_upper.i_lo is at chord=chord_upper_far with normal range
        # [normal_env_upper, normal_max]; channel.k along normal — matches
        # chan_outflow_top_upper.i_lo exactly (nj * nk_chan).  Identity
        # AxisMap, same-axis i_lo <-> i_lo.
        chan_upper._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_outflow_top_upper.face("i_lo"),
            axis_map=AxisMap.from_pairs(("j", True), ("k", False)),  # j reversed (ADR 0057)
        )
        chan_outflow_top_upper._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_upper.face("i_lo"),
            axis_map=AxisMap.from_pairs(("j", True), ("k", False)),
        )
        # chan_lower.i_hi is at chord=chord_lower_far (since chan_lower's i
        # runs lower-TE -> lower-far-tip).  Pair to chan_outflow_top_lower.i_lo.
        chan_lower._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=chan_outflow_top_lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        chan_outflow_top_lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_lower.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )

        # Wake-half seams to the wrap upper / lower blocks.
        # upper.i_lo is at chord=chord_upper_far with normal range [0,
        # normal_env_upper] (or [+te/2, normal_env_upper] for blunt);
        # nk = nk_wrap; same as chan_outflow_wake_upper.i_lo.  Identity
        # AxisMap, same-axis i_lo <-> i_lo.
        upper._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=chan_outflow_wake_upper.face("i_lo"),
            axis_map=AxisMap.from_pairs(("j", True), ("k", False)),  # j reversed (ADR 0057)
        )
        chan_outflow_wake_upper._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=upper.face("i_lo"),
            axis_map=AxisMap.from_pairs(("j", True), ("k", False)),
        )
        # lower.i_hi is at chord=chord_lower_far (lower's i runs TE ->
        # far-tip).  Pair to chan_outflow_wake_lower.i_lo (i_hi <-> i_lo
        # identity AxisMap).
        lower._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=chan_outflow_wake_lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        chan_outflow_wake_lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=lower.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )

        # Blunt: wake_base.i_lo <-> chan_outflow_wb.i_lo.
        if _is_blunt:
            wake_base._set_face_bc(
                "i_lo",
                FaceBC.SHARED,
                partner=chan_outflow_wb.face("i_lo"),
                axis_map=AxisMap.from_pairs(("j", True), ("k", False)),  # j reversed (ADR 0057)
            )
            chan_outflow_wb._set_face_bc(
                "i_lo",
                FaceBC.SHARED,
                partner=wake_base.face("i_lo"),
                axis_map=AxisMap.from_pairs(("j", True), ("k", False)),
            )

        # Outermost outflow faces are at the bottom / top frame edges.
        # With the convention "channel-k points outward from wake-cut":
        #   top_upper.k_hi = top frame    (HexBlock default FREE — OK)
        #   top_lower.k_hi = bottom frame (default is FREE — OK; lower
        #     stack's k_hi face is the outer-most slot, not k_lo).
        # No override needed.

        outflow_channels = [
            chan_outflow_top_upper,
            chan_outflow_wake_upper,
            chan_outflow_wake_lower,
            chan_outflow_top_lower,
        ]
        if _is_blunt:
            outflow_channels.append(chan_outflow_wb)

        return [*wrap_blocks, *khi_channels, *outflow_channels]

    def fill_channel_blocks(self, blocks: list[HexBlock]) -> None:
        """Populate every channel block's interior after wrap marching.

        Walks the channel blocks in order, computes each channel's
        outer quad on the user-specified frame, and fills the
        ``k=1..nk_channel-1`` interior by ruled-surface frustum
        interpolation.  C1 seam pinning (per-(i, j) first step =
        wrap's last marching step) is applied to the 5 k_hi-side
        channels when ``self.c1_seam`` is True.
        """
        # Identify wrap blocks and channel blocks by label.  Robust to
        # caller-injected ordering changes.
        labels = [b.label for b in blocks]
        upper = blocks[labels.index("upper")]
        body = blocks[labels.index("body")]
        lower = blocks[labels.index("lower")]
        wake_base = blocks[labels.index("wake_base")] if self._is_blunt else None

        chan_upper = blocks[labels.index("chan_upper")]
        chan_body_top = blocks[labels.index("chan_body_top")]
        chan_body_front = blocks[labels.index("chan_body_front")]
        chan_body_bot = blocks[labels.index("chan_body_bot")]
        chan_lower = blocks[labels.index("chan_lower")]
        chan_outflow_top_upper = blocks[labels.index("chan_outflow_top_upper")]
        chan_outflow_wake_upper = blocks[labels.index("chan_outflow_wake_upper")]
        chan_outflow_wake_lower = blocks[labels.index("chan_outflow_wake_lower")]
        chan_outflow_top_lower = blocks[labels.index("chan_outflow_top_lower")]
        chan_outflow_wb = blocks[labels.index("chan_outflow_wake_base")] if self._is_blunt else None

        # Span extents from the wrap envelope's j_lo / j_hi rows.
        sp_idx = self.frame.span_axis_idx
        span_lo = float(body.nodes[:, 0, 0, sp_idx].mean())
        span_hi = float(body.nodes[:, -1, 0, sp_idx].mean())

        # k_hi envelopes at the marched outer layer.
        upper_env = upper.nodes[:, :, -1, :]  # (ni_wake, nj, 3)
        body_env = body.nodes[:, :, -1, :]
        lower_env = lower.nodes[:, :, -1, :]
        wake_base_env = wake_base.nodes if self._is_blunt and wake_base is not None else None

        # Validate that the marched envelope lies inside the frame
        # (chord, normal axes only).  Tolerance: the channel can
        # interpolate inward as long as the envelope is on the wrap
        # side of the frame edges; reject when envelope crosses out.
        cax = self.frame.chord_axis
        nax = self.frame.normal_axis
        envelope_pts = np.concatenate(
            [upper_env.reshape(-1, 3), body_env.reshape(-1, 3), lower_env.reshape(-1, 3)],
            axis=0,
        )
        env_chord = envelope_pts[:, cax]
        env_normal = envelope_pts[:, nax]
        if not (
            self.frame.chord_min <= env_chord.min()
            and env_chord.max() <= self.frame.chord_max
            and self.frame.normal_min <= env_normal.min()
            and env_normal.max() <= self.frame.normal_max
        ):
            raise ParameterError(
                f"CHTopology.fill_channel_blocks: marched wrap envelope is "
                f"outside the frame [chord {self.frame.chord_min},"
                f"{self.frame.chord_max}] x [normal {self.frame.normal_min},"
                f"{self.frame.normal_max}].  Envelope chord range "
                f"[{env_chord.min()}, {env_chord.max()}], normal range "
                f"[{env_normal.min()}, {env_normal.max()}].  Either "
                f"enlarge the frame or reduce nk / wrap step sizes."
            )

        # Compute outer quads for all 7 (or 8) channels.
        quads = outer_quads_for_cgrid(
            frame=self.frame,
            upper_envelope=upper_env,
            body_envelope=body_env,
            lower_envelope=lower_env,
            i_split_body_upper=self._i_split_upper,
            i_split_body_lower=self._i_split_lower,
            span_lo=span_lo,
            span_hi=span_hi,
            is_blunt=self._is_blunt,
            wake_base_envelope=wake_base_env,
        )

        # ---- 5 k_hi-side channels: bilinear frustum + optional C1 ----
        khi_specs: list[tuple[HexBlock, np.ndarray, str]] = [
            (chan_upper, upper.nodes, "upper"),
            (chan_body_top, body.nodes, "body_top"),
            (chan_body_front, body.nodes, "body_front"),
            (chan_body_bot, body.nodes, "body_bot"),
            (chan_lower, lower.nodes, "lower"),
        ]

        # The 5 k_hi-side channels share their k=0 with a wrap.k_hi
        # slice.  C1 first-step is the wrap's last layer thickness
        # at every (i, j) within the channel's i-extent.
        # Sub-block i-ranges within the wrap's nodes for the 3 body
        # channels:
        body_i_ranges = {
            "body_top": (0, self._i_split_upper + 1),
            "body_front": (self._i_split_upper, self._i_split_lower + 1),
            "body_bot": (self._i_split_lower, body.nodes.shape[0]),
        }

        for chan, wrap_nodes, key in khi_specs:
            inner = chan.nodes[:, :, 0, :]  # SHARED-by-view with wrap.k_hi slice
            outer_quad = quads[key]
            if self.c1_seam:
                # last-step on the wrap, per-(i, j), within the
                # channel's i-extent.
                if key in body_i_ranges:
                    i_lo, i_hi = body_i_ranges[key]
                    last_step = np.linalg.norm(
                        wrap_nodes[i_lo:i_hi, :, -1, :] - wrap_nodes[i_lo:i_hi, :, -2, :],
                        axis=-1,
                    )
                else:
                    last_step = np.linalg.norm(
                        wrap_nodes[:, :, -1, :] - wrap_nodes[:, :, -2, :], axis=-1
                    )
                c1 = last_step
            else:
                c1 = None
            frustum = build_channel_block(
                inner_face=inner,
                outer_quad=outer_quad,
                nk_channel=self.nk_channel,
                k_law=self.k_law_channel,
                k_beta=self.k_beta_channel,
                c1_first_step=c1,
            )
            chan.nodes[:, :, 1:, :] = frustum[:, :, 1:, :]

        # ---- 4 (or 5) outflow channels: conformal SHARED to wrap ----
        # The outflow region downstream of the wrap is partitioned per
        # side into:
        #   - top_X:   conformal SHARED to chan_X.i_lo/i_hi
        #              (channel-k count = nk_chan; normal range
        #               [envelope_X, frame_top/bot])
        #   - wake_X:  conformal SHARED to wrap.X.i_lo/i_hi
        #              (channel-k count = nk_wrap; normal range
        #               [wake_cut, envelope_X])
        # plus, for blunt TE, a wake_base outflow bridging the TE gap.
        #
        # Filling strategy: each outflow block's i_lo face (at chord =
        # chord_X_far) is COPIED bit-exact from the partner wrap or
        # k_hi-side channel face, so the conformal SHARED is honoured.
        # The i_hi face (at chord = chord_max) inherits the same
        # (j, normal) pattern but with chord = chord_max.  Interior is
        # bilinear in chord between the two faces, with the chord
        # distribution shaped by k_law_channel.
        from shore.volume.spacing import build_i_coords

        cax = self.frame.chord_axis
        rear = self.frame.chord_max

        def _fill_outflow_block(
            chan: HexBlock,
            *,
            inner_face: np.ndarray,
        ) -> None:
            """Fill a rectangular outflow block from a copied inner face.

            ``inner_face`` is the (nj, nk_face, 3) data to plant on the
            channel's i_lo face (chord = chord_far).  The outer i_hi
            face has the same span/normal coordinates but chord =
            chord_max.  Interior columns are bilinear in chord under
            ``k_law_channel`` spacing.
            """
            ni_c = chan.nodes.shape[0]
            nj_c, nk_c, _ = inner_face.shape
            assert chan.nodes.shape[1:3] == (nj_c, nk_c), (
                f"outflow block {chan.label}: inner face shape mismatch "
                f"{inner_face.shape} vs {chan.nodes.shape}"
            )

            # i_lo face = inner_face copy.
            chan.nodes[0, :, :, :] = inner_face

            # i_hi face = inner with chord swapped to rear.
            outer_face = inner_face.copy()
            outer_face[..., cax] = rear
            chan.nodes[-1, :, :, :] = outer_face

            # Interior: blend between inner and outer along channel-i
            # under k_law_channel.
            t = build_i_coords(
                ni_c, 0.0, 1.0, spacing=self.k_law_channel, beta=self.k_beta_channel
            )  # (ni_c,)
            t_arr = t[1:-1, None, None, None]  # (ni_c-2, 1, 1, 1)
            inner4 = inner_face[None, ...]
            outer4 = outer_face[None, ...]
            chan.nodes[1:-1, :, :, :] = (1.0 - t_arr) * inner4 + t_arr * outer4

        # ---- Inner-face data sources (conformal SHARED partners) ----
        # top_upper inner face = chan_upper.i_lo face (= chan_upper.nodes[0, :, :, :])
        # top_lower inner face = chan_lower.i_hi face (= chan_lower.nodes[-1, :, :, :])
        # wake_upper inner face = upper.i_lo face (= upper.nodes[0, :, :, :])
        # wake_lower inner face = lower.i_hi face (= lower.nodes[-1, :, :, :])
        # wake_base inner face (blunt) = wake_base.i_lo face
        # The upper-side blocks and wake_base have channel-k pointing up
        # (+normal): with i along the chord and j along the span they would
        # be left-handed, so their j runs against the partner's (ADR 0057).
        _fill_outflow_block(
            chan_outflow_top_upper,
            inner_face=chan_upper.nodes[0, ::-1, :, :].copy(),
        )
        _fill_outflow_block(
            chan_outflow_top_lower,
            inner_face=chan_lower.nodes[-1, :, :, :].copy(),
        )
        _fill_outflow_block(
            chan_outflow_wake_upper,
            inner_face=upper.nodes[0, ::-1, :, :].copy(),
        )
        _fill_outflow_block(
            chan_outflow_wake_lower,
            inner_face=lower.nodes[-1, :, :, :].copy(),
        )
        if self._is_blunt and chan_outflow_wb is not None and wake_base is not None:
            _fill_outflow_block(
                chan_outflow_wb,
                inner_face=wake_base.nodes[0, ::-1, :, :].copy(),
            )

        # The conformal-SHARED partners are i_lo/i_hi face-graph entries.
        # We just wrote the channel's i_lo face from a COPY of the wrap
        # face data; the wrap face itself is unchanged, so geometric
        # coincidence is exact (bit-identical) at the seam.

        _ = quads  # outflow quads are now documentation-only; computed for
        # side-effect validation in outer_quads_for_cgrid.


def _projector_from_mesh(mesh: object) -> object:
    """Build a temporary Projector using mesh.centroid (legacy mesh= path)."""
    from shore.surface.projector import Projector

    return Projector(mesh, mesh.centroid)  # type: ignore[attr-defined]
