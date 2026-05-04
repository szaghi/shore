"""Topology ABC and concrete topologies (cubed-sphere, O-grid, C-grid, OH-hybrid)."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from shore.exceptions import ParameterError
from shore.mesh.block import HexBlock
from shore.mesh.face import AxisMap, FaceBC, FaceName
from shore.mesh.spacing import Spacing1D
from shore.surface.cap import build_cap_k0
from shore.surface.cgrid import build_cgrid_k0
from shore.surface.equator import build_equator
from shore.topology.oh_channel import build_channel_block
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
        theta_cap=theta_cap, nc=nc, pole="south", projector=cap_proj,
        j_spacing=j_spacing, j_beta=j_beta, phi_offset=phi_offset,
    )
    cap_north_k0 = build_cap_k0(
        theta_cap=theta_cap, nc=nc, pole="north", projector=cap_proj,
        j_spacing=j_spacing, j_beta=j_beta, phi_offset=phi_offset,
    )

    # Equator k=0 — great-circle slerp + ray-cast.
    equator_k0 = build_equator(
        cap_south=cap_south_k0, cap_north=cap_north_k0, ni=ni,
        theta_cap=theta_cap, projector=cap_proj,
        i_spacing=i_spacing, i_beta=i_beta,
    )

    # Equatorial ring buffer of depth total_k.  Sub-block sees [..., :nk, :].
    ring_nj = 4 * (nc - 1) + 1
    ring = np.empty((ni, ring_nj, total_k, 3), dtype=np.float64)
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
    cap_north_full = np.empty((nc, nc, total_k, 3), dtype=np.float64)
    cap_south_full = np.empty((nc, nc, total_k, 3), dtype=np.float64)
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
                label=f"sub{q}", nodes=sub_wrap_views[q],
                periodic_i=False, periodic_j=False, spacing_k=_sp_k,
            )
        )
    cap_north_blk = HexBlock(
        label="cap_north", nodes=cap_north_wrap,
        periodic_i=False, periodic_j=False, spacing_k=_sp_k,
    )
    cap_south_blk = HexBlock(
        label="cap_south", nodes=cap_south_wrap,
        periodic_i=False, periodic_j=False, spacing_k=_sp_k,
    )

    # Sub-block adjacency: sub_q.j_hi ↔ sub_{q+1 mod 4}.j_lo (SHARED).
    for q in range(4):
        q_next = (q + 1) % 4
        sub_blocks[q]._set_face_bc(
            "j_hi", FaceBC.SHARED, partner=sub_blocks[q_next].face("j_lo"),
            axis_map=AxisMap.identity("j_hi"),
        )
        sub_blocks[q_next]._set_face_bc(
            "j_lo", FaceBC.SHARED, partner=sub_blocks[q].face("j_hi"),
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
                cap_face, FaceBC.DIRICHLET,
                partner=sub_blk.face(sub_face),
                axis_map=sub_axis_map.inverse(sub_face, cap_face),
            )
            sub_blk._set_face_bc(
                sub_face, FaceBC.DIRICHLET,
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
        theta_cap: float = 0.05 * np.pi,
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
        cap_te_frac: float = 0.95,
        cap_tip_inset: float | None = None,
        cap_n_cap_j: int = 12,
        cap_j_spacing: str = "tanh",
        cap_j_beta: float = 4.0,
        cap_nk: int | None = None,
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
            ``"none"`` (default), ``"jlo"``, ``"jhi"``, or ``"both"``.
            When set to anything other than ``"none"``, the C-grid
            mesh is augmented with **conformal** wing-tip cap blocks
            at the requested spanwise extreme(s).  The cap blocks
            close the wing-tip airfoil cross-section as a wall, with
            their airfoil-wall nodes coincident bit-exactly with the
            body's k=0 wall at the spanwise tip.

            For sharp-TE bodies the cap-tips path **auto-rounds the
            body's TE over its full span** by ``cap_te_thickness``
            (or ``cap_te_thickness_frac * chord`` if not given), so
            the body wrap and the cap meshes match without any
            chimera-fringe handover.  ``blocks=3`` is automatically
            promoted to ``blocks=4`` in that case (a stderr note is
            emitted).  ``blocks=1`` (legacy) + ``cap_tips != "none"``
            is rejected.

            For already-blunt bodies the existing TE thickness is used
            and no rounding is applied.
        cap_te_thickness :
            Absolute TE thickness for body-rounding when ``cap_tips !=
            "none"`` and the body is sharp.  Mutually exclusive with
            ``cap_te_thickness_frac``.
        cap_te_thickness_frac :
            Same, as a fraction of chord.  Default ``0.002`` (0.2%).
        cap_le_frac, cap_te_frac :
            Chord fractions for the cap's LE-cut and TE-cut.  Default
            ``0.05`` and ``0.95``.
        cap_tip_inset :
            Cap spanwise wall thickness.  Default = ``1e-3 * chord``.
        cap_n_cap_j :
            Cap j-axis (airfoil-thickness) node count.  Default ``12``.
        cap_j_spacing, cap_j_beta :
            Cap j-axis spacing law (``"uniform"`` / ``"tanh"`` /
            ``"tanh2"``) and clustering strength.  Default ``"tanh"``,
            beta ``4.0``.
        cap_nk :
            Cap k-axis (spanwise) node count.  Default = ``nk`` (the
            body's wall-normal layer count).

        Returns
        -------
        list[HexBlock]
            For ``blocks=3``: ``[upper, body, lower]``.
            For ``blocks=4``: ``[upper, body, lower, wake_base]``.
            For ``blocks=1``: ``[cgrid]``.

            When ``cap_tips != "none"``: the C-grid blocks are
            followed by 5 cap blocks per requested tip
            (``cap_jlo_le``, ``cap_jlo_upper``, ``cap_jlo_lower``,
            ``cap_jlo_te_upper``, ``cap_jlo_te_lower`` for jlo;
            mirror for jhi), plus 1 extra ``cap_jX_wake_base`` block
            per tip when the body wrap is blunt (``blocks=4``).
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
                cap_te_frac=cap_te_frac,
                cap_tip_inset=cap_tip_inset,
                cap_n_cap_j=cap_n_cap_j,
                cap_j_spacing=cap_j_spacing,
                cap_j_beta=cap_j_beta,
                cap_nk=cap_nk,
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
        cap_te_frac: float,
        cap_tip_inset: float | None,
        cap_n_cap_j: int,
        cap_j_spacing: str,
        cap_j_beta: float,
        cap_nk: int | None,
    ) -> list[HexBlock]:
        """Build the conformal-cap C-grid: 4-block C-grid + cap blocks per tip.

        See :meth:`build` for the user-facing semantics of every parameter.
        """
        import sys

        from shore.surface.cgrid import round_body_te_in_surface
        from shore.topology.cgrid_cap import (
            build_cap_blocks,
            find_symmetric_le_cuts,
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
        upper, body, lower, wake_base = cgrid_blocks

        # ----- Build cap blocks per requested tip -----
        tips_to_build = (
            ("jlo", "jhi") if cap_tips == "both" else (cap_tips,)
        )
        cap_blocks: list[HexBlock] = []

        n_stations = body.nodes.shape[1]
        nk_cap = int(cap_nk) if cap_nk is not None else int(nk)

        if cap_tip_inset is None:
            tip_inset = 1e-3 * chord_len
        else:
            tip_inset = float(cap_tip_inset)

        for tip in tips_to_build:
            j_extreme = 0 if tip == "jlo" else n_stations - 1
            j_inboard = 1 if tip == "jlo" else n_stations - 2

            # Cap input: body cross-section at the spanwise tip, on the
            # airfoil wall (k=0).  Body has been rounded already, so this
            # cross-section has a non-zero TE thickness — no auto-rounding
            # needed inside the cap.
            cross_section = body.nodes[:, j_extreme, 0, :].copy()
            inboard_centroid = body.nodes[:, j_inboard, 0, :].mean(axis=0)
            extreme_centroid = cross_section.mean(axis=0)
            inboard_dir = inboard_centroid - extreme_centroid
            inboard_norm = float(np.linalg.norm(inboard_dir))
            if inboard_norm < 1e-12:
                raise ParameterError(
                    f"CGridTopology.build: at tip {tip!r}, the body's "
                    f"adjacent spanwise stations have coincident centroids — "
                    f"cannot infer outward span direction."
                )
            span_dir = -inboard_dir / inboard_norm

            # LE-cut and TE-cut indices on the (now-blunt) body cross-section.
            le_up, le_lo = find_symmetric_le_cuts(
                cross_section, le_frac=cap_le_frac, span_axis_world=span_axis_world
            )
            te_up, te_lo = find_symmetric_le_cuts(
                cross_section, le_frac=cap_te_frac, span_axis_world=span_axis_world
            )

            # Wake-base cross-section: spanwise-extreme face of the wake_base
            # block at this tip.  Shape (ni_wake, nk_w, 3) — exactly what
            # build_cap_blocks expects.
            wake_base_xs = wake_base.nodes[:, j_extreme, :, :].copy()

            cap_arrays = build_cap_blocks(
                cross_section=cross_section,
                le_cut_upper_idx=le_up,
                le_cut_lower_idx=le_lo,
                te_cut_upper_idx=te_up,
                te_cut_lower_idx=te_lo,
                span_dir=span_dir,
                tip_inset=tip_inset,
                n_cap_j=cap_n_cap_j,
                nk_cap=nk_cap,
                j_spacing=cap_j_spacing,
                j_beta=cap_j_beta,
                is_blunt=True,
                wake_base_cross_section=wake_base_xs,
            )

            # Restitch the airfoil-wall nodes (cap.j_hi at k=nk_cap-1) to be
            # bit-exact copies of the body's k=0 wall nodes at j_extreme.
            # build_cap_blocks already produces matching values from the
            # passed cross_section, but we re-copy here to guarantee bit-
            # exact agreement against the body block's actual storage even
            # if floating-point arithmetic in the builder introduced any
            # drift along i.
            cap_arrays["le"][:, -1, -1, :] = cross_section[
                np.arange(le_up, le_lo + 1, dtype=np.int64), :
            ]
            # cap_upper sweeps body i = le_up → te_up (reversed).
            cap_arrays["upper"][:, -1, -1, :] = cross_section[
                np.arange(le_up, te_up - 1, -1, dtype=np.int64), :
            ]
            # cap_lower sweeps body i = le_lo → te_lo (forward).
            cap_arrays["lower"][:, -1, -1, :] = cross_section[
                np.arange(le_lo, te_lo + 1, dtype=np.int64), :
            ]
            # cap_te_upper sweeps body i = te_up → 0 (reversed).
            cap_arrays["te_upper"][:, -1, -1, :] = cross_section[
                np.arange(te_up, -1, -1, dtype=np.int64), :
            ]
            # cap_te_lower sweeps body i = te_lo → ni-1 (forward).
            cap_arrays["te_lower"][:, -1, -1, :] = cross_section[
                np.arange(te_lo, cross_section.shape[0], dtype=np.int64), :
            ]

            for label_key, label in (
                ("le", "le"),
                ("upper", "upper"),
                ("lower", "lower"),
                ("te_upper", "te_upper"),
                ("te_lower", "te_lower"),
                ("wake_base", "wake_base"),
            ):
                blk = HexBlock(
                    label=f"cap_{tip}_{label}",
                    nodes=cap_arrays[label_key],
                    periodic_i=False,
                    periodic_j=False,
                    spacing_k=spacing_k,
                )
                blk._set_face_bc("k_lo", FaceBC.WALL)
                blk._set_face_bc("j_hi", FaceBC.WALL)
                cap_blocks.append(blk)

            # Wire the cap's internal SHARED seams.  Block layout in the
            # last 6 entries of cap_blocks (per tip): le, upper, lower,
            # te_upper, te_lower, wake_base.
            cap_le_blk = cap_blocks[-6]
            cap_upper_blk = cap_blocks[-5]
            cap_lower_blk = cap_blocks[-4]
            cap_te_upper_blk = cap_blocks[-3]
            cap_te_lower_blk = cap_blocks[-2]

            # cap_le.i_lo ↔ cap_upper.i_lo  (LE-cut, suction side)
            cap_le_blk._set_face_bc(
                "i_lo", FaceBC.SHARED,
                partner=cap_upper_blk.face("i_lo"),
                axis_map=AxisMap.identity("i_lo"),
            )
            cap_upper_blk._set_face_bc(
                "i_lo", FaceBC.SHARED,
                partner=cap_le_blk.face("i_lo"),
                axis_map=AxisMap.identity("i_lo"),
            )
            # cap_le.i_hi ↔ cap_lower.i_lo  (LE-cut, pressure side)
            cap_le_blk._set_face_bc(
                "i_hi", FaceBC.SHARED,
                partner=cap_lower_blk.face("i_lo"),
                axis_map=AxisMap.identity("i_hi"),
            )
            cap_lower_blk._set_face_bc(
                "i_lo", FaceBC.SHARED,
                partner=cap_le_blk.face("i_hi"),
                axis_map=AxisMap.identity("i_lo"),
            )
            # cap_upper.j_lo ↔ cap_lower.j_lo  (chord cut)
            cap_upper_blk._set_face_bc(
                "j_lo", FaceBC.SHARED,
                partner=cap_lower_blk.face("j_lo"),
                axis_map=AxisMap.identity("j_lo"),
            )
            cap_lower_blk._set_face_bc(
                "j_lo", FaceBC.SHARED,
                partner=cap_upper_blk.face("j_lo"),
                axis_map=AxisMap.identity("j_lo"),
            )
            # cap_upper.i_hi ↔ cap_te_upper.i_lo  (TE-cut, suction side)
            cap_upper_blk._set_face_bc(
                "i_hi", FaceBC.SHARED,
                partner=cap_te_upper_blk.face("i_lo"),
                axis_map=AxisMap.identity("i_hi"),
            )
            cap_te_upper_blk._set_face_bc(
                "i_lo", FaceBC.SHARED,
                partner=cap_upper_blk.face("i_hi"),
                axis_map=AxisMap.identity("i_lo"),
            )
            # cap_lower.i_hi ↔ cap_te_lower.i_lo  (TE-cut, pressure side)
            cap_lower_blk._set_face_bc(
                "i_hi", FaceBC.SHARED,
                partner=cap_te_lower_blk.face("i_lo"),
                axis_map=AxisMap.identity("i_hi"),
            )
            cap_te_lower_blk._set_face_bc(
                "i_lo", FaceBC.SHARED,
                partner=cap_lower_blk.face("i_hi"),
                axis_map=AxisMap.identity("i_lo"),
            )
            # cap_te_upper.j_lo ↔ cap_te_lower.j_lo  (chord-cut continuation)
            cap_te_upper_blk._set_face_bc(
                "j_lo", FaceBC.SHARED,
                partner=cap_te_lower_blk.face("j_lo"),
                axis_map=AxisMap.identity("j_lo"),
            )
            cap_te_lower_blk._set_face_bc(
                "j_lo", FaceBC.SHARED,
                partner=cap_te_upper_blk.face("j_lo"),
                axis_map=AxisMap.identity("j_lo"),
            )

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
        nodes = np.empty((ni_total, nj, nk, 3), dtype=np.float64)
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
        nodes_all = np.empty((ni_total, nj, nk, 3), dtype=np.float64)
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
        nodes_all = np.empty((ni_total, nj, nk, 3), dtype=np.float64)
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
        wake_base_nodes = np.empty((ni_wake, nj, nk_w, 3), dtype=np.float64)
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
    """C-grid wing-tip cap topology — 3 (sharp) or 4 (blunt) overset blocks.

    Builds the H-style tip closure for one tip of a C-grid wing as an
    independent chimera component.  The cap is composed of 3 (sharp TE)
    or 4 (blunt TE) hexahedral sub-blocks that fill the airfoil-interior
    cross-section at the wing-tip plane in an "H3 T-fill":

    - ``cap_le``    — wraps the leading edge from the upper-LE-cut to the
                       lower-LE-cut.
    - ``cap_upper`` — chord-band on the suction side, between LE-cut and
                       upper-TE-cut, above the chord cut.
    - ``cap_lower`` — mirror of ``cap_upper`` on the pressure side.
    - ``cap_wake_base`` (blunt only) — closes the wake_base block's
                       spanwise face downstream of the TE base segment.

    All sub-blocks share the same axis convention: i = chordwise along
    the airfoil profile (or LE wrap), j = airfoil-thickness direction
    (chord cut → wall), k = spanwise.  k=0 is the outboard wing-tip
    wall (offset from the cross-section by ``tip_inset * span_dir``);
    k=nk-1 is the inboard chimera-overlap plane (coincident with the
    body wrap's spanwise-extreme cross-section).

    Default face BCs:

    - **k_lo**: WALL (the actual wing-tip wall closure).
    - **k_hi**: FREE (chimera fringe receiving from the body wrap).
    - **j_hi**: WALL (the airfoil suction / pressure / LE / TE-base
      surface, mirrored into the cap).
    - **j_lo**: SHARED between cap_upper and cap_lower (chord cut),
      otherwise FREE.
    - **i_lo / i_hi**: SHARED at the LE-cuts between cap_le and
      cap_upper / cap_lower; FREE otherwise (TE-cut chord-perp
      sliver).

    The cap is **fully volumetric at build time** — calling
    :meth:`shore.mesh.Mesh.march` on a CGridTipCap mesh is a no-op.
    The cap is meant to be assembled against a C-grid body wrap as a
    chimera component via ``shore cc-par``.

    Parameters
    ----------
    le_cut_upper_idx, le_cut_lower_idx :
        i-indices in the ``surface`` cross-section polyline at the
        LE-cut on the suction and pressure sides.  Must satisfy
        ``0 < le_cut_upper_idx < le_cut_lower_idx < ni - 1``.
    span_dir :
        Length-3 unit vector pointing **outward** along the span axis
        at this tip.
    tip_inset :
        Spanwise thickness of the cap slab.  > 0.
    n_cap_j :
        Node count along the cap's j-axis (airfoil-thickness).  >= 2.
    cap_j_spacing :
        ``"uniform"`` / ``"tanh"`` / ``"tanh2"``.  ``"tanh"`` (default)
        clusters near the airfoil wall.
    cap_j_beta :
        Clustering strength.
    is_blunt :
        Toggles the 4th sub-block ``cap_wake_base``.
    wake_base_cross_section :
        Required when ``is_blunt=True``.  ``(ni_wake, nk_w, 3)`` —
        the wake_base block's spanwise-extreme face from the C-grid.

    Future input modes (not implemented yet)
    ----------------------------------------
    The current ``build()`` accepts a 1D cross-section polyline as
    ``surface``.  In future, ``build()`` will additionally accept a 2D
    sliced STL/STEP cross-section ``(ni, nj_thin, 3)`` so the cap can
    be driven directly by a CAD/STL surface input — the same way the
    body wrap is.  The signature remains the same.
    """

    def __init__(
        self,
        le_cut_upper_idx: int,
        le_cut_lower_idx: int,
        span_dir: tuple[float, float, float],
        tip_inset: float,
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
    ) -> None:
        # Validate span_dir up front.
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
        if tip_inset <= 0:
            raise ParameterError(
                f"CGridTipCap: tip_inset must be > 0, got {tip_inset}"
            )
        if n_cap_j < 2:
            raise ParameterError(
                f"CGridTipCap: n_cap_j must be >= 2, got {n_cap_j}"
            )
        if cap_j_spacing not in ("uniform", "tanh", "tanh2"):
            raise ParameterError(
                f"CGridTipCap: cap_j_spacing must be 'uniform' / 'tanh' / "
                f"'tanh2', got {cap_j_spacing!r}"
            )
        if is_blunt and wake_base_cross_section is None:
            raise ParameterError(
                "CGridTipCap: is_blunt=True requires wake_base_cross_section."
            )
        if te_thickness is not None and te_thickness_frac is not None:
            raise ParameterError(
                "CGridTipCap: pass either te_thickness (absolute) or "
                "te_thickness_frac (fraction of chord), not both."
            )
        if te_thickness is not None and te_thickness <= 0:
            raise ParameterError(
                f"CGridTipCap: te_thickness must be > 0, got {te_thickness}"
            )
        if te_thickness_frac is not None and te_thickness_frac <= 0:
            raise ParameterError(
                f"CGridTipCap: te_thickness_frac must be > 0, "
                f"got {te_thickness_frac}"
            )
        self.le_cut_upper_idx = int(le_cut_upper_idx)
        self.le_cut_lower_idx = int(le_cut_lower_idx)
        self.te_cut_upper_idx = (
            None if te_cut_upper_idx is None else int(te_cut_upper_idx)
        )
        self.te_cut_lower_idx = (
            None if te_cut_lower_idx is None else int(te_cut_lower_idx)
        )
        self.span_dir = span_arr
        self.tip_inset = float(tip_inset)
        self.n_cap_j = int(n_cap_j)
        self.cap_j_spacing = cap_j_spacing
        self.cap_j_beta = float(cap_j_beta)
        self.is_blunt = bool(is_blunt)
        self.wake_base_cross_section = (
            np.asarray(wake_base_cross_section, dtype=np.float64)
            if wake_base_cross_section is not None
            else None
        )
        self.auto_round_te = bool(auto_round_te)
        self.te_thickness = (
            None if te_thickness is None else float(te_thickness)
        )
        self.te_thickness_frac = (
            None if te_thickness_frac is None else float(te_thickness_frac)
        )

    def build(
        self,
        surface: np.ndarray,
        nk: int,
        spacing_k: Spacing1D | None = None,
        **kwargs,
    ) -> list[HexBlock]:
        """Build the tip cap blocks from a cross-section polyline.

        Parameters
        ----------
        surface :
            ``(ni_body, 3)`` cross-section polyline at the wing tip.
            Sweeps upper-TE-cut → suction → LE → pressure →
            lower-TE-cut.
        nk :
            Spanwise node count for the cap (``nk_cap``).  Must be >= 2.
        spacing_k :
            Reserved for future per-block k-spacing overrides.
            Currently the cap's k-axis is uniformly spaced from the
            outboard wall to the cross-section plane.

        Returns
        -------
        list[HexBlock]
            ``[cap_le, cap_upper, cap_lower]`` for sharp TE, or
            ``[cap_le, cap_upper, cap_lower, cap_wake_base]`` for blunt.
        """
        import sys

        from shore.topology.cgrid_cap import (
            auto_round_sharp_te,
            build_cap_blocks,
            detect_te_thickness,
        )

        if nk < 2:
            raise ParameterError(
                f"CGridTipCap.build: nk must be >= 2, got {nk}"
            )

        # ------------------------------------------------------------------
        # Sharp/blunt detection + auto-round.
        # ------------------------------------------------------------------
        # Identify the spanwise world axis from the *outward* span_dir
        # (largest absolute component).
        span_axis_world = int(np.argmax(np.abs(self.span_dir)))

        existing_te = detect_te_thickness(surface, span_axis_world)
        cs_for_cap = surface
        if existing_te < 1e-9:
            # Sharp input.  Auto-round (default) or refuse.
            if not self.auto_round_te:
                raise ParameterError(
                    "CGridTipCap.build: cross-section has a sharp TE "
                    "(i=0 and i=ni-1 coincide) but auto_round_te=False; "
                    "the H-fill cap topology has a singular axis on a "
                    "sharp TE.  Either enable auto_round_te (default) "
                    "or pre-round the TE in the source body."
                )
            # Resolve the synthetic TE thickness.
            if self.te_thickness is not None:
                te = self.te_thickness
            else:
                # te_thickness_frac of the airfoil's chord extent.
                extents = surface.max(axis=0) - surface.min(axis=0)
                extents[span_axis_world] = -np.inf
                chord_len = float(extents.max())
                frac = self.te_thickness_frac if self.te_thickness_frac is not None else 0.002
                te = frac * chord_len
            cs_for_cap = auto_round_sharp_te(
                surface, te_thickness=te, span_axis_world=span_axis_world,
            )
            sys.stderr.write(
                f"CGridTipCap: sharp TE detected; auto-rounding cap with "
                f"te_thickness={te:.6g} (the body wrap is unchanged).\n"
            )
        # Otherwise: input is already blunt → use surface as-is.

        cap_arrays = build_cap_blocks(
            cross_section=cs_for_cap,
            le_cut_upper_idx=self.le_cut_upper_idx,
            le_cut_lower_idx=self.le_cut_lower_idx,
            te_cut_upper_idx=self.te_cut_upper_idx,
            te_cut_lower_idx=self.te_cut_lower_idx,
            span_dir=self.span_dir,
            tip_inset=self.tip_inset,
            n_cap_j=self.n_cap_j,
            nk_cap=int(nk),
            j_spacing=self.cap_j_spacing,
            j_beta=self.cap_j_beta,
            is_blunt=self.is_blunt,
            wake_base_cross_section=self.wake_base_cross_section,
        )
        _sp_k = spacing_k or Spacing1D(law="geometric", ds=1e-3)

        cap_le = HexBlock(
            label="cap_le",
            nodes=cap_arrays["le"],
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        cap_upper = HexBlock(
            label="cap_upper",
            nodes=cap_arrays["upper"],
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        cap_lower = HexBlock(
            label="cap_lower",
            nodes=cap_arrays["lower"],
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        cap_te_upper = HexBlock(
            label="cap_te_upper",
            nodes=cap_arrays["te_upper"],
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        cap_te_lower = HexBlock(
            label="cap_te_lower",
            nodes=cap_arrays["te_lower"],
            periodic_i=False,
            periodic_j=False,
            spacing_k=_sp_k,
        )
        blocks: list[HexBlock] = [
            cap_le, cap_upper, cap_lower, cap_te_upper, cap_te_lower,
        ]

        if self.is_blunt:
            cap_wake_base = HexBlock(
                label="cap_wake_base",
                nodes=cap_arrays["wake_base"],
                periodic_i=False,
                periodic_j=False,
                spacing_k=_sp_k,
            )
            blocks.append(cap_wake_base)

        # ------------------------------------------------------------
        # Set face BCs.  Defaults (FREE / WALL) per the class docstring;
        # then wire SHARED seams between sub-blocks.
        # ------------------------------------------------------------
        # k_lo: outboard wing-tip wall (WALL) on every cap block.
        # k_hi: inboard chimera-overlap plane (FREE) on every cap block).
        # j_hi: airfoil wall (WALL) on every chord-band cap block.
        for blk in (cap_le, cap_upper, cap_lower, cap_te_upper, cap_te_lower):
            blk._set_face_bc("k_lo", FaceBC.WALL)
            blk._set_face_bc("j_hi", FaceBC.WALL)
            # k_hi remains FREE (default).

        if self.is_blunt:
            cap_wake_base._set_face_bc("k_lo", FaceBC.WALL)
            # cap_wake_base.j_hi: the airfoil's TE-base wall (a flat
            # rectangular strip).  WALL.
            cap_wake_base._set_face_bc("j_hi", FaceBC.WALL)

        # ------------------------------------------------------------
        # SHARED seams between cap sub-blocks.
        # ------------------------------------------------------------
        # 1. cap_le.i_lo ↔ cap_upper.i_lo (LE-cut, suction side).
        #    Both columns are at the upper-LE-cut chord position with j
        #    running from chord-perp / chord-cut (j=0) to airfoil wall
        #    (j=nj-1) and k spanning the cap thickness.  Identity AxisMap.
        cap_le._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=cap_upper.face("i_lo"),
            axis_map=AxisMap.identity("i_lo"),
        )
        cap_upper._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=cap_le.face("i_lo"),
            axis_map=AxisMap.identity("i_lo"),
        )
        # 2. cap_le.i_hi ↔ cap_lower.i_lo (LE-cut, pressure side).
        cap_le._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=cap_lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        cap_lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=cap_le.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        # 3. cap_upper.j_lo ↔ cap_lower.j_lo (chord cut at y=0,
        #    chord-band region).  Both blocks have j=0 = chord cut,
        #    j=nj-1 = wall (suction for upper, pressure for lower).
        #    i-axes both run LE-cut → TE-cut.  Identity AxisMap on (i, k).
        cap_upper._set_face_bc(
            "j_lo",
            FaceBC.SHARED,
            partner=cap_lower.face("j_lo"),
            axis_map=AxisMap.identity("j_lo"),
        )
        cap_lower._set_face_bc(
            "j_lo",
            FaceBC.SHARED,
            partner=cap_upper.face("j_lo"),
            axis_map=AxisMap.identity("j_lo"),
        )
        # 4. cap_upper.i_hi ↔ cap_te_upper.i_lo (TE-cut, suction side).
        cap_upper._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=cap_te_upper.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        cap_te_upper._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=cap_upper.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        # 5. cap_lower.i_hi ↔ cap_te_lower.i_lo (TE-cut, pressure side).
        cap_lower._set_face_bc(
            "i_hi",
            FaceBC.SHARED,
            partner=cap_te_lower.face("i_lo"),
            axis_map=AxisMap.identity("i_hi"),
        )
        cap_te_lower._set_face_bc(
            "i_lo",
            FaceBC.SHARED,
            partner=cap_lower.face("i_hi"),
            axis_map=AxisMap.identity("i_lo"),
        )
        # 6. cap_te_upper.j_lo ↔ cap_te_lower.j_lo (chord-cut continuation
        #    in the TE region).  Both have j=0 = chord cut, j=nj-1 = wall.
        cap_te_upper._set_face_bc(
            "j_lo",
            FaceBC.SHARED,
            partner=cap_te_lower.face("j_lo"),
            axis_map=AxisMap.identity("j_lo"),
        )
        cap_te_lower._set_face_bc(
            "j_lo",
            FaceBC.SHARED,
            partner=cap_te_upper.face("j_lo"),
            axis_map=AxisMap.identity("j_lo"),
        )
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
            raise ParameterError(
                f"OHTopology: nk_channel must be >= 2, got {nk_channel}"
            )
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
                raise ParameterError(
                    "OHTopology.build needs ni and nj (or a surface array)."
                )
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
            ni=ni, nj=nj, nk=nk_wrap, total_k=total_k,
            cap_proj=cap_proj, spacing_k=spacing_k,
            theta_cap=theta_cap_oh,
            i_spacing=i_spacing, i_beta=i_beta,
            j_spacing=j_spacing, j_beta=j_beta,
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
            label = (
                f"chan_sub{q}" if q < 4 else
                "chan_north" if q == 4 else
                "chan_south"
            )
            chan_blk = HexBlock(
                label=label, nodes=chan_view,
                periodic_i=False, periodic_j=False,
                spacing_k=spacing_k or Spacing1D(law="geometric", ds=1e-3),
            )
            channel_blocks.append(chan_blk)

        # Wire wrap.k_hi <-> channel.k_lo (SHARED, view-aliased, identity).
        for w_blk, c_blk in zip(wrap_blocks, channel_blocks, strict=True):
            w_blk._set_face_bc(
                "k_hi", FaceBC.SHARED, partner=c_blk.face("k_lo"),
                axis_map=AxisMap.identity("k_hi"),
            )
            c_blk._set_face_bc(
                "k_lo", FaceBC.SHARED, partner=w_blk.face("k_hi"),
                axis_map=AxisMap.identity("k_lo"),
            )

        # Wire channel-channel lateral seams to mirror the wrap topology.
        # (sub_q.j_hi ↔ sub_(q+1).j_lo for q in 0..3 → same on chan side)
        chan_subs = channel_blocks[:4]
        for q in range(4):
            q_next = (q + 1) % 4
            chan_subs[q]._set_face_bc(
                "j_hi", FaceBC.SHARED, partner=chan_subs[q_next].face("j_lo"),
                axis_map=AxisMap.identity("j_hi"),
            )
            chan_subs[q_next]._set_face_bc(
                "j_lo", FaceBC.SHARED, partner=chan_subs[q].face("j_hi"),
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
                    cap_face, FaceBC.DIRICHLET,
                    partner=csub.face(sub_face_lit),
                    axis_map=sub_axis_map.inverse(sub_face_lit, cap_face),
                )
                csub._set_face_bc(
                    sub_face_lit, FaceBC.DIRICHLET,
                    partner=chan_cap_blk.face(cap_face),
                    axis_map=sub_axis_map,
                )

        # k_hi on every channel block defaults to FREE (HexBlock constructor).
        # No further wiring needed.

        return [*wrap_blocks, *channel_blocks]

    def fill_channel_blocks(self, blocks: list[HexBlock]) -> None:
        """Build channel block interiors after the wrap has been marched.

        Walks the 6 wrap–channel pairs, computes each channel's outer face
        on the user-supplied bounding box, and fills the channel's
        ``k=1..nk_channel-1`` interior by ruled-surface frustum interpolation.

        The wrap–channel seam (channel k=0 = wrap k_hi) is already populated
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

        for w_blk, c_blk in zip(wrap_blocks, channel_blocks, strict=True):
            inner = w_blk.nodes[:, :, -1, :]  # (ni, nj, 3) — wrap envelope
            outer_quad = self._outer_quad_for_block(inner)
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

    def _outer_quad_for_block(self, inner_face: np.ndarray) -> np.ndarray:
        """Find the box face this wrap block points at; return its 4 corners.

        Algorithm: pick the cube face from the **full** envelope centroid
        (mean of all i, j envelope nodes — robust because corners sit at
        cap-equator junctions where the bulk is not).  Then snap each of
        the 4 envelope corners onto the cube by:

        1. Setting the face-normal component to the picked plane value.
        2. For each in-plane component, taking the sign of the corner's
           coordinate (around the box centre) and pushing it to the
           corresponding ``vmin`` / ``vmax``.

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
        if not np.all(
            (full_centroid >= self.box_vmin) & (full_centroid <= self.box_vmax)
        ):
            raise ParameterError(
                f"OHTopology: wrap envelope centroid {full_centroid} is "
                f"outside the bounding box [{self.box_vmin}, {self.box_vmax}].  "
                f"The box must contain the marched wrap envelope; either "
                f"enlarge the box or reduce nk / wrap step sizes."
            )

        dists = np.array([
            self.box_vmax[0] - full_centroid[0],  # +x face
            full_centroid[0] - self.box_vmin[0],  # -x face
            self.box_vmax[1] - full_centroid[1],  # +y face
            full_centroid[1] - self.box_vmin[1],  # -y face
            self.box_vmax[2] - full_centroid[2],  # +z face
            full_centroid[2] - self.box_vmin[2],  # -z face
        ])
        face_idx = int(np.argmin(dists))
        axis_idx = face_idx // 2
        plane_val = (self.box_vmax if face_idx % 2 == 0 else self.box_vmin)[axis_idx]

        box_centre = 0.5 * (self.box_vmin + self.box_vmax)

        out_corners = np.stack([c00, c10, c11, c01], axis=0).copy()
        out_corners[:, axis_idx] = plane_val
        for ax in range(3):
            if ax == axis_idx:
                continue
            # Push each corner's in-plane component to vmin or vmax based on
            # its sign relative to the box centre.  This snaps every corner
            # to a cube vertex of the picked face.
            sign = np.sign(out_corners[:, ax] - box_centre[ax])
            # Treat exactly-on-axis corners (sign == 0) as positive; in
            # practice they don't occur for cubed-sphere envelopes but the
            # snap must be deterministic.
            sign = np.where(sign == 0, 1.0, sign)
            out_corners[:, ax] = np.where(
                sign > 0, self.box_vmax[ax], self.box_vmin[ax]
            )
        return out_corners


def _projector_from_mesh(mesh: object) -> object:
    """Build a temporary Projector using mesh.centroid (legacy mesh= path)."""
    from shore.surface.projector import Projector

    return Projector(mesh, mesh.centroid)  # type: ignore[attr-defined]
