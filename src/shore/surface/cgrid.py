"""C-grid surface construction: TE detection, span slicing, C-curve build.

A C-grid is a single-block structured mesh that wraps around a body
with a sharp trailing edge — typically an airfoil or hydrofoil — by
cutting the topological cylinder along the wake.  Three pieces of
geometry processing turn an STL into the C-grid k=0 layer:

1. :func:`detect_te_edges` — find the TE polyline via dihedral-angle
   classification on triangle-pair edges.
2. :func:`slice_by_span` — slice the STL with planes perpendicular
   to the span axis, producing one cross-section polyline per
   station, ordered upper-TE → LE → lower-TE.
3. :func:`build_c_curve` — add upper + lower wake branches to the
   per-station body polyline, returning the (ni_total, 2) C-curve.

The orchestrator :func:`build_cgrid_k0` chains them together to
produce the (ni, n_stations, 3) k=0 layer that
:class:`shore.mesh.topology.CGridTopology` consumes.

Conventions
-----------
- ``span_axis`` defaults to ``"z"``; the body extends from
  ``span_axis_min`` to ``span_axis_max`` along this axis.
- The C-curve i-direction starts at the **upper TE** (i=0), sweeps
  forward over the suction surface, around the LE, back over the
  pressure surface, ends at the **lower TE** (i=ni_body-1).  The
  upper wake branch then extends i from ``-ni_wake`` to ``0``
  (concatenated below the body wrap in i-space, so total i runs
  ``[upper-wake-far, upper-TE, ..., LE, ..., lower-TE, lower-wake-far]``).
- Wake branches extend in the **global freestream direction**
  (default ``+x_world``), regardless of local TE orientation.
- Per-station polylines are re-parameterised by **arclength fraction**
  so i = 0.5 sits at the LE on every station, regardless of chord
  variation along the span (tapered / swept wings).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import trimesh

from shore.exceptions import MeshInputError, ParameterError

_AXIS_NAMES = ("x", "y", "z")
_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------


@dataclass
class TEDetectionResult:
    """Outcome of :func:`detect_te_edges`.

    Attributes
    ----------
    polylines :
        List of ordered ``(n_pts, 3)`` arrays of TE vertex
        coordinates, each sorted by ``span_axis``.  A **sharp** TE
        produces one polyline; a **blunt** TE produces two
        (the upper-edge and lower-edge feature lines of the
        finite-thickness TE strip).  Order: by descending span
        extent.
    is_blunt :
        ``True`` iff ``len(polylines) >= 2`` and the two strands are
        within ``blunt_te_max_thickness * chord`` of each other (i.e.
        they're the upper / lower of the same TE, not unrelated
        features).  Consumers use this flag to decide whether to
        treat the section as sharp- or blunt-TE.
    te_thickness :
        For blunt TEs (``is_blunt == True``), the average distance
        between the two strands at matched span coordinates.  ``0.0``
        for sharp TEs.
    span_axis :
        The axis label (``"x"`` / ``"y"`` / ``"z"``) the polylines
        are sorted on.
    max_dihedral_deg :
        Largest dihedral angle observed in the input mesh, in degrees.
        Useful for diagnostics when no TE is found.
    n_candidate_groups :
        Number of disconnected span-aligned candidate groups
        discovered.
    """

    polylines: list[np.ndarray]
    is_blunt: bool
    te_thickness: float
    span_axis: str
    max_dihedral_deg: float
    n_candidate_groups: int

    @property
    def polyline(self) -> np.ndarray:
        """Convenience: the longest TE polyline.

        For sharp TEs this is the only polyline; for blunt TEs it's
        either the upper or lower strand (consumers that don't care
        about blunt TEs can use this directly).
        """
        return self.polylines[0]


# ---------------------------------------------------------------------------
# TE detection
# ---------------------------------------------------------------------------


def detect_te_edges(
    mesh: trimesh.Trimesh,
    *,
    dihedral_deg: float = 30.0,
    span_axis: str = "z",
) -> TEDetectionResult:
    """Find the trailing-edge polyline of a body STL by dihedral classification.

    For every edge shared by exactly two triangles, compute the
    dihedral angle between the two face normals.  An edge whose
    dihedral exceeds ``dihedral_deg`` (i.e. the body folds by more
    than ``180° - dihedral_deg`` flat-line at that edge) is a TE
    candidate.  Candidates are clustered by vertex connectivity; the
    largest connected group is taken as the trailing edge.  The
    resulting polyline is ordered by ``span_axis`` coordinate.

    Parameters
    ----------
    mesh :
        Watertight ``trimesh.Trimesh`` of the body.
    dihedral_deg :
        Sharpness threshold in degrees.  An edge is a TE candidate
        when ``angle(face_normals) >= dihedral_deg``.  Default
        ``30°`` (≈ feature line at body curvature).  Lower (e.g.
        ``20°``) for blunter TEs; higher (e.g. ``45°``) to reject
        spurious tessellation seams.
    span_axis :
        ``"x"``, ``"y"``, or ``"z"``.  The polyline is sorted along
        this axis at the end.  Default ``"z"``.

    Returns
    -------
    TEDetectionResult
        ``polyline`` is the ordered ``(n_te_pts, 3)`` array of TE
        vertices.

    Raises
    ------
    MeshInputError
        If no edge exceeds the threshold (body is too smooth — try
        a smaller ``dihedral_deg`` or supply the TE manually), or if
        the largest candidate group is too small (< 2 points) to
        form a polyline, or if the resulting polyline isn't
        approximately aligned with ``span_axis``.
    ParameterError
        On invalid arguments.
    """
    if span_axis not in _AXIS_INDEX:
        raise ParameterError(
            f"detect_te_edges: span_axis must be 'x' / 'y' / 'z', got {span_axis!r}"
        )
    if not (0.0 < dihedral_deg < 180.0):
        raise ParameterError(
            f"detect_te_edges: dihedral_deg must be in (0, 180), got {dihedral_deg}"
        )
    span_idx = _AXIS_INDEX[span_axis]

    # Build the edge → adjacent-faces map.  trimesh exposes face_adjacency
    # which already gives us pairs of faces sharing exactly one edge.
    face_adj = mesh.face_adjacency  # (m, 2): pairs of face indices
    if len(face_adj) == 0:
        raise MeshInputError(
            "detect_te_edges: mesh has no shared-edge face pairs; cannot detect a TE"
        )
    adj_edges = mesh.face_adjacency_edges  # (m, 2): vertex pairs along each shared edge
    face_normals = mesh.face_normals  # (n_faces, 3)

    # Dihedral angle per shared edge: angle between the two face normals.
    n0 = face_normals[face_adj[:, 0]]
    n1 = face_normals[face_adj[:, 1]]
    cos_th = np.einsum("ij,ij->i", n0, n1).clip(-1.0, 1.0)
    angles_deg = np.degrees(np.arccos(cos_th))  # 0 = coplanar, 180 = back-to-back
    max_dihedral = float(angles_deg.max())

    # TE candidate edges: angle >= threshold.
    te_mask = angles_deg >= dihedral_deg
    if not te_mask.any():
        raise MeshInputError(
            f"detect_te_edges: no edge dihedral exceeds {dihedral_deg}° "
            f"(max found: {max_dihedral:.2f}°).  The body has no sharp feature; "
            "lower --te-dihedral-deg, or supply the TE polyline manually via "
            "--te-curve."
        )

    candidate_edges = adj_edges[te_mask]  # (k, 2): vertex index pairs

    # Filter candidates by edge-direction: a TE runs *along the span
    # axis*, while spurious feature lines (LE-projection ridges,
    # cap-fan rays, coordinate-axis tessellation seams) run
    # perpendicular to span.  Keep only edges whose direction is
    # within ~30° of the span axis.  The threshold is wide enough to
    # tolerate a slightly skew TE (swept wing) but tight enough to
    # reject LE / cap edges.
    edge_dirs = mesh.vertices[candidate_edges[:, 1]] - mesh.vertices[candidate_edges[:, 0]]
    edge_lengths = np.linalg.norm(edge_dirs, axis=1)
    edge_dirs_unit = edge_dirs / np.maximum(edge_lengths, 1e-15)[:, None]
    span_dir = np.zeros(3)
    span_dir[span_idx] = 1.0
    cos_span = np.abs(edge_dirs_unit @ span_dir)  # |cos|: direction-agnostic
    # Edges within ~30° of the span axis: |cos| >= cos(30°) ≈ 0.866.
    span_aligned = cos_span >= np.cos(np.radians(30.0))
    if not span_aligned.any():
        raise MeshInputError(
            f"detect_te_edges: no candidate edge is aligned with span axis "
            f"{span_axis!r} (best alignment {np.degrees(np.arccos(cos_span.max())):.1f}°).  "
            "The TE may not be parallel to the span axis — check --span-axis "
            "or supply --te-curve manually."
        )
    candidate_edges = candidate_edges[span_aligned]

    # Cluster the (now span-aligned) candidates by vertex connectivity
    # (union-find).  Vertices that share a candidate edge belong to
    # the same group.
    n_verts = len(mesh.vertices)
    parent = np.arange(n_verts, dtype=np.int64)

    def _find(v: int) -> int:
        # Path-compressing find.
        root = v
        while parent[root] != root:
            root = int(parent[root])
        while parent[v] != root:
            v, parent[v] = int(parent[v]), root
        return root

    def _union(a: int, b: int) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[ra] = rb

    for v0, v1 in candidate_edges:
        _union(int(v0), int(v1))

    # Group candidate vertices by their union-find root.
    candidate_verts = np.unique(candidate_edges.flatten())
    groups: dict[int, list[int]] = {}
    for v in candidate_verts:
        r = _find(int(v))
        groups.setdefault(r, []).append(int(v))

    n_groups = len(groups)

    # Sort groups by descending span extent.  Real TE strands span the
    # full body; spurious tessellation seams typically span a short
    # range.
    def _span_extent(verts: np.ndarray) -> float:
        coords = mesh.vertices[verts]
        return float(coords[:, span_idx].max() - coords[:, span_idx].min())

    sorted_groups = sorted(
        groups.values(),
        key=lambda v: -_span_extent(np.asarray(v, dtype=np.int64)),
    )

    # Keep only groups whose span extent is at least 80% of the longest's.
    # This admits the matching pair of strands on a blunt TE (both span
    # the whole body) while rejecting short tessellation seams.
    longest_extent = _span_extent(np.asarray(sorted_groups[0], dtype=np.int64))
    if longest_extent < 1e-12:
        raise MeshInputError(
            "detect_te_edges: longest span-aligned group has zero span "
            "extent — the TE detection found no candidates spanning the body."
        )
    keep_threshold = 0.8 * longest_extent
    selected_groups = [
        np.asarray(g, dtype=np.int64)
        for g in sorted_groups
        if _span_extent(np.asarray(g, dtype=np.int64)) >= keep_threshold
    ]

    # Of the kept groups, the TE proper sits at one **chord-extreme**:
    # for an airfoil, this is the maximum-x position (or whichever
    # non-span axis the freestream points along).  Only keep groups
    # whose mean chord position is within 0.5% of body extent of the
    # maximum.  This rejects LE-curvature smoothing-edge artefacts on
    # cosine-spaced sections (which appear at the *second-most-extreme*
    # chord position) while admitting genuine sharp TEs at the
    # literal chord extremum.
    body_spreads = mesh.bounds[1] - mesh.bounds[0]
    non_span_axes = [a for a in (0, 1, 2) if a != span_idx]
    chord_axis = (
        non_span_axes[0]
        if body_spreads[non_span_axes[0]] >= body_spreads[non_span_axes[1]]
        else non_span_axes[1]
    )

    def _chord_position(verts: np.ndarray) -> float:
        return float(mesh.vertices[verts][:, chord_axis].mean())

    # Filter by literal maximum chord coordinate.  Per group, take the
    # **maximum** vertex chord coord (not the mean — the TE strand on
    # a slightly tilted body might have one or two vertices off the
    # extremum but still be the right strand).  Among groups, only
    # keep those whose max chord is within 0.05% of the global max:
    # this rejects LE-curvature smoothing-edge artefacts at slightly
    # lower chord positions (cosine-spaced section: artefact at the
    # second-most-extreme chord, ~0.2% of body extent in).
    def _max_chord(verts: np.ndarray) -> float:
        return float(mesh.vertices[verts][:, chord_axis].max())

    if len(selected_groups) > 1:
        max_chords = [_max_chord(g) for g in selected_groups]
        global_max = max(max_chords)
        chord_tolerance = 5e-4 * float(body_spreads[chord_axis])
        selected_groups = [
            g
            for g, mc in zip(selected_groups, max_chords, strict=True)
            if abs(mc - global_max) <= chord_tolerance
        ]

    # Build sorted polylines per kept group.
    polylines: list[np.ndarray] = []
    for g in selected_groups:
        coords = mesh.vertices[g]
        order = np.argsort(coords[:, span_idx])
        polylines.append(coords[order])

    if any(len(p) < 2 for p in polylines):
        raise MeshInputError("detect_te_edges: a kept TE candidate has fewer than 2 vertices.")

    # Sanity: the longest polyline must be roughly parallel to ``span_axis``.
    p0 = polylines[0]
    spreads = p0.max(axis=0) - p0.min(axis=0)
    dominant_axis = int(np.argmax(spreads))
    if dominant_axis != span_idx:
        raise MeshInputError(
            f"detect_te_edges: detected TE is most extended along axis "
            f"{_AXIS_NAMES[dominant_axis]!r} (spreads {spreads.tolist()}), "
            f"not the requested span axis {span_axis!r}.  Re-run with "
            f"--span-axis {_AXIS_NAMES[dominant_axis]} or rotate the STL."
        )

    # Blunt-TE detection: at least 2 strands, both spanning the body,
    # separated by a small (< 1% of the longest body extent) distance
    # in the cross-section.  The two strands are roughly co-located in
    # the chord direction (= the dominant *non*-span direction).
    is_blunt = False
    te_thickness = 0.0
    if len(polylines) >= 2:
        # Distance between the longest two polylines, measured at
        # matched span coordinates.  Use min-distance per p0-point to
        # the closest p1-point.
        p1 = polylines[1]
        dists = np.linalg.norm(p0[:, None, :] - p1[None, :, :], axis=-1)  # (n0, n1)
        avg_min = float(dists.min(axis=1).mean())
        body_extent = float(np.linalg.norm(mesh.bounds[1] - mesh.bounds[0]))
        if avg_min < 0.05 * body_extent:
            is_blunt = True
            te_thickness = avg_min

    if n_groups > 1 and len(selected_groups) == 1:
        warnings.warn(
            f"detect_te_edges: found {n_groups} disconnected span-aligned "
            f"groups; using the longest ({len(p0)} vertices).  Other "
            "groups had span extent below 80% of the longest and were "
            "discarded.  Inspect the STL if this seems wrong.",
            UserWarning,
            stacklevel=2,
        )

    return TEDetectionResult(
        polylines=polylines,
        is_blunt=is_blunt,
        te_thickness=te_thickness,
        span_axis=span_axis,
        max_dihedral_deg=max_dihedral,
        n_candidate_groups=n_groups,
    )


# ---------------------------------------------------------------------------
# Span slicing
# ---------------------------------------------------------------------------


@dataclass
class SpanSliceResult:
    """Per-station cross-section polyline plus its TE anchor.

    Attributes
    ----------
    polyline_2d :
        ``(ni_body, 2)`` arclength-resampled cross-section polyline in
        the local in-plane coordinate frame (chord axis x thickness
        axis).  Ordered upper-TE -> LE -> lower-TE.
    te_position_2d :
        ``(2,)`` location of the upper-TE anchor on this station, in
        the same local 2D frame.  For sharp TE this coincides with
        ``polyline_2d[0]`` and ``polyline_2d[-1]``; for blunt TE it
        sits between them.
    span_coord :
        Span-axis coordinate of this station.
    """

    polyline_2d: np.ndarray
    te_position_2d: np.ndarray
    span_coord: float


def slice_by_span(
    mesh: trimesh.Trimesh,
    *,
    n_stations: int,
    ni_body: int,
    span_axis: str = "z",
    te_result: TEDetectionResult | None = None,
) -> list[SpanSliceResult]:
    """Slice the STL at ``n_stations`` planes perpendicular to span.

    Each slice produces one closed cross-section polyline.  The
    polyline is rotated into a local 2D frame (the span axis dropped),
    re-ordered so it starts at the upper TE and walks around the LE,
    and resampled by **arclength fraction** to a common ``ni_body``
    node count.  This common parametrisation lets stations with
    different chord lengths (tapered / swept wings) line up at the
    same parametric i-position — i = 0 is upper-TE on every station,
    i = ni_body//2 is approximately the LE on every station,
    i = ni_body-1 is lower-TE.

    Parameters
    ----------
    mesh :
        Watertight body STL.
    n_stations :
        Number of span stations to slice.  Stations are placed
        uniformly between the body's span-axis min and max, with the
        first and last stations strictly inside the body (offset by
        2% of the span extent) to avoid degenerate slices on the
        end-cap faces.
    ni_body :
        Common point count of the resampled cross-section polyline.
        Must be ``>= 8``.  Even values give symmetric LE/TE
        placement; odd values place a single node exactly at the LE.
    span_axis :
        ``"x"``, ``"y"``, or ``"z"``.  Default ``"z"``.
    te_result :
        Optional ``TEDetectionResult`` from
        :func:`detect_te_edges`.  When supplied, the per-station
        polyline ordering anchors at the TE; without it, the
        function picks the polyline endpoint closest to the
        chord-axis maximum as the TE — which works for symmetric
        airfoils but is fragile for cambered sections.

    Returns
    -------
    list[SpanSliceResult]
        One entry per station, sorted by span coordinate.

    Raises
    ------
    MeshInputError
        If a slice produces no polyline (station outside the body),
        a non-closed polyline (body has a hole at this station), or
        multiple disconnected polylines (body is multi-component).
    ParameterError
        On invalid arguments.
    """
    if span_axis not in _AXIS_INDEX:
        raise ParameterError(f"slice_by_span: span_axis must be 'x' / 'y' / 'z', got {span_axis!r}")
    if n_stations < 2:
        raise ParameterError(f"slice_by_span: n_stations must be >= 2, got {n_stations}")
    if ni_body < 8:
        raise ParameterError(f"slice_by_span: ni_body must be >= 8, got {ni_body}")

    span_idx = _AXIS_INDEX[span_axis]
    non_span_axes = [a for a in (0, 1, 2) if a != span_idx]
    chord_axis = (
        non_span_axes[0]
        if (mesh.bounds[1, non_span_axes[0]] - mesh.bounds[0, non_span_axes[0]])
        >= (mesh.bounds[1, non_span_axes[1]] - mesh.bounds[0, non_span_axes[1]])
        else non_span_axes[1]
    )
    thick_axis = non_span_axes[1] if chord_axis == non_span_axes[0] else non_span_axes[0]

    span_min = float(mesh.bounds[0, span_idx])
    span_max = float(mesh.bounds[1, span_idx])
    span_extent = span_max - span_min
    # Place stations strictly inside the body (offset by 2% of the span
    # extent at each end) so the slicer doesn't intersect the end-cap
    # triangles at glancing angles.
    margin = 0.02 * span_extent
    span_coords = np.linspace(span_min + margin, span_max - margin, n_stations)

    # Per-station TE position lookup: linear-interpolate the
    # te_result's polyline (sorted by span_axis) to each station.
    # For blunt TE, also build a lookup on the second strand (lower TE)
    # so the body wrap can terminate cleanly at the lower-TE node
    # without traversing the TE-base segment.
    is_blunt = te_result is not None and te_result.is_blunt
    if te_result is not None:
        te_polyline = te_result.polyline  # (n_te, 3) — upper strand
        te_span = te_polyline[:, span_idx]
        te_chord = te_polyline[:, chord_axis]
        te_thick = te_polyline[:, thick_axis]
        if is_blunt:
            lower_polyline = te_result.polylines[1]  # (n_te, 3) — lower strand
            lower_span = lower_polyline[:, span_idx]
            lower_chord = lower_polyline[:, chord_axis]
            lower_thick = lower_polyline[:, thick_axis]
        else:
            lower_span = lower_chord = lower_thick = None  # type: ignore[assignment]
    else:
        te_span = te_chord = te_thick = None  # type: ignore[assignment]
        lower_span = lower_chord = lower_thick = None  # type: ignore[assignment]

    plane_normal = np.zeros(3)
    plane_normal[span_idx] = 1.0

    results: list[SpanSliceResult] = []
    for s_coord in span_coords:
        plane_origin = np.zeros(3)
        plane_origin[span_idx] = s_coord
        # mesh_plane returns a list of segments (line strips).  We
        # then stitch them into a closed polyline.
        try:
            section_2d_lines = trimesh.intersections.mesh_plane(
                mesh, plane_normal=plane_normal, plane_origin=plane_origin
            )
        except Exception as exc:
            raise MeshInputError(
                f"slice_by_span: trimesh.mesh_plane failed at span={s_coord:.4f}: {exc}"
            ) from exc

        if len(section_2d_lines) == 0:
            raise MeshInputError(
                f"slice_by_span: no intersection at span={s_coord:.4f} — "
                "station is outside the body's span extent."
            )

        # mesh_plane gives a (m, 2, 3) array of line segments in 3D.
        # Drop the span axis to get the 2D polyline.
        segs_2d = section_2d_lines[:, :, [chord_axis, thick_axis]]  # (m, 2, 2)
        polyline_2d = _stitch_segments(segs_2d)
        if polyline_2d is None:
            raise MeshInputError(
                f"slice_by_span: cross-section at span={s_coord:.4f} "
                "is non-closed or multi-component (body has a hole or is "
                "multi-component at this station)."
            )

        # TE anchor(s) for this station: interpolate from te_result.
        lower_anchor: np.ndarray | None = None
        if te_result is not None:
            te_chord_at_s = float(np.interp(s_coord, te_span, te_chord))
            te_thick_at_s = float(np.interp(s_coord, te_span, te_thick))
            te_anchor = np.array([te_chord_at_s, te_thick_at_s])
            if is_blunt:
                lower_chord_at_s = float(np.interp(s_coord, lower_span, lower_chord))
                lower_thick_at_s = float(np.interp(s_coord, lower_span, lower_thick))
                lower_anchor = np.array([lower_chord_at_s, lower_thick_at_s])
        else:
            # Fallback: the TE is the point on the polyline at maximum chord.
            te_idx = int(np.argmax(polyline_2d[:, 0]))
            te_anchor = polyline_2d[te_idx]

        # Reorder the polyline so it starts at the upper-TE and walks
        # forward (around the LE) to the lower-TE.  For blunt TE, the
        # walk terminates at the lower-TE anchor (excluding the TE-base
        # segment from the body wrap).
        ordered = _reorder_polyline_around_te(polyline_2d, te_anchor, lower_anchor)

        # Resample by arclength fraction to ``ni_body`` points.
        resampled = _resample_by_arclength(ordered, ni_body)

        results.append(
            SpanSliceResult(
                polyline_2d=resampled,
                te_position_2d=te_anchor,
                span_coord=float(s_coord),
            )
        )

    return results


def _stitch_segments(segs_2d: np.ndarray) -> np.ndarray | None:
    """Stitch (m, 2, 2) line segments into a single closed polyline.

    Returns ``None`` if the segments don't form a single closed loop
    (multiple components, or an open curve).
    """
    if len(segs_2d) == 0:
        return None
    # Build a vertex set with merged near-duplicates.
    flat = segs_2d.reshape(-1, 2)
    # Tolerance for vertex matching: a tiny fraction of the local
    # extent of segs.  Using 1e-6 absolute is safe for normalised
    # geometries; for body coordinates we scale by the bounding diag.
    tol = max(1e-9, 1e-6 * float(np.linalg.norm(flat.max(axis=0) - flat.min(axis=0))))

    # Round vertices for matching.  Trimesh segments often share
    # endpoints with bit-exact equality already, but the safety net.
    vid_map: dict[tuple[int, int], int] = {}
    unique_verts: list[np.ndarray] = []
    seg_endpoints: list[tuple[int, int]] = []  # (vid_start, vid_end) per seg
    grid = np.maximum(tol, 1e-12)
    for s in segs_2d:
        endpts = []
        for p in s:
            key = (round(p[0] / grid), round(p[1] / grid))
            vid = vid_map.get(key)
            if vid is None:
                vid = len(unique_verts)
                vid_map[key] = vid
                unique_verts.append(p.copy())
            endpts.append(vid)
        seg_endpoints.append(tuple(endpts))  # type: ignore[arg-type]

    # Dedupe colinear duplicate segments.  ``trimesh.mesh_plane`` may
    # emit two cuts on each side of a triangle-pair diagonal when the
    # plane intersects the shared diagonal exactly — both cuts share
    # the diagonal vertex.  A "duplicate" is two segments with the
    # same unordered vertex pair.  Keep only one.
    seen_pairs: dict[tuple[int, int], int] = {}
    deduped_endpoints: list[tuple[int, int]] = []
    for a, b in seg_endpoints:
        if a == b:
            continue  # zero-length segment (vertex collapsed under rounding)
        key = (a, b) if a < b else (b, a)
        if key in seen_pairs:
            continue
        seen_pairs[key] = len(deduped_endpoints)
        deduped_endpoints.append((a, b))
    seg_endpoints = deduped_endpoints
    if len(seg_endpoints) == 0:
        return None

    # Build adjacency: each vertex should appear in exactly 2 segments
    # for a single closed loop.
    adj: dict[int, list[int]] = {}  # vid → list of segment indices
    for si, (a, b) in enumerate(seg_endpoints):
        adj.setdefault(a, []).append(si)
        adj.setdefault(b, []).append(si)
    if any(len(s) != 2 for s in adj.values()):
        return None  # open curve, T-junction, or duplicated edges

    # Walk the loop starting from the first segment.
    loop_vids: list[int] = []
    visited_segs: set[int] = set()
    cur_seg = 0
    cur_v, next_v = seg_endpoints[0]
    loop_vids.append(cur_v)
    while True:
        loop_vids.append(next_v)
        visited_segs.add(cur_seg)
        # Find the next segment incident to next_v that isn't cur_seg.
        cands = [s for s in adj[next_v] if s != cur_seg]
        if len(cands) != 1:
            return None  # malformed
        cur_seg = cands[0]
        a, b = seg_endpoints[cur_seg]
        cur_v, next_v = (next_v, b) if a == next_v else (next_v, a)
        if cur_seg in visited_segs:
            break
        if len(visited_segs) > len(seg_endpoints) + 5:
            return None  # loop bug guard

    if len(visited_segs) != len(seg_endpoints):
        return None  # multiple components

    # The walk closed back to the start; loop_vids has length n_seg + 1
    # with first == last.  Drop the duplicate.
    if loop_vids[0] != loop_vids[-1]:
        return None
    loop_vids.pop()

    polyline = np.array([unique_verts[v] for v in loop_vids])
    return polyline


def _reorder_polyline_around_te(
    polyline: np.ndarray,
    te_anchor: np.ndarray,
    lower_te_anchor: np.ndarray | None = None,
) -> np.ndarray:
    """Reorder a closed cross-section polyline to start at the upper-TE.

    The input polyline is closed (loop) but unordered relative to the
    TE.  Find the polyline node closest to ``te_anchor`` (= upper-TE),
    rotate the array so that node is at index 0, and orient the walk
    so that the next step climbs onto the suction (positive-thickness)
    surface.

    The rotation/flip ensures the canonical SHORE convention:
    ``i = 0`` is upper-TE, sweep over the suction surface to the LE,
    then back over the pressure surface to lower-TE at ``i = ni_body - 1``.

    Parameters
    ----------
    polyline :
        ``(N, 2)`` closed cross-section polyline (LE / TE order arbitrary).
    te_anchor :
        ``(2,)`` upper-TE position; the walk starts at the polyline
        node closest to this point.
    lower_te_anchor :
        Optional ``(2,)`` lower-TE position.  When supplied (blunt TE),
        the walk **stops** at the polyline node closest to this point,
        excluding the TE-base segment from the body wrap.  When omitted
        (sharp TE), the walk traverses the entire closed polyline and
        re-appends the start vertex so arclength parametrisation
        includes the closing edge — both endpoints land at the same
        point.
    """
    # Closest node to te_anchor.
    dists = np.linalg.norm(polyline - te_anchor[None, :], axis=1)
    te_idx = int(np.argmin(dists))

    # Rotate so te_idx is first.
    rotated = np.roll(polyline, -te_idx, axis=0)

    # Determine winding: walk one step forward.  If thickness becomes
    # positive (we move up onto the suction surface), the winding is
    # already correct.  Otherwise flip.
    next_thick = rotated[1, 1]
    cur_thick = rotated[0, 1]
    if next_thick - cur_thick < 0:
        # Walking down — flip the loop so we walk up first.
        rotated = np.concatenate([[rotated[0]], rotated[:0:-1]])

    if lower_te_anchor is None:
        # Sharp TE: closed polyline, re-append the start vertex at
        # the end so arclength parametrisation includes the closing
        # edge (both endpoints coincide in space).
        return np.vstack([rotated, rotated[0:1]])

    # Blunt TE: truncate the walk at the polyline node closest to
    # the lower-TE anchor.  This excludes the TE-base segment of the
    # closed polyline from the body wrap — the wake-base block W
    # meshes the TE base in the 4-block topology.
    lower_dists = np.linalg.norm(rotated - lower_te_anchor[None, :], axis=1)
    lower_idx = int(np.argmin(lower_dists))
    if lower_idx == 0:
        raise MeshInputError(
            "_reorder_polyline_around_te: lower-TE anchor closest to "
            "the upper-TE node — degenerate blunt TE detection."
        )
    return rotated[: lower_idx + 1]


# ---------------------------------------------------------------------------
# C-curve construction (body + wake branches)
# ---------------------------------------------------------------------------


@dataclass
class CCurveResult:
    """Per-station C-curve: body wrap with upper + lower wake branches.

    Attributes
    ----------
    polyline_2d :
        ``(ni_total, 2)`` array in the local 2D frame
        (chord axis x thickness axis).  Layout, in i-order:

        - ``[0 : ni_wake_upper]`` — upper wake, far-field → upper-TE
        - ``[ni_wake_upper : ni_wake_upper + ni_body]`` — body, upper-TE → LE → lower-TE
        - ``[ni_wake_upper + ni_body : ni_total]`` — lower wake, lower-TE → far-field

        Equivalently: ``i = 0`` is the upper-wake far-field tip; ``i =
        ni_total - 1`` is the lower-wake far-field tip.  The TE sits
        at ``i = ni_wake_upper`` (upper) and ``i = ni_wake_upper +
        ni_body - 1`` (lower).
    ni_wake_upper, ni_wake_lower :
        Node counts of each wake branch (including the TE node, which
        is shared with the body).  When the C-curve is flattened, the
        shared TE node appears twice — once as the upper-wake terminus,
        once as the body i=0 — but the consumer typically uses
        ``polyline_2d`` directly, where each TE is a single node.
    span_coord :
        Span-axis coordinate of this station.
    """

    polyline_2d: np.ndarray
    ni_wake_upper: int
    ni_wake_lower: int
    span_coord: float


def build_c_curve(
    slice_result: SpanSliceResult,
    *,
    ni_wake: int,
    wake_length: float,
    freestream_2d: np.ndarray,
    wake_growth: float = 1.15,
) -> CCurveResult:
    """Add upper + lower wake branches to a station body polyline.

    The body polyline (output of :func:`slice_by_span`) is a wrap from
    upper-TE around the LE to lower-TE.  This function appends two
    wake branches in the global freestream direction:

    - **Upper branch**: starts at the upper-wake far-field tip,
      contracts toward the upper-TE; concatenated **before** the
      body wrap.
    - **Lower branch**: starts at the lower-TE, expands toward the
      lower-wake far-field tip; concatenated **after** the body wrap.

    Wake i-spacing is **geometric**: starting from the body's TE-cell
    edge length and growing by ``wake_growth`` per cell, until the
    wake reaches ``wake_length`` (truncated to ``ni_wake`` nodes).

    Parameters
    ----------
    slice_result :
        Output of :func:`slice_by_span` for one station — a body wrap
        polyline ordered upper-TE → LE → lower-TE.
    ni_wake :
        Number of nodes per wake branch.  Each branch has this many
        nodes, **including** the TE node it shares with the body.
        Must be ``>= 2``.  Total C-curve length is ``2 * (ni_wake - 1)
        + ni_body``.
    wake_length :
        Length of each wake branch from the TE, in body-coordinate
        units.  Typical value: 20 * chord.
    freestream_2d :
        Length-2 unit vector ``(u_chord, u_thick)`` giving the
        freestream direction in the station's local frame.  Should
        approximate ``[+1, 0]`` for unswept wings; for swept wings,
        the projection of the global freestream onto the station
        plane.
    wake_growth :
        Geometric growth ratio for wake i-spacing.  ``1.0`` gives
        uniform spacing; ``1.15`` (default) gives mild stretching.

    Returns
    -------
    CCurveResult
        ``polyline_2d`` is the full C-curve in the local 2D frame.

    Raises
    ------
    ParameterError
        On invalid arguments.
    """
    if ni_wake < 2:
        raise ParameterError(f"build_c_curve: ni_wake must be >= 2, got {ni_wake}")
    if wake_length <= 0.0:
        raise ParameterError(f"build_c_curve: wake_length must be > 0, got {wake_length}")
    if wake_growth <= 0.0:
        raise ParameterError(f"build_c_curve: wake_growth must be > 0, got {wake_growth}")
    fs = np.asarray(freestream_2d, dtype=float).reshape(-1)
    if fs.size != 2:
        raise ParameterError(f"build_c_curve: freestream_2d must have length 2, got {fs.size}")
    fs_norm = float(np.linalg.norm(fs))
    if fs_norm < 1e-12:
        raise ParameterError("build_c_curve: freestream_2d has zero magnitude")
    fs_unit = fs / fs_norm

    body = slice_result.polyline_2d  # (ni_body, 2)
    if len(body) < 4:
        raise ParameterError(f"build_c_curve: body polyline has too few points ({len(body)} < 4)")

    upper_te = body[0]
    lower_te = body[-1]

    # Body cell length at the TE — first segment for upper, last for lower.
    ds_upper = float(np.linalg.norm(body[1] - body[0]))
    ds_lower = float(np.linalg.norm(body[-1] - body[-2]))

    upper_wake = _build_wake_branch(
        anchor=upper_te,
        direction=fs_unit,
        ni_wake=ni_wake,
        wake_length=wake_length,
        ds_root=ds_upper,
        growth=wake_growth,
    )
    lower_wake = _build_wake_branch(
        anchor=lower_te,
        direction=fs_unit,
        ni_wake=ni_wake,
        wake_length=wake_length,
        ds_root=ds_lower,
        growth=wake_growth,
    )

    # Concatenate: upper-wake (far → TE), body (TE → LE → TE), lower-wake (TE → far).
    # Each wake branch has ni_wake nodes ending at the TE; we drop the
    # TE-coincident endpoint when concatenating to avoid duplication.
    # Upper branch comes in TE-first order from _build_wake_branch
    # (anchor first, far last); reverse it so far-field is i=0.
    upper_far_to_te = upper_wake[::-1]  # (ni_wake, 2): far → TE
    # The body's first node is upper-TE — already coincides with
    # upper_far_to_te[-1].  Drop the body's first node.
    body_no_upper_te = body[1:]
    # Lower branch is in TE-first order — its first node coincides
    # with lower_te (= body[-1]).  Drop the lower branch's first node.
    lower_te_to_far = lower_wake[1:]

    polyline = np.concatenate([upper_far_to_te, body_no_upper_te, lower_te_to_far], axis=0)

    return CCurveResult(
        polyline_2d=polyline,
        ni_wake_upper=ni_wake,
        ni_wake_lower=ni_wake,
        span_coord=slice_result.span_coord,
    )


def _build_wake_branch(
    *,
    anchor: np.ndarray,
    direction: np.ndarray,
    ni_wake: int,
    wake_length: float,
    ds_root: float,
    growth: float,
) -> np.ndarray:
    """Build a single wake branch starting at ``anchor``.

    Returns ``(ni_wake, 2)`` array.  ``out[0] == anchor``;
    ``out[-1]`` is at ``anchor + wake_length * direction``.  Spacing
    is geometric from ``ds_root`` at the root, but globally rescaled
    so the cumulative arclength sums to exactly ``wake_length``.
    """
    if growth == 1.0:
        # Uniform spacing.
        ts = np.linspace(0.0, wake_length, ni_wake)
    else:
        # Geometric: cell sizes ds_root * growth**k, k = 0..ni_wake-2.
        # Sum = ds_root * (growth**(ni_wake-1) - 1) / (growth - 1).
        # Rescale so total = wake_length.
        cells = ds_root * growth ** np.arange(ni_wake - 1)
        scale = wake_length / cells.sum()
        cells = cells * scale
        ts = np.concatenate([[0.0], np.cumsum(cells)])
    return anchor[None, :] + ts[:, None] * direction[None, :]


# ---------------------------------------------------------------------------
# k=0 orchestrator
# ---------------------------------------------------------------------------


@dataclass
class CGridK0Result:
    """3D k=0 surface layer of a C-grid topology.

    Attributes
    ----------
    nodes :
        ``(ni_total, n_stations, 3)`` array of world-frame node
        coordinates.  ``ni_total = 2 * (ni_wake - 1) + ni_body``.
    ni_body, ni_wake :
        Counts used to build the curve.  ``ni_body`` is the body wrap
        node count; ``ni_wake`` is the per-branch wake node count
        (each wake branch shares its TE node with the body — see
        :class:`CCurveResult` for the exact slice layout).
    span_axis :
        ``"x"``, ``"y"``, or ``"z"``.
    te_result :
        The :class:`TEDetectionResult` produced during k=0 build.
    """

    nodes: np.ndarray
    ni_body: int
    ni_wake: int
    span_axis: str
    te_result: TEDetectionResult

    @property
    def ni_total(self) -> int:
        """Total i-direction node count of the concatenated layer."""
        return 2 * (self.ni_wake - 1) + self.ni_body

    @property
    def n_stations(self) -> int:
        """Span-direction station count."""
        return int(self.nodes.shape[1])

    def per_block_layers(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return per-block k=0 layers for the 3-block C-grid topology.

        Splits the concatenated ``nodes`` array along i into three
        contiguous segments, one per :class:`CGridTopology` block:

        - **upper wake**: ``nodes[0 : ni_wake]`` — sweeps
          ``upper-wake-far-tip → upper-TE``.
        - **body**: ``nodes[ni_wake - 1 : ni_wake - 1 + ni_body]`` —
          sweeps ``upper-TE → suction → LE → pressure → lower-TE``.
        - **lower wake**: ``nodes[ni_wake - 1 + ni_body - 1 : ni_total]`` —
          sweeps ``lower-TE → lower-wake-far-tip``.

        The slices overlap by one node at each TE column (upper wake
        last = body first, body last = lower wake first), so the
        returned arrays share TE columns by **numpy view aliasing**
        — modifying a TE node through one view updates the partner
        view automatically.  This is the bit-exact node coincidence
        the 3-block topology's SHARED seams rely on.

        Returns
        -------
        (upper, body, lower) : tuple of np.ndarray
            Each ``(ni_block, n_stations, 3)`` view into ``self.nodes``.
            Shapes: ``(ni_wake, n_stations, 3)``,
            ``(ni_body, n_stations, 3)``, ``(ni_wake, n_stations, 3)``.
        """
        i_te_upper = self.ni_wake - 1
        i_te_lower = i_te_upper + self.ni_body - 1
        upper = self.nodes[: i_te_upper + 1, :, :]
        body = self.nodes[i_te_upper : i_te_lower + 1, :, :]
        lower = self.nodes[i_te_lower:, :, :]
        return upper, body, lower


def build_cgrid_k0(
    mesh: trimesh.Trimesh,
    *,
    ni_body: int,
    ni_wake: int,
    n_stations: int,
    wake_length: float,
    span_axis: str = "z",
    freestream_world: np.ndarray | None = None,
    te_dihedral_deg: float = 30.0,
    wake_growth: float = 1.15,
) -> CGridK0Result:
    """Build the C-grid k=0 surface layer from a body STL.

    Pipeline:

    1. :func:`detect_te_edges` — find the TE polyline.
    2. :func:`slice_by_span` — slice the STL into ``n_stations``
       cross-section polylines, each resampled to ``ni_body`` nodes
       and ordered upper-TE → LE → lower-TE.
    3. :func:`build_c_curve` — for each station, append upper +
       lower wake branches in the (per-station projected) freestream
       direction.
    4. Lift the resulting 2D C-curves back into 3D by re-inserting
       the span coordinate, producing the world-frame
       ``(ni_total, n_stations, 3)`` k=0 layer.

    Parameters
    ----------
    mesh :
        Watertight body STL.
    ni_body :
        Per-station body wrap node count (upper-TE → LE → lower-TE).
    ni_wake :
        Per-branch wake node count (shared TE node included).
    n_stations :
        Number of span-axis stations.  ``>= 2``.
    wake_length :
        Length of each wake branch from the TE, in body-coordinate
        units.
    span_axis :
        ``"x"`` / ``"y"`` / ``"z"``.  Default ``"z"``.
    freestream_world :
        Length-3 vector giving the freestream direction in world
        coordinates.  If ``None``, defaults to the unit vector along
        the body's dominant chord axis (the largest non-span body
        extent), pointing **toward** the TE.
    te_dihedral_deg :
        Dihedral threshold for TE detection.  Default ``30°``.
    wake_growth :
        Geometric growth ratio for wake i-spacing.

    Returns
    -------
    CGridK0Result
        ``nodes`` is the ``(ni_total, n_stations, 3)`` world-frame
        array, ready to be marched outward by a :class:`HexBlock`.

    Raises
    ------
    MeshInputError, ParameterError
        Propagated from the pipeline stages.
    """
    if span_axis not in _AXIS_INDEX:
        raise ParameterError(
            f"build_cgrid_k0: span_axis must be 'x' / 'y' / 'z', got {span_axis!r}"
        )
    span_idx = _AXIS_INDEX[span_axis]
    non_span_axes = [a for a in (0, 1, 2) if a != span_idx]
    chord_axis = (
        non_span_axes[0]
        if (mesh.bounds[1, non_span_axes[0]] - mesh.bounds[0, non_span_axes[0]])
        >= (mesh.bounds[1, non_span_axes[1]] - mesh.bounds[0, non_span_axes[1]])
        else non_span_axes[1]
    )
    thick_axis = non_span_axes[1] if chord_axis == non_span_axes[0] else non_span_axes[0]

    # Default freestream: along the chord axis, pointing toward the TE
    # (i.e. toward maximum chord coordinate).  This is the natural
    # convention for airfoils: the body's nose is upstream, the TE is
    # downstream, freestream blows from nose to tail.
    if freestream_world is None:
        fs = np.zeros(3)
        fs[chord_axis] = 1.0
    else:
        fs = np.asarray(freestream_world, dtype=float).reshape(-1)
        if fs.size != 3:
            raise ParameterError(
                f"build_cgrid_k0: freestream_world must have length 3, got {fs.size}"
            )
        if np.linalg.norm(fs) < 1e-12:
            raise ParameterError("build_cgrid_k0: freestream_world has zero magnitude")

    # Project the freestream into the station's local 2D frame
    # (chord, thickness).  The span component is dropped — wakes
    # extend in-plane, perpendicular to span.
    freestream_2d = np.array([fs[chord_axis], fs[thick_axis]])
    if np.linalg.norm(freestream_2d) < 1e-12:
        raise ParameterError(
            "build_cgrid_k0: freestream_world is parallel to the span axis — "
            "wake direction is undefined.  Supply a freestream with a "
            "non-zero in-plane component."
        )

    te_result = detect_te_edges(mesh, dihedral_deg=te_dihedral_deg, span_axis=span_axis)
    slices = slice_by_span(
        mesh,
        n_stations=n_stations,
        ni_body=ni_body,
        span_axis=span_axis,
        te_result=te_result,
    )
    if len(slices) != n_stations:
        raise MeshInputError(f"build_cgrid_k0: expected {n_stations} slices, got {len(slices)}")

    # Build a C-curve per station.
    c_curves = [
        build_c_curve(
            s,
            ni_wake=ni_wake,
            wake_length=wake_length,
            freestream_2d=freestream_2d,
            wake_growth=wake_growth,
        )
        for s in slices
    ]

    # Sanity: every C-curve must have the same ni_total (otherwise the
    # k=0 layer can't be a regular structured grid).
    ni_total = c_curves[0].polyline_2d.shape[0]
    for k, cc in enumerate(c_curves):
        if cc.polyline_2d.shape[0] != ni_total:
            raise MeshInputError(
                f"build_cgrid_k0: C-curve at station {k} has "
                f"{cc.polyline_2d.shape[0]} nodes, expected {ni_total}.  "
                "All stations must share i-resolution."
            )

    # Lift the 2D C-curves back to 3D.  The local 2D frame uses
    # (chord_axis, thick_axis); we re-insert the span_axis coordinate
    # at each station.
    nodes = np.empty((ni_total, n_stations, 3), dtype=float)
    for k, cc in enumerate(c_curves):
        nodes[:, k, chord_axis] = cc.polyline_2d[:, 0]
        nodes[:, k, thick_axis] = cc.polyline_2d[:, 1]
        nodes[:, k, span_idx] = cc.span_coord

    return CGridK0Result(
        nodes=nodes,
        ni_body=ni_body,
        ni_wake=ni_wake,
        span_axis=span_axis,
        te_result=te_result,
    )


def _resample_by_arclength(polyline: np.ndarray, n: int) -> np.ndarray:
    """Resample a polyline at *n* equally-spaced arclength fractions.

    Parameters
    ----------
    polyline : (m, 2)
        Open polyline (treated as a sequence of points).
    n :
        Output point count.

    Returns
    -------
    (n, 2) array
        Resampled polyline.  ``out[0] == polyline[0]``, ``out[-1] ==
        polyline[-1]``.
    """
    if len(polyline) < 2:
        raise ValueError("_resample_by_arclength: polyline has fewer than 2 points")
    seg_lens = np.linalg.norm(np.diff(polyline, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_lens)])
    total = cum[-1]
    if total <= 0:
        raise ValueError("_resample_by_arclength: polyline has zero total length")
    targets = np.linspace(0.0, total, n)
    # Per-target segment index via searchsorted.
    out = np.empty((n, polyline.shape[1]))
    for i, t in enumerate(targets):
        idx = int(np.searchsorted(cum, t, side="right") - 1)
        idx = max(0, min(idx, len(seg_lens) - 1))
        seg_t = (t - cum[idx]) / max(seg_lens[idx], 1e-15)
        out[i] = (1.0 - seg_t) * polyline[idx] + seg_t * polyline[idx + 1]
    return out


def round_body_te_in_surface(
    surface: np.ndarray,
    ni_wake: int,
    ni_body: int,
    te_thickness: float,
    span_axis: str,
) -> None:
    """Round the body's TE in-place over the full span of a C-grid k=0 surface.

    Splits the body's i=0 (upper-TE) and i=ni_body-1 (lower-TE) columns
    apart by ``te_thickness`` along the airfoil-thickness direction at
    every spanwise station.  Also shifts each wake block's k=0 strip
    in bulk so the upper wake sits at ``+te_thickness/2`` and the lower
    wake at ``-te_thickness/2`` along their full chord-axis extent —
    not just at the TE-side endpoints.

    Without the wake-strip shift the wake blocks' k=0 strips converge
    back to ``y=0`` immediately downstream of the TE (the wake builder
    in :func:`build_cgrid_k0` runs each branch from a sharp anchor at
    y=0), producing a degenerate near-zero-thickness ``wake_base``
    block in the 4-block C-grid layout — a spurious singular line just
    downstream of the TE.

    The body sub-region of the concatenated k=0 array is
    ``surface[ni_wake-1 : ni_wake-1 + ni_body, :, :]``.  The upper-wake
    region is ``surface[: ni_wake, :, :]`` (i=0 far-tip → i=ni_wake-1
    upper-TE; the latter aliases the body's i=0 column).  The lower-
    wake region is ``surface[ni_wake-1 + ni_body - 1 :, :, :]``
    (i=0 lower-TE → i=ni_wake-1 far-tip; the former aliases the body's
    i=ni_body-1 column).

    Parameters
    ----------
    surface :
        ``(ni_total, n_stations, 3)`` k=0 surface array as produced by
        :func:`build_cgrid_k0`.  Must have a sharp TE at every station.
    ni_wake, ni_body :
        Block-count integers from the C-grid construction.
    te_thickness :
        World-frame TE thickness to apply.  Must be > 0.
    span_axis :
        ``'x'`` / ``'y'`` / ``'z'``.

    Raises
    ------
    ParameterError
        On invalid arguments or non-sharp input.
    """
    if te_thickness <= 0:
        raise ParameterError(
            f"round_body_te_in_surface: te_thickness must be > 0, got {te_thickness}"
        )
    if span_axis not in ("x", "y", "z"):
        raise ParameterError(
            f"round_body_te_in_surface: span_axis must be 'x'/'y'/'z', got {span_axis!r}"
        )
    span_axis_world = "xyz".index(span_axis)
    n_stations = surface.shape[1]
    i_te_upper = ni_wake - 1
    i_te_lower = i_te_upper + ni_body - 1

    # Sanity: the body must be sharp at every station before rounding.
    for j in range(n_stations):
        delta = surface[i_te_upper, j, :] - surface[i_te_lower, j, :]
        if float(np.linalg.norm(delta)) > 1e-9:
            raise ParameterError(
                f"round_body_te_in_surface: station j={j} already has a "
                f"non-sharp TE (separation {float(np.linalg.norm(delta)):.6g}); "
                f"this helper only rounds bodies with a sharp TE."
            )

    # Per-station: identify the chord-axis (largest in-plane extent of the
    # body's i-axis) and the thickness-axis (the remaining in-plane axis).
    # Apply the anti-symmetric split to the body's TE columns AND to the
    # full wake blocks' k=0 strips so the wake stays at y = ±te/2 over
    # its full length (otherwise the wake blocks' k=0 strips collapse
    # back to y=0 just downstream of the TE — see docstring).
    half = 0.5 * te_thickness
    for j in range(n_stations):
        body_xs = surface[i_te_upper : i_te_lower + 1, j, :]
        extents = body_xs.max(axis=0) - body_xs.min(axis=0)
        extents[span_axis_world] = -np.inf
        chord_axis = int(np.argmax(extents))
        thickness_axis = next(
            a for a in (0, 1, 2) if a != span_axis_world and a != chord_axis
        )
        # 1. Body TE columns: anti-symmetric split.  These columns are
        #    aliased to the wake blocks' TE-side endpoints, so shifting
        #    them moves those endpoints too.
        surface[i_te_upper, j, thickness_axis] += +half
        surface[i_te_lower, j, thickness_axis] += -half
        # 2. Upper wake's k=0 strip (i = 0 .. ni_wake - 2): shift
        #    +half along thickness_axis so the whole strip sits at
        #    y = +te/2.  Index ni_wake - 1 was already shifted in step 1.
        surface[0 : ni_wake - 1, j, thickness_axis] += +half
        # 3. Lower wake's k=0 strip (i = i_te_lower + 1 .. ni_total - 1):
        #    shift -half so the whole strip sits at y = -te/2.  Index
        #    i_te_lower was already shifted in step 1.
        surface[i_te_lower + 1 :, j, thickness_axis] += -half
