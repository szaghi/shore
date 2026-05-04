"""Cap-fill builder for the C-grid topology's wing-tip closures.

The C-grid topology wraps a wing's outer surface but leaves the spanwise
boundaries (j=0 and j=n_stations-1) FREE — the wing is implicitly assumed
to extend to infinity in span, or to terminate at a Chimera-fringe boundary.
For a finite wing with closed tips, the j_lo / j_hi faces of the body wrap
need a **wall** closure: a thin spanwise slab covering the airfoil interior
cross-section at the tip plane.

This module builds that closure as **3 cap sub-blocks** (sharp TE) or **4
cap sub-blocks** (blunt TE) per tip, in an H-style "T-fill" layout that has
no singular axis (Xall-compatible).

Cap-block axis convention
-------------------------
Each cap sub-block is shaped ``(ni_cap, nj_cap, 2, 3)``:

- **i-axis**: along the airfoil profile (chordwise-ish).  For ``cap_upper``
  i runs from upper-LE-cut (i=0) to upper-TE-cut (i=ni_cap-1) — *reversed*
  w.r.t. the body's i-axis.  For ``cap_lower`` i runs from lower-LE-cut to
  lower-TE-cut — forward w.r.t. body i.  For ``cap_le`` i wraps from
  upper-LE-cut around the LE point to lower-LE-cut.
- **j-axis**: airfoil-thickness direction (perpendicular to chord, in the
  cross-section's 2D plane).  For ``cap_upper`` j runs from chord-cut at
  y=0 (j=0) to airfoil suction wall (j=nj_cap-1).  ``cap_lower`` is the
  mirror.  ``cap_le`` j runs from chord-perp inner line (j=0) to airfoil
  profile wrapping the LE (j=nj_cap-1).
- **k-axis**: spanwise.  **2 nodes, 1 cell** thick.  k=0 is the wing-tip
  wall (offset by ``tip_inset`` along ``span_dir`` from the body cross-section),
  k=1 is the body cross-section plane (SHARED with the body wrap's j_lo).

The resulting block is **not marched**.  It is a pure parametric fill,
fully determined by the body's j=0 cross-section polyline plus
``tip_inset`` and the cap j-spacing law.
"""

from __future__ import annotations

import numpy as np

from shore.exceptions import ParameterError
from shore.volume.spacing import build_i_coords


def build_cap_blocks(
    cross_section: np.ndarray,
    le_cut_upper_idx: int,
    le_cut_lower_idx: int,
    span_dir: np.ndarray,
    tip_inset: float,
    n_cap_j: int,
    nk_cap: int = 2,
    j_spacing: str = "tanh",
    j_beta: float = 4.0,
    is_blunt: bool = False,
    wake_base_cross_section: np.ndarray | None = None,
    te_cut_upper_idx: int | None = None,
    te_cut_lower_idx: int | None = None,
) -> dict[str, np.ndarray]:
    """Build the cap sub-blocks at one wing tip of a C-grid.

    K-axis convention
    -----------------
    The cap is a 3D structured slab spanning the wing-tip region.
    Its k-axis is **spanwise**:

    - **k=0**:        outboard wing-tip wall, offset from
                      ``cross_section`` by ``tip_inset * span_dir``.
                      Wired to ``WALL`` by the topology layer.
    - **k=nk_cap-1**: inboard chimera-overlap plane, coincident with
                      ``cross_section``.  Wired to ``FREE`` (chimera
                      fringe receiving from the C-grid body wrap) by
                      the topology layer.

    Interior layers (``1 <= k <= nk_cap - 2``) interpolate linearly
    between the two end faces along ``span_dir``.

    Parameters
    ----------
    cross_section :
        ``(ni_body_post_ext, 3)`` array.  The body's spanwise-extreme
        cross-section polyline.  i sweeps upper-TE-cut → suction →
        LE → pressure → lower-TE-cut.
    le_cut_upper_idx, le_cut_lower_idx :
        i-indices in ``cross_section`` at the LE-cut on the suction and
        pressure sides respectively.  Must satisfy
        ``0 < le_cut_upper_idx < le_cut_lower_idx < ni_body - 1``.
    te_cut_upper_idx, te_cut_lower_idx :
        Optional i-indices in ``cross_section`` at the TE-cut on each
        side.  When ``None`` (default), the cap chord-band blocks
        extend all the way to the body's i=0 / i=ni-1 endpoints —
        which **collapses to a singular axis at a sharp TE** (where
        the wall point coincides with the chord cut).  For sharp TE
        always pass valid TE-cut indices that sit *strictly inboard*
        of the TE point so the cap_upper.i_hi / cap_lower.i_hi
        columns are non-degenerate.  The fluid sliver between the
        TE-cut and the actual TE point is left unmeshed by the cap;
        chimera-fringe overlap with the C-grid wake blocks closes
        the gap.  Must satisfy
        ``0 <= te_cut_upper_idx < le_cut_upper_idx < le_cut_lower_idx
        < te_cut_lower_idx <= ni_body - 1`` and the chord-band
        symmetry constraint
        ``le_cut_upper_idx - te_cut_upper_idx == te_cut_lower_idx -
        le_cut_lower_idx`` (equal i-station counts for cap_upper and
        cap_lower so the chord-cut SHARED face is 1:1).
    span_dir :
        Unit vector pointing **outward** along the span axis at this tip
        — i.e. away from the body interior.  Used to offset the wing-tip
        wall plane from the body cross-section by ``tip_inset``.
    tip_inset :
        Spanwise thickness of the cap slab (k-axis extent in world units).
        Must be > 0.
    n_cap_j :
        Node count along the cap's j-axis (airfoil-thickness direction).
        Must be >= 2.  All cap sub-blocks share this resolution along
        their j-axis.
    nk_cap :
        Node count along the cap's k-axis (spanwise).  Must be >= 2.
        Default ``2`` keeps a single-cell-thick slab; larger values
        give finer spanwise resolution at the wing tip.
    j_spacing :
        Spacing law for the cap j-axis: ``'uniform'`` / ``'tanh'`` /
        ``'tanh2'``.  ``'tanh'`` (default) clusters nodes near the
        airfoil wall (j_hi side) for boundary-layer resolution.
    j_beta :
        Clustering strength for tanh / tanh2.
    is_blunt :
        When True, includes a 4th cap sub-block ``'wake_base'`` for the
        TE-base region.  ``wake_base_cross_section`` must be supplied.
    wake_base_cross_section :
        Required when ``is_blunt=True``.  Shape ``(ni_wake, nk_w, 3)``:
        the wake_base block's j=0 face nodes (the wake_base's spanwise
        end face, which the cap closes).  ``ni_wake`` is the wake's
        node count, ``nk_w`` is the wake_base's k-axis resolution
        across the TE thickness.

    Returns
    -------
    dict[str, np.ndarray]
        Maps cap sub-block label to its node array of shape
        ``(ni_cap, nj_cap, nk_cap, 3)``.  Keys are always
        ``'le'``, ``'upper'``, ``'lower'``, ``'te_upper'``, and
        ``'te_lower'`` (the 5-block H-fill that closes the airfoil
        interior without singular axes — provided the cross-section
        has a non-zero TE thickness).  When ``is_blunt=True`` *and*
        the user's body wrap has a real ``wake_base`` block, an
        extra ``'wake_base'`` cap sub-block is included.

    Raises
    ------
    ParameterError
        On invalid shapes, indices, or unsupported spacing law.
    """
    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    if not (
        isinstance(cross_section, np.ndarray)
        and cross_section.ndim == 2
        and cross_section.shape[-1] == 3
    ):
        raise ParameterError(
            f"build_cap_blocks: cross_section must have shape (ni, 3); "
            f"got {None if not isinstance(cross_section, np.ndarray) else cross_section.shape}"
        )
    ni_body = cross_section.shape[0]
    if not (0 < le_cut_upper_idx < le_cut_lower_idx < ni_body - 1):
        raise ParameterError(
            f"build_cap_blocks: le_cut_upper_idx ({le_cut_upper_idx}) must be "
            f"strictly less than le_cut_lower_idx ({le_cut_lower_idx}), and "
            f"both must be strictly inside [1, ni_body - 2] (ni_body={ni_body})."
        )
    span_arr = np.asarray(span_dir, dtype=np.float64)
    if span_arr.shape != (3,):
        raise ParameterError(
            f"build_cap_blocks: span_dir must have shape (3,); got {span_arr.shape}"
        )
    span_norm = float(np.linalg.norm(span_arr))
    if not np.isclose(span_norm, 1.0, atol=1e-9):
        raise ParameterError(
            f"build_cap_blocks: span_dir must be a unit vector; "
            f"got |span_dir|={span_norm:.6g}.  Normalise before calling."
        )
    if tip_inset <= 0:
        raise ParameterError(f"build_cap_blocks: tip_inset must be > 0, got {tip_inset}")
    if n_cap_j < 2:
        raise ParameterError(f"build_cap_blocks: n_cap_j must be >= 2, got {n_cap_j}")
    if nk_cap < 2:
        raise ParameterError(f"build_cap_blocks: nk_cap must be >= 2, got {nk_cap}")
    if j_spacing not in ("uniform", "tanh", "tanh2"):
        raise ParameterError(
            f"build_cap_blocks: j_spacing must be 'uniform' / 'tanh' / 'tanh2', got {j_spacing!r}"
        )
    if is_blunt and wake_base_cross_section is None:
        raise ParameterError(
            "build_cap_blocks: is_blunt=True requires wake_base_cross_section "
            "(the wake_base block's j=0 face nodes)."
        )

    # Default TE-cut indices to the body endpoints (preserves legacy
    # behaviour, keeps the singular-TE failure mode for sharp TE — but
    # callers should always pass real TE-cut indices in practice).
    if te_cut_upper_idx is None:
        te_cut_upper_idx = 0
    if te_cut_lower_idx is None:
        te_cut_lower_idx = ni_body - 1
    if not (
        0
        <= te_cut_upper_idx
        < le_cut_upper_idx
        < le_cut_lower_idx
        < te_cut_lower_idx
        <= ni_body - 1
    ):
        raise ParameterError(
            f"build_cap_blocks: TE-cut and LE-cut indices must satisfy "
            f"0 <= te_cut_upper ({te_cut_upper_idx}) < le_cut_upper "
            f"({le_cut_upper_idx}) < le_cut_lower ({le_cut_lower_idx}) "
            f"< te_cut_lower ({te_cut_lower_idx}) <= ni_body - 1 "
            f"({ni_body - 1})."
        )

    # ------------------------------------------------------------------
    # Build cap_upper, cap_lower, cap_le from the cross_section polyline.
    # ------------------------------------------------------------------
    # The k-axis offset: each cap node at k=0 is its k=1 counterpart minus
    # tip_inset * span_dir (k_hi sits on the body cross-section plane,
    # k_lo is the wing-tip wall offset outward).
    k_offset = tip_inset * span_arr  # (3,)

    # j parameter in [0, 1] under the chosen spacing law.  The cap j-axis
    # runs from chord-cut (j=0) to airfoil wall (j=nj-1); for boundary-
    # layer resolution we want clustering near the wall.  build_i_coords
    # clusters near its start (low end), so for 'tanh' / 'tanh2' we
    # invert the parameter so the cluster lands at the wall side.
    if j_spacing == "tanh":
        t_raw = build_i_coords(n_cap_j, 0.0, 1.0, spacing="tanh", beta=j_beta)
        t_j = 1.0 - t_raw[::-1]  # cluster near j=nj-1 (wall)
    else:
        # 'uniform' is symmetric; 'tanh2' is two-sided so already symmetric.
        t_j = build_i_coords(n_cap_j, 0.0, 1.0, spacing=j_spacing, beta=j_beta)

    blocks: dict[str, np.ndarray] = {}

    # cap_upper -------------------------------------------------------
    # i: 0 (upper-LE-cut) → ni_chord-1 (upper-TE-cut, body i = te_cut_upper_idx).
    # Wall sits at cross_section[le_cut_upper_idx - i] for i = 0..ni_chord-1.
    # j: 0 (chord cut at y=0 in airfoil-thickness frame) → nj-1 (wall).
    blocks["upper"] = _build_chordwise_cap(
        cross_section=cross_section,
        i_indices=np.arange(le_cut_upper_idx, te_cut_upper_idx - 1, -1, dtype=np.int64),
        cut_to_wall_t=t_j,
        k_offset=k_offset,
        nk_cap=nk_cap,
    )

    # cap_lower -------------------------------------------------------
    # i: 0 (lower-LE-cut) → ni_chord-1 (lower-TE-cut, body i = te_cut_lower_idx).
    # By symmetry of the H3 fill, cap_lower must have the same ni and nj
    # as cap_upper for the chord-cut SHARED to be 1:1.
    ni_chord_upper = le_cut_upper_idx - te_cut_upper_idx + 1
    ni_chord_lower = te_cut_lower_idx - le_cut_lower_idx + 1
    if ni_chord_upper != ni_chord_lower:
        raise ParameterError(
            f"build_cap_blocks: chord-cut SHARED requires symmetric LE/TE-cut "
            f"placement, but cap_upper has {ni_chord_upper} i-stations "
            f"while cap_lower has {ni_chord_lower}.  Both cuts must be "
            f"equidistant from the LE point in i-station count."
        )
    blocks["lower"] = _build_chordwise_cap(
        cross_section=cross_section,
        i_indices=np.arange(le_cut_lower_idx, te_cut_lower_idx + 1, dtype=np.int64),
        cut_to_wall_t=t_j,
        k_offset=k_offset,
        nk_cap=nk_cap,
    )

    # cap_le ----------------------------------------------------------
    # i: 0 (upper-LE-cut) → ni_le-1 (lower-LE-cut), wraps around LE.
    # j: 0 (chord-perp inner line) → nj-1 (airfoil profile wrapping LE).
    # The chord-perp inner line is the line perpendicular to chord at the
    # cap_le_frac chord position, spanning from the chord cut at y=0 up
    # to the corresponding wall point.  By construction it terminates on
    # the SHARED LE-cut columns of cap_upper and cap_lower.
    blocks["le"] = _build_le_cap(
        cross_section=cross_section,
        le_cut_upper_idx=le_cut_upper_idx,
        le_cut_lower_idx=le_cut_lower_idx,
        cut_to_wall_t=t_j,
        k_offset=k_offset,
        nk_cap=nk_cap,
    )

    # cap_te_upper / cap_te_lower -------------------------------------
    # The TE-region cap blocks fill the airfoil interior between the
    # TE-cut chord-perp line (chord = te_frac) and the (possibly
    # synthesised blunt) TE base at chord = 1.  They mirror cap_upper /
    # cap_lower's structure: i = chord-band, j = chord_cut → wall.
    #
    # cap_te_upper: i sweeps TE-cut (i=0, body i = te_cut_upper_idx) →
    #               TE-base-upper (i=ni-1, body i = 0).
    # cap_te_lower: i sweeps TE-cut (i=0, body i = te_cut_lower_idx) →
    #               TE-base-lower (i=ni-1, body i = ni_body-1).
    #
    # On a blunt-or-rounded airfoil cross-section the body's i=0 and
    # i=ni-1 nodes sit at distinct points (TE-base segment has finite
    # thickness), so cap_te_upper.i_hi and cap_te_lower.i_hi columns
    # are non-degenerate.  On a sharp cross-section i=0 and i=ni-1
    # coincide → cap_te_*.i_hi columns collapse and you get the
    # singular axis we are trying to avoid; ``CGridTipCap`` rejects
    # sharp inputs unless auto-rounding is enabled.
    ni_te_upper = te_cut_upper_idx + 1
    ni_te_lower = ni_body - te_cut_lower_idx
    if ni_te_upper != ni_te_lower:
        raise ParameterError(
            f"build_cap_blocks: TE-base SHARED requires symmetric TE-cut "
            f"placement, but cap_te_upper has {ni_te_upper} i-stations "
            f"while cap_te_lower has {ni_te_lower}.  TE cuts must be "
            f"equidistant from the body's TE endpoints in i-station count."
        )
    blocks["te_upper"] = _build_chordwise_cap(
        cross_section=cross_section,
        i_indices=np.arange(te_cut_upper_idx, -1, -1, dtype=np.int64),
        cut_to_wall_t=t_j,
        k_offset=k_offset,
        nk_cap=nk_cap,
    )
    blocks["te_lower"] = _build_chordwise_cap(
        cross_section=cross_section,
        i_indices=np.arange(te_cut_lower_idx, ni_body, dtype=np.int64),
        cut_to_wall_t=t_j,
        k_offset=k_offset,
        nk_cap=nk_cap,
    )

    # cap_wake_base (body-blunt only) ---------------------------------
    if is_blunt:
        wb = np.asarray(wake_base_cross_section, dtype=np.float64)
        if wb.ndim != 3 or wb.shape[-1] != 3:
            raise ParameterError(
                f"build_cap_blocks: wake_base_cross_section must have "
                f"shape (ni_wake, nk_w, 3); got {wb.shape}"
            )
        # The wake_base's spanwise face is already a 2D structured grid
        # (ni_wake, nk_w).  The cap is that grid extruded outward
        # along span_dir for nk_cap layers.  Axis convention matches
        # the chordwise cap blocks: ni_wake → i, nk_w → j, span → k.
        # k=0 is the outboard wing-tip wall, k=nk_cap-1 is the inboard
        # chimera-overlap plane (= wb).
        ni_wake, nk_w = wb.shape[:2]
        wake_base_blk = np.empty((ni_wake, nk_w, nk_cap, 3), dtype=np.float64)
        for k in range(nk_cap):
            t = k / (nk_cap - 1)  # 0 → outboard wall, 1 → cross-section plane
            wake_base_blk[:, :, k, :] = wb + (1.0 - t) * k_offset[None, None, :]
        blocks["wake_base"] = wake_base_blk

    return blocks


def _build_chordwise_cap(
    cross_section: np.ndarray,
    i_indices: np.ndarray,
    cut_to_wall_t: np.ndarray,
    k_offset: np.ndarray,
    nk_cap: int,
) -> np.ndarray:
    """Build cap_upper or cap_lower (the chord-band cap blocks).

    Each i corresponds to a body cross-section node (the wall point).
    The cap fills the strip from the chord-cut line (y=0 in airfoil
    frame, j=0 in block frame) to the wall (j=nj-1) along ``cut_to_wall_t``,
    extruded along ``span_dir`` for ``nk_cap`` layers between the
    outboard wing-tip wall (k=0) and the cross-section plane
    (k=nk_cap-1).

    The "chord-cut line" is the 2D y=0 line in the airfoil's local
    frame, recovered from the cross-section's centroid in the
    thickness direction.

    Parameters
    ----------
    cross_section :
        Body cross-section polyline.
    i_indices :
        Indices into ``cross_section`` for this cap block's i-axis.
        The wall-side nodes are ``cross_section[i_indices, :]``.
    cut_to_wall_t :
        ``(n_cap_j,)`` parameter in [0, 1] from chord cut (j=0) to wall (j=nj-1).
    k_offset :
        ``(3,)`` outward span offset (= tip_inset * span_dir).
    nk_cap :
        Spanwise node count.  k=0 is the outboard wall, k=nk_cap-1
        is the inboard chimera-overlap plane.

    Returns
    -------
    ndarray, shape ``(ni_cap, n_cap_j, nk_cap, 3)``.
    """
    n_cap_j = cut_to_wall_t.shape[0]
    ni_cap = i_indices.shape[0]

    # Wall points on the cross-section plane (k = nk_cap-1, j = nj-1).
    wall = cross_section[i_indices, :]  # (ni_cap, 3)

    # Chord-cut line: project the wall onto the centroid's thickness
    # coordinate (= the airfoil's y=0 line, in the typical case where
    # the airfoil cross-section is symmetric about y=0).
    centroid = cross_section.mean(axis=0)  # (3,)
    thickness_axis = _detect_thickness_axis(cross_section, k_offset)
    chord_cut = wall.copy()
    chord_cut[:, thickness_axis] = centroid[thickness_axis]

    # Build the in-plane (i, j) layer at the cross-section plane.
    plane = np.empty((ni_cap, n_cap_j, 3), dtype=np.float64)
    for j in range(n_cap_j):
        t = cut_to_wall_t[j]
        plane[:, j, :] = (1.0 - t) * chord_cut + t * wall

    # Extrude along +k_offset (= +tip_inset * span_dir, outward) for
    # nk_cap layers.  k=0 is the outboard wall (offset by +k_offset
    # from the plane); k=nk_cap-1 is the cross-section plane itself.
    blk = np.empty((ni_cap, n_cap_j, nk_cap, 3), dtype=np.float64)
    for k in range(nk_cap):
        s = k / (nk_cap - 1)  # 0 at outboard wall, 1 at cross-section plane
        blk[:, :, k, :] = plane + (1.0 - s) * k_offset[None, None, :]
    return blk


def _build_le_cap(
    cross_section: np.ndarray,
    le_cut_upper_idx: int,
    le_cut_lower_idx: int,
    cut_to_wall_t: np.ndarray,
    k_offset: np.ndarray,
    nk_cap: int,
) -> np.ndarray:
    """Build cap_le, the LE-region cap that wraps around the leading edge.

    i: 0 (upper-LE-cut) → ni_le-1 (lower-LE-cut), one node per body
    cross-section node along the LE wrap.
    j: 0 (chord-perp inner line) → nj-1 (airfoil profile wrap).
    k: 0 (outboard wall) → nk_cap-1 (cross-section plane).

    The chord-perp inner line is the vertical line at chord = cap_le_frac
    spanning from the upper-LE-cut wall point down to the lower-LE-cut
    wall point.  Its endpoints are SHARED with cap_upper.i_lo and
    cap_lower.i_lo — bit-exact node coincidence.
    """
    n_cap_j = cut_to_wall_t.shape[0]
    ni_le = le_cut_lower_idx - le_cut_upper_idx + 1

    wall = cross_section[le_cut_upper_idx : le_cut_lower_idx + 1, :]  # (ni_le, 3)

    thickness_axis = _detect_thickness_axis(cross_section, k_offset)
    centroid = cross_section.mean(axis=0)
    inner = wall.copy()
    inner[:, thickness_axis] = centroid[thickness_axis]

    plane = np.empty((ni_le, n_cap_j, 3), dtype=np.float64)
    for j in range(n_cap_j):
        t = cut_to_wall_t[j]
        plane[:, j, :] = (1.0 - t) * inner + t * wall

    blk = np.empty((ni_le, n_cap_j, nk_cap, 3), dtype=np.float64)
    for k in range(nk_cap):
        s = k / (nk_cap - 1)
        blk[:, :, k, :] = plane + (1.0 - s) * k_offset[None, None, :]
    return blk


def _detect_thickness_axis(
    cross_section: np.ndarray,
    k_offset: np.ndarray,
) -> int:
    """Identify which world axis (0, 1, 2) is the airfoil-thickness direction.

    The cross-section is a 2D curve embedded in 3D.  Two of the three
    world axes are in-plane (chord and thickness); the third is out-of-
    plane (span).  ``k_offset`` is the span direction.  We pick the
    in-plane axis with the smaller spread as the thickness axis (chord
    direction has the larger extent for a typical airfoil).
    """
    span_axis = int(np.argmax(np.abs(k_offset)))
    extents = cross_section.max(axis=0) - cross_section.min(axis=0)
    extents[span_axis] = np.inf  # exclude span axis from the argmin
    return int(np.argmin(extents))


def find_symmetric_le_cuts(
    cross_section: np.ndarray,
    le_frac: float,
    span_axis_world: int,
) -> tuple[int, int]:
    """Find the upper / lower LE-cut indices in a body cross-section.

    The cross-section sweeps upper-TE-cut → suction → LE → pressure →
    lower-TE-cut.  This locates the i-indices on each side of the LE
    point that correspond to chord position ``le_frac``, snapped to be
    **symmetric** around the LE point in i-station count (so
    ``upper_idx + lower_idx == 2 * le_idx``).  Symmetry is required by
    the cap-fill builder for the chord-cut SHARED face.

    Parameters
    ----------
    cross_section :
        ``(ni, 3)`` polyline.
    le_frac :
        Chord-fraction of the LE-cut, in ``(0, 1)``.
    span_axis_world :
        Which world axis (0 / 1 / 2) is the span direction; the
        chord axis is the in-plane axis with the largest extent.

    Returns
    -------
    (le_cut_upper_idx, le_cut_lower_idx) : tuple of int

    Raises
    ------
    ParameterError
        On invalid ``le_frac`` or when the i-axis is too coarse to
        resolve the requested cut.
    """
    if not (0.0 < le_frac < 1.0):
        raise ParameterError(f"find_symmetric_le_cuts: le_frac must be in (0, 1), got {le_frac}")
    if span_axis_world not in (0, 1, 2):
        raise ParameterError(
            f"find_symmetric_le_cuts: span_axis_world must be 0/1/2, got {span_axis_world}"
        )

    ni = cross_section.shape[0]
    extents = cross_section.max(axis=0) - cross_section.min(axis=0)
    extents[span_axis_world] = -np.inf
    chord_axis = int(np.argmax(extents))

    x = cross_section[:, chord_axis]
    x_le = float(x.min())
    x_te = float(x.max())
    chord_len = x_te - x_le
    le_idx = int(np.argmin(x))
    x_target = x_le + le_frac * chord_len

    upper_idx = le_idx
    while upper_idx > 0 and x[upper_idx] < x_target:
        upper_idx -= 1
    # Symmetry around the body's i-axis midpoint (which corresponds to
    # the LE for a typical airfoil i-distribution).  This guarantees
    # equal i-station counts for cap_upper and cap_lower regardless of
    # whether ``ni_body`` is even or odd.
    lower_idx = ni - 1 - upper_idx

    if lower_idx <= upper_idx:
        raise ParameterError(
            f"find_symmetric_le_cuts: chord fraction {le_frac} resolves to "
            f"a single i-station ({upper_idx}); pick a larger fraction or "
            f"use a finer i-axis."
        )
    if upper_idx == le_idx or lower_idx == le_idx:
        raise ParameterError(
            f"find_symmetric_le_cuts: chord fraction {le_frac} resolves to "
            f"the LE point itself; pick a larger fraction or use a finer i-axis."
        )
    return upper_idx, lower_idx


def detect_te_thickness(
    cross_section: np.ndarray,
    span_axis_world: int,
) -> float:
    """Return the TE thickness of a cross-section polyline.

    Sharp TE → ``0.0`` (i=0 and i=ni-1 coincide).  Blunt TE → distance
    between i=0 and i=ni-1 in the airfoil's thickness direction.
    """
    delta = cross_section[-1, :] - cross_section[0, :]
    in_plane = delta.copy()
    in_plane[span_axis_world] = 0.0
    return float(np.linalg.norm(in_plane))


def auto_round_sharp_te(
    cross_section: np.ndarray,
    te_thickness: float,
    span_axis_world: int,
) -> np.ndarray:
    """Synthesise a blunt cross-section by splitting the sharp TE node.

    Returns a copy of ``cross_section`` with the i=0 and i=ni-1 nodes
    moved apart by ``te_thickness`` along the airfoil-thickness
    direction (perpendicular to chord, in the airfoil's 2D plane), so
    the TE-base segment has finite extent.  Used by ``CGridTipCap`` to
    auto-round sharp-TE inputs *only inside the cap*; the user's
    body wrap is left untouched.

    The thickness direction is the in-plane axis with the smaller
    spread (the same heuristic used by ``_detect_thickness_axis``).
    The split is anti-symmetric: i=0 moves by +te_thickness/2 and
    i=ni-1 by -te_thickness/2 so the TE-base midpoint stays at the
    original sharp-TE location.

    Parameters
    ----------
    cross_section :
        ``(ni, 3)`` polyline.  Must have ``cross_section[0]`` and
        ``cross_section[-1]`` coincident (sharp TE).
    te_thickness :
        Desired TE thickness.  Must be > 0.
    span_axis_world :
        World axis (0 / 1 / 2) that is the spanwise direction; used
        to identify the airfoil-thickness in-plane axis.

    Returns
    -------
    ndarray, shape ``(ni, 3)``, same dtype as input.

    Raises
    ------
    ParameterError
        On invalid inputs or when the cross-section is not actually
        sharp.
    """
    if te_thickness <= 0:
        raise ParameterError(f"auto_round_sharp_te: te_thickness must be > 0, got {te_thickness}")
    if span_axis_world not in (0, 1, 2):
        raise ParameterError(
            f"auto_round_sharp_te: span_axis_world must be 0/1/2, got {span_axis_world}"
        )
    delta = cross_section[-1, :] - cross_section[0, :]
    if float(np.linalg.norm(delta)) > 1e-9:
        raise ParameterError(
            f"auto_round_sharp_te: cross-section is not sharp (TE endpoints "
            f"differ by {float(np.linalg.norm(delta)):.6g}); auto-round "
            f"only operates on sharp TE inputs."
        )

    extents = cross_section.max(axis=0) - cross_section.min(axis=0)
    extents[span_axis_world] = -np.inf
    chord_axis = int(np.argmax(extents))
    candidates = [a for a in (0, 1, 2) if a != span_axis_world and a != chord_axis]
    thickness_axis = candidates[0]

    out = cross_section.copy()
    half = 0.5 * te_thickness
    out[0, thickness_axis] += +half
    out[-1, thickness_axis] += -half
    return out
