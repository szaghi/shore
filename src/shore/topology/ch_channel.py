"""Outer-quad helpers for the CH (C-grid + H-channel) topology.

The CH topology bridges a marched C-grid wrap (3-block sharp-TE or
4-block blunt-TE) to an axis-aligned far-field rectangular frame in
the chord/normal plane (extruded across span).  Each channel block is
a hexahedral frustum with:

- ``k=0``     bit-aliased view of one wrap face (e.g. ``upper.k_hi``).
- ``k=nk-1``  bilinear sweep of 4 corners on the frame.
- interior   ruled-surface interpolation by
             :func:`shore.topology.oh_channel.build_channel_block`.

This module computes the **planar outer quads** (4 world-frame corners
each) for the 7 channels of CH-3 and 8 channels of CH-4, given the
wrap envelope's already-marched k_hi corners and the user's
chord/normal box extents.

The corner ordering ``(c00, c10, c11, c01)`` of every quad matches
:func:`build_channel_block`'s ``(a=0,b=0), (a=1,b=0), (a=1,b=1),
(a=0,b=1)`` convention, where ``a`` indexes the channel's i-axis and
``b`` indexes its j-axis (= the wrap's span-axis).  This guarantees
that channel-channel lateral seams (j↔j) and chord-direction seams
(i↔i) inherit the wrap's seam topology with identity AxisMaps.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from shore.exceptions import ParameterError


@dataclass(frozen=True)
class CHFrame:
    """Cached chord/normal/span axis indices and box extents for a CH frame.

    Attributes
    ----------
    chord_axis, normal_axis, span_axis_idx :
        World-axis indices (0=x, 1=y, 2=z) of the chord, normal, and
        span directions.  All three differ.
    chord_sign :
        ``+1`` if the freestream's chord component is positive (downstream
        = +chord), ``-1`` if negative.  Channel ``chan_outflow_*`` lives
        on the +chord side regardless: the convention is "outflow = far
        end of wake", and the wake direction is freestream-aligned.
    chord_min, chord_max, normal_min, normal_max :
        Frame extents in world coordinates (already chord/normal sign
        sorted: ``chord_min < chord_max``, ``normal_min < normal_max``).
    """

    chord_axis: int
    normal_axis: int
    span_axis_idx: int
    chord_sign: int
    chord_min: float
    chord_max: float
    normal_min: float
    normal_max: float


def resolve_frame(
    *,
    freestream_world: np.ndarray | None,
    span_axis: str,
    chord_min: float,
    chord_max: float,
    normal_min: float,
    normal_max: float,
) -> CHFrame:
    """Resolve world-axis indices and validate the box extents.

    The chord axis is the larger-magnitude non-span component of
    ``freestream_world``; the normal axis is the remaining non-span
    axis.  ``freestream_world=None`` defaults to ``+chord`` along
    whichever non-span world axis comes first in the (x, y, z) order.

    Raises
    ------
    ParameterError
        On invalid inputs.
    """
    if span_axis not in ("x", "y", "z"):
        raise ParameterError(f"resolve_frame: span_axis must be 'x' / 'y' / 'z', got {span_axis!r}")
    span_idx = "xyz".index(span_axis)
    non_span = [a for a in (0, 1, 2) if a != span_idx]

    if freestream_world is None:
        chord_axis = non_span[0]
        chord_sign = 1
    else:
        fs = np.asarray(freestream_world, dtype=float).reshape(-1)
        if fs.size != 3:
            raise ParameterError(
                f"resolve_frame: freestream_world must have length 3, got {fs.size}"
            )
        if np.linalg.norm(fs) < 1e-12:
            raise ParameterError("resolve_frame: freestream_world has zero magnitude")
        # Pick the non-span axis with the largest |fs| component as chord.
        comps = [abs(float(fs[a])) for a in non_span]
        if max(comps) < 1e-12:
            raise ParameterError(
                "resolve_frame: freestream_world is parallel to the span axis — "
                "wake direction is undefined.  Supply a freestream with a "
                "non-zero in-plane component."
            )
        chord_axis = non_span[int(np.argmax(comps))]
        chord_sign = 1 if float(fs[chord_axis]) >= 0 else -1
    normal_axis = next(a for a in non_span if a != chord_axis)

    if not (chord_max > chord_min):
        raise ParameterError(
            f"resolve_frame: chord_max must exceed chord_min; got "
            f"chord_min={chord_min}, chord_max={chord_max}"
        )
    if not (normal_max > normal_min):
        raise ParameterError(
            f"resolve_frame: normal_max must exceed normal_min; got "
            f"normal_min={normal_min}, normal_max={normal_max}"
        )

    return CHFrame(
        chord_axis=chord_axis,
        normal_axis=normal_axis,
        span_axis_idx=span_idx,
        chord_sign=chord_sign,
        chord_min=float(chord_min),
        chord_max=float(chord_max),
        normal_min=float(normal_min),
        normal_max=float(normal_max),
    )


def _make_world_point(
    *,
    chord: float,
    normal: float,
    span: float,
    frame: CHFrame,
) -> np.ndarray:
    """Assemble a 3D world point from chord/normal/span coordinates."""
    p = np.zeros(3, dtype=np.float64)
    p[frame.chord_axis] = chord
    p[frame.normal_axis] = normal
    p[frame.span_axis_idx] = span
    return p


def _quad(
    *,
    c0: tuple[float, float],
    c1: tuple[float, float],
    span_lo: float,
    span_hi: float,
    frame: CHFrame,
) -> np.ndarray:
    """Build a (4, 3) outer quad in world coords.

    A channel block's i-axis runs from frame point ``c0=(chord, normal)``
    to ``c1=(chord, normal)`` (in the chord/normal plane), and its
    j-axis runs from ``span_lo`` to ``span_hi``.  The 4 corners are in
    ``(a=0,b=0), (a=1,b=0), (a=1,b=1), (a=0,b=1)`` order.
    """
    c00 = _make_world_point(chord=c0[0], normal=c0[1], span=span_lo, frame=frame)
    c10 = _make_world_point(chord=c1[0], normal=c1[1], span=span_lo, frame=frame)
    c11 = _make_world_point(chord=c1[0], normal=c1[1], span=span_hi, frame=frame)
    c01 = _make_world_point(chord=c0[0], normal=c0[1], span=span_hi, frame=frame)
    return np.stack([c00, c10, c11, c01], axis=0)


def outer_quads_for_cgrid(
    *,
    frame: CHFrame,
    upper_envelope: np.ndarray,
    body_envelope: np.ndarray,
    lower_envelope: np.ndarray,
    i_split_body_upper: int,
    i_split_body_lower: int,
    span_lo: float,
    span_hi: float,
    is_blunt: bool,
    wake_base_envelope: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Compute the 7 (or 8) planar outer quads for a CH topology.

    Each entry is a ``(4, 3)`` world-coord array in the
    ``(a=0,b=0), (a=1,b=0), (a=1,b=1), (a=0,b=1)`` corner order
    expected by :func:`build_channel_block`.  ``a`` is the channel's
    i-axis (along the frame edge) and ``b`` is the channel's j-axis
    (along span).

    Parameters
    ----------
    frame :
        Resolved :class:`CHFrame` (chord/normal axes + box extents).
    upper_envelope, body_envelope, lower_envelope :
        ``(ni_block, n_stations, 3)`` arrays of the wrap blocks' k_hi
        slice **after marching**.  Used to anchor the in-plane chord
        coordinate of each channel along the frame edge.
    i_split_body_upper, i_split_body_lower :
        Body-block i-indices at which ``body.k_hi`` is split into the
        ``chan_body_top`` / ``chan_body_front`` / ``chan_body_bot``
        channels.  Must satisfy ``0 < i_split_body_upper <
        i_split_body_lower < ni_body - 1``.
    span_lo, span_hi :
        Spanwise extents of the channel blocks (= the wrap's span
        extents at ``j_lo`` / ``j_hi``, taken from the same envelope).
    is_blunt :
        ``True`` for CH-4 blunt-TE wrap (yields an extra
        ``chan_outflow_wake_base`` quad).
    wake_base_envelope :
        ``(ni_wake, n_stations, nk_w, 3)`` wake-base block.  Required
        when ``is_blunt`` to anchor the third outflow quad's normal
        extents.

    Returns
    -------
    dict[str, ndarray]
        Keyed by channel label.  CH-3 keys: ``"upper"``, ``"body_top"``,
        ``"body_front"``, ``"body_bot"``, ``"lower"``,
        ``"outflow_upper"``, ``"outflow_lower"``.  CH-4 adds
        ``"outflow_wake_base"``.

    Notes
    -----
    The "chord positions" along the top / bottom frame edges are taken
    from the actual envelope's chord component at each split point.
    Adjacent channels share a quad-corner: ``chan_upper`` ends at
    chord = ``upper_envelope[-1, :, chord_axis]`` (= upper-TE
    envelope), and ``chan_body_top`` starts at the same chord — bit-
    exactly when ``body_envelope[0, :, :]`` is the SHARED upper-TE
    column with ``upper_envelope[-1, :, :]`` (which it is via numpy
    view aliasing in the wrap allocation).
    """
    if not (0 < i_split_body_upper < i_split_body_lower < body_envelope.shape[0] - 1):
        raise ParameterError(
            f"outer_quads_for_cgrid: invalid body splits "
            f"(0 < {i_split_body_upper} < {i_split_body_lower} < "
            f"{body_envelope.shape[0] - 1}) violated"
        )
    if is_blunt and wake_base_envelope is None:
        raise ParameterError(
            "outer_quads_for_cgrid: wake_base_envelope is required when is_blunt=True"
        )

    cax = frame.chord_axis
    nax = frame.normal_axis

    # Per-channel chord anchor positions taken from the envelope.  We
    # collapse along span by averaging — the frame is span-uniform, and
    # the wrap envelope's chord coordinate varies negligibly across span
    # at the same i-index (stations are independent C-curves but they
    # share a freestream and TE).
    def _chord_at(env: np.ndarray, i: int) -> float:
        return float(env[i, :, cax].mean())

    # Body-side chord positions at the 3 internal splits + the 2 TEs.
    chord_upper_te = _chord_at(body_envelope, 0)  # = upper.k_hi[-1] chord
    chord_body_split_upper = _chord_at(body_envelope, i_split_body_upper)
    chord_body_split_lower = _chord_at(body_envelope, i_split_body_lower)
    chord_lower_te = _chord_at(body_envelope, body_envelope.shape[0] - 1)

    # Wake far-tip chord positions on the upper / lower wakes.  For a
    # standard +x freestream these sit near chord_max.  We DON'T pin
    # them to chord_max — instead we let the channel's outer quad
    # anchor at the actual envelope chord; the channel's i-extent
    # along the frame edge then matches the wrap's i-extent and the
    # outflow quad picks up the residual.
    chord_upper_far = _chord_at(upper_envelope, 0)
    chord_lower_far = _chord_at(lower_envelope, lower_envelope.shape[0] - 1)

    # The frame edges live at:
    #   top    : normal = frame.normal_max, chord ∈ [chord_min, chord_max]
    #   bottom : normal = frame.normal_min, chord ∈ [chord_min, chord_max]
    #   front  : chord  = frame.chord_min,  normal ∈ [normal_min, normal_max]
    #   rear   : chord  = frame.chord_max,  normal ∈ [normal_min, normal_max]
    # (with chord_max on the +freestream side).
    top = frame.normal_max
    bot = frame.normal_min
    front = frame.chord_min
    rear = frame.chord_max

    quads: dict[str, np.ndarray] = {}

    # chan_upper: upper.k_hi → top frame edge.
    # i runs upper-far-tip (envelope-chord chord_upper_far) → upper-TE
    # (chord_upper_te).  Outer face on top edge between the same chord
    # interval.
    quads["upper"] = _quad(
        c0=(chord_upper_far, top),
        c1=(chord_upper_te, top),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )

    # chan_body_top: body.k_hi[0:i_split_upper+1] → top frame edge.
    # Outer chord interval [chord_upper_te, chord_min] on the top edge.
    # The shared column with chan_body_front sits at the top-front
    # frame corner (chord_min, normal_max), so chan_body_front can
    # start from the same corner — both planar quads agree at the
    # overlap.  The i_split_body_upper position only controls the wrap-
    # side i-extent; the outer face's endpoint is pinned to the corner
    # regardless of where on the body the split lands.
    quads["body_top"] = _quad(
        c0=(chord_upper_te, top),
        c1=(front, top),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )
    # Suppress unused-variable lint: split chord positions retained for
    # documentation purposes (and possible future use in non-corner
    # split modes).
    _ = chord_body_split_upper, chord_body_split_lower

    # chan_body_front: body.k_hi[i_split_upper:i_split_lower+1] → front
    # frame edge.  Outer face sweeps from the top-front corner down to
    # the bottom-front corner, both on chord = chord_min.
    quads["body_front"] = _quad(
        c0=(front, top),
        c1=(front, bot),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )

    # chan_body_bot: body.k_hi[i_split_lower:] → bottom frame edge.
    # Outer chord interval [chord_min, chord_lower_te] on the bottom
    # edge.  Shared column with chan_body_front at the bottom-front
    # corner (chord_min, normal_min).
    quads["body_bot"] = _quad(
        c0=(front, bot),
        c1=(chord_lower_te, bot),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )

    # chan_lower: lower.k_hi → bottom frame edge.
    # i runs lower-TE → lower-far-tip; outer chord interval
    # [chord_lower_te, chord_lower_far] on the bottom edge.
    quads["lower"] = _quad(
        c0=(chord_lower_te, bot),
        c1=(chord_lower_far, bot),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )

    # chan_outflow_upper: a rectangular block in (chord, normal) sitting
    # downstream of the wrap far-tip, ABOVE the wake-cut.  Channel axes:
    #
    #     i (channel)  → chord  (chord_upper_far  → chord_max)
    #     j (channel)  → span
    #     k (channel)  → normal (wake_seam_upper  → normal_max)
    #
    # Inner face (k=0) lies on the wake-cut plane (normal = wake_seam_upper);
    # outer face (k=last) lies on the top frame edge (normal = normal_max).
    # Left face (i=0) is at chord = chord_upper_far (FREE — adjacent to but
    # not coupled with the wrap's upper.i_lo face; downstream cc.par can
    # add chimera fringes there).  Right face (i=last) is at chord =
    # chord_max (rear frame edge, FREE).
    #
    # outer_quads_for_cgrid only needs to express the outer face quad
    # (i=last); the inner face for the outflow channels is built freshly
    # by fill_channel_blocks since there is no wrap face to alias.
    # We still emit one quad per outflow block here, namely the channel's
    # k_hi face (the top / bottom frame edge between far-tip and rear).
    # The channel itself, when assembled, is the bilinear sweep of the
    # 8 corners (k_lo and k_hi quads).  We pack k_lo / k_hi into ONE
    # synthetic quad whose "i-axis" runs chord and whose "j-axis" runs
    # along channel-k from inner to outer.  build_channel_block can't
    # consume that, so for outflow channels we instead use a simpler
    # approach: fill_channel_blocks builds these channels by direct
    # bilinear interpolation in (chord, normal) using the 4 corners
    # below as the channel's bounding rectangle in the chord/normal
    # plane.  The "outer_quad" for the outflow channels is therefore
    # the **rear-frame edge** (top half / bottom half), and the
    # "inner_face" is constructed in fill_channel_blocks from the
    # chord_far-tip line at normal interval (0, normal_max).
    wake_seam_upper = 0.0
    wake_seam_lower = 0.0
    if is_blunt:
        # Blunt TE: the upper wake k=0 strip sits at +te/2 in normal,
        # lower at -te/2.  Use the wake_base envelope to pick these up.
        wake_seam_upper = float(wake_base_envelope[0, :, -1, nax].mean())
        wake_seam_lower = float(wake_base_envelope[0, :, 0, nax].mean())

    # chan_outflow_upper: rectangular block bounded by
    # (chord_upper_far, wake_seam_upper) ↔ (chord_max, normal_max).
    # The "outer_quad" stored here is the **k_hi face** = top frame edge.
    # fill_channel_blocks reads this + the rectangle bounds from the
    # block label to build the full block.
    quads["outflow_upper"] = _quad(
        c0=(chord_upper_far, top),
        c1=(rear, top),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )

    # chan_outflow_lower: bounded by
    # (chord_lower_far, normal_min) ↔ (chord_max, wake_seam_lower).
    # k_hi = bottom frame edge.
    quads["outflow_lower"] = _quad(
        c0=(chord_lower_far, bot),
        c1=(rear, bot),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )

    # Stash the wake-seam normals + chord-far-tip positions in extra
    # quads keyed by sentinel labels so fill_channel_blocks can find them
    # without re-deriving from the envelope.
    quads["__outflow_upper_inner__"] = _quad(
        c0=(chord_upper_far, wake_seam_upper),
        c1=(rear, wake_seam_upper),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )
    quads["__outflow_lower_inner__"] = _quad(
        c0=(chord_lower_far, wake_seam_lower),
        c1=(rear, wake_seam_lower),
        span_lo=span_lo,
        span_hi=span_hi,
        frame=frame,
    )

    if is_blunt:
        # chan_outflow_wake_base: thin block bridging the TE gap
        # downstream of the wrap, between (chord_far_tip, wake_seam_lower)
        # and (chord_max, wake_seam_upper).
        chord_wb = float(wake_base_envelope[0, :, 0, cax].mean())
        quads["outflow_wake_base"] = _quad(
            c0=(chord_wb, wake_seam_upper),
            c1=(rear, wake_seam_upper),
            span_lo=span_lo,
            span_hi=span_hi,
            frame=frame,
        )
        quads["__outflow_wake_base_inner__"] = _quad(
            c0=(chord_wb, wake_seam_lower),
            c1=(rear, wake_seam_lower),
            span_lo=span_lo,
            span_hi=span_hi,
            frame=frame,
        )

    return quads
