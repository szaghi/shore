---
title: shore.topology.ch_channel
outline: [2, 3]
---

# `shore.topology.ch_channel`

Outer-quad helpers for the CH (C-grid + H-channel) topology.

The CH topology bridges a marched C-grid wrap (3-block sharp-TE or
4-block blunt-TE) to an axis-aligned far-field rectangular frame in
the chord/normal plane (extruded across span).  Each channel block is
a hexahedral frustum with:

- `k=0`     bit-aliased view of one wrap face (e.g. `upper.k_hi`).
- `k=nk-1`  bilinear sweep of 4 corners on the frame.
- interior   ruled-surface interpolation by
             :func:`shore.topology.oh_channel.build_channel_block`.

This module computes the **planar outer quads** (4 world-frame corners
each) for the 7 channels of CH-3 and 8 channels of CH-4, given the
wrap envelope's already-marched k_hi corners and the user's
chord/normal box extents.

The corner ordering `(c00, c10, c11, c01)` of every quad matches
:func:`build_channel_block`'s ``(a=0,b=0), (a=1,b=0), (a=1,b=1),
(a=0,b=1)` convention, where `a`` indexes the channel's i-axis and
`b` indexes its j-axis (= the wrap's span-axis).  This guarantees
that channel-channel lateral seams (j↔j) and chord-direction seams
(i↔i) inherit the wrap's seam topology with identity AxisMaps.

## `CHFrame` <Badge type="tip" text="dataclass" /> {#CHFrame}

```python
CHFrame(
    chord_axis: int,
    normal_axis: int,
    span_axis_idx: int,
    chord_sign: int,
    chord_min: float,
    chord_max: float,
    normal_min: float,
    normal_max: float,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/ch_channel.py#L36">source</a>

Cached chord/normal/span axis indices and box extents for a CH frame.

**Attributes**

- **`chord_axis, normal_axis, span_axis_idx`** — World-axis indices (0=x, 1=y, 2=z) of the chord, normal, and span directions.  All three differ.
- **`chord_sign`** — `+1` if the freestream's chord component is positive (downstream = +chord), `-1` if negative.  Channel `chan_outflow_*` lives on the +chord side regardless: the convention is "outflow = far end of wake", and the wake direction is freestream-aligned.
- **`chord_min, chord_max, normal_min, normal_max`** — Frame extents in world coordinates (already chord/normal sign sorted: `chord_min < chord_max`, `normal_min < normal_max`).

## `resolve_frame` <Badge type="info" text="function" /> {#resolve_frame}

```python
resolve_frame(
    *,
    freestream_world: np.ndarray | None,
    span_axis: str,
    chord_min: float,
    chord_max: float,
    normal_min: float,
    normal_max: float,
) -> CHFrame
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/ch_channel.py#L64">source</a>

Resolve world-axis indices and validate the box extents.

The chord axis is the larger-magnitude non-span component of
`freestream_world`; the normal axis is the remaining non-span
axis.  `freestream_world=None` defaults to `+chord` along
whichever non-span world axis comes first in the (x, y, z) order.

**Raises**

- `ParameterError` — On invalid inputs.

## `outer_quads_for_cgrid` <Badge type="info" text="function" /> {#outer_quads_for_cgrid}

```python
outer_quads_for_cgrid(
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
) -> dict[str, np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/ch_channel.py#L173">source</a>

Compute the 7 (or 8) planar outer quads for a CH topology.

Each entry is a `(4, 3)` world-coord array in the
`(a=0,b=0), (a=1,b=0), (a=1,b=1), (a=0,b=1)` corner order
expected by :func:`build_channel_block`.  `a` is the channel's
i-axis (along the frame edge) and `b` is the channel's j-axis
(along span).

**Parameters**

- **`frame`** — Resolved :class:`CHFrame` (chord/normal axes + box extents).
- **`upper_envelope, body_envelope, lower_envelope`** — `(ni_block, n_stations, 3)` arrays of the wrap blocks' k_hi slice **after marching**.  Used to anchor the in-plane chord coordinate of each channel along the frame edge.
- **`i_split_body_upper, i_split_body_lower`** — Body-block i-indices at which `body.k_hi` is split into the `chan_body_top` / `chan_body_front` / `chan_body_bot` channels.  Must satisfy `0 < i_split_body_upper < i_split_body_lower < ni_body - 1`.
- **`span_lo, span_hi`** — Spanwise extents of the channel blocks (= the wrap's span extents at `j_lo` / `j_hi`, taken from the same envelope).
- **`is_blunt`** — `True` for CH-4 blunt-TE wrap (yields an extra `chan_outflow_wake_base` quad).
- **`wake_base_envelope`** — `(ni_wake, n_stations, nk_w, 3)` wake-base block.  Required when `is_blunt` to anchor the third outflow quad's normal extents.

**Returns**

- `dict[str, ndarray]` — Keyed by channel label.  CH-3 keys: `"upper"`, `"body_top"`, `"body_front"`, `"body_bot"`, `"lower"`, `"outflow_upper"`, `"outflow_lower"`.  CH-4 adds `"outflow_wake_base"`.

**Notes**

The "chord positions" along the top / bottom frame edges are taken
from the actual envelope's chord component at each split point.
Adjacent channels share a quad-corner: `chan_upper` ends at
chord = `upper_envelope[-1, :, chord_axis]` (= upper-TE
envelope), and `chan_body_top` starts at the same chord — bit-
exactly when `body_envelope[0, :, :]` is the SHARED upper-TE
column with `upper_envelope[-1, :, :]` (which it is via numpy
view aliasing in the wrap allocation).
