---
title: shore.topology.oh_channel
outline: [2, 3]
---

# `shore.topology.oh_channel`

Channel-block frustum builder for the OH (cubed-sphere + H-channel) topology.

A *channel block* sits between one wrap block's outer envelope (k=nk_w-1) and
one face of the user-supplied bounding box.  The two end faces are:

- `k=0`     `inner_face` — bit-copy of the wrap's `k_hi` slice.
- `k=nk-1`  bilinear (a, b) sweep of `outer_quad` (the box face's 4 corners).

Interior k-layers are produced by ruled-surface interpolation between the
two end faces under a 1D spacing law (`uniform` / `tanh` / `tanh2`).
The outer-face parameterisation is index-for-index aligned with the inner
face so that the channel-channel lateral seams (j_hi ↔ j_lo on adjacent
sub-blocks, and i-seams to the cap channels) inherit the wrap's seam
topology with identity AxisMaps.

C1 continuity at the wrap-channel seam is optional via `c1_first_step`:
when set, the first interior layer is placed at distance `c1_first_step`
from k=0 along the local k-direction, matching the wrap's last marching
step exactly so that the spacing has no jump across the seam.

## `build_channel_block` <Badge type="info" text="function" /> {#build_channel_block}

```python
build_channel_block(
    inner_face: np.ndarray,
    outer_quad: np.ndarray,
    nk_channel: int,
    k_law: str = 'uniform',
    k_beta: float = 3.0,
    c1_first_step: float | np.ndarray | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/topology/oh_channel.py#L30">source</a>

Build a channel frustum block, shape `(ni, nj, nk_channel, 3)`.

**Parameters**

- **`inner_face`** — Wrap block's k_hi outer envelope, shape `(ni, nj, 3)`.  This becomes the channel's k=0 row, bit-exactly.
- **`outer_quad`** — Four 3D corners of the box face this channel block "owns", shape `(4, 3)`.  Order: `(a=0, b=0)`, `(a=1, b=0)`, `(a=1, b=1)`, `(a=0, b=1)`.  The same parametric (a, b) used by `surface/cap.py`, so that diagonal sub-block junctions in the wrap map onto box-face corners.
- **`nk_channel`** — Total number of layers in the channel block (k=0 + interior + k=nk-1). Must be `>= 2`.
- **`k_law`** — Spacing law along k: `'uniform'`, `'tanh'` (cluster near inner face / wrap seam), or `'tanh2'` (cluster at both ends).
- **`k_beta`** — Clustering strength for tanh / tanh2.
- **`c1_first_step`** — When supplied as a scalar, the distance from k=0 to k=1 at every (i, j) equals this value exactly.  When supplied as an `(ni, nj)` array, the distance at each (i, j) equals the corresponding entry — useful for matching the wrap's per-(i, j) last marching step exactly so the seam is C1 at every node, not just on average.  Must be positive everywhere it's supplied.

**Returns**

- `numpy.ndarray, shape ``(ni, nj, nk_channel, 3)``` — Frustum block in block-local (i, j, k) layout.

**Raises**

- `ParameterError` — On invalid shapes, unknown spacing laws, or non-positive `c1_first_step`.
