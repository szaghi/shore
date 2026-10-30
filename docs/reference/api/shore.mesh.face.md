---
title: shore.mesh.face
outline: [2, 3]
---

# `shore.mesh.face`

Face and FaceBC: face descriptors for structured blocks.

## `FaceBC` <Badge type="tip" text="class" /> {#FaceBC}

```python
class FaceBC(enum.Enum)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L34">source</a>

Boundary condition type for a block face.

## `AxisMap` <Badge type="tip" text="dataclass" /> {#AxisMap}

```python
AxisMap(along: tuple[tuple[AxisLabel, bool], tuple[AxisLabel, bool]])
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L46">source</a>

Mapping from a face's along-seam axes onto its partner's along-seam axes.

A face perpendicular to one block axis has two *along-seam* axes
lying in the face plane (e.g. an `i_lo` face has along-axes
`(j, k)`).  When the face is SHARED / DIRICHLET / PERIODIC with
another face on a partner block, every consumer (cc.par writer,
overset connectivity, debug viewers, etc.) needs to know how the
parametric (i, j, k) of one side maps onto the partner's (i, j, k):
same axis or different axis, increasing or reversed.

`along` is a 2-tuple of `(partner_axis, flipped)` pairs, in the
order given by :data:`FACE_ALONG_AXES` for *this face's* name.
For an `i_lo` face whose `j`-row coincides with the partner's
`i`-row in the *same* direction, and whose `k`-row coincides
with the partner's `k`-row in the same direction, the AxisMap is
`along=(("i", False), ("k", False))`.

Two consistency requirements are enforced at construction:

- The two `partner_axis` entries must differ (a 2D face cannot
  degenerate onto a single 1D partner axis).
- The two `partner_axis` entries must be drawn from the partner
  face's own along-axes — but that is verified at wiring time by
  :meth:`AxisMap.validate_against_partner`, since :class:`AxisMap`
  itself is decoupled from any specific face name.

**Attributes**

- **`along`** — `((partner_axis_for_first_along, flipped), (partner_axis_for_second_along, flipped))`.

**Examples**

Same-axis seam, no flips (the most common case):

```
AxisMap.identity("i_lo")  # -> along=(("j", False), ("k", False))

```
Cross-axis seam with one reversal (e.g. an `i_lo` face whose
`j`-row maps to the partner's reversed `i`-row, `k` preserved):

```
AxisMap.from_pairs(("i", True), ("k", False))
```

### `AxisMap.from_pairs` <Badge type="info" text="classmethod" /> {#AxisMap-from_pairs}

```python
from_pairs(first: tuple[AxisLabel, bool], second: tuple[AxisLabel, bool]) -> AxisMap
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L109">source</a>

Construct an :class:`AxisMap` from two `(axis, flipped)` pairs.

### `AxisMap.identity` <Badge type="info" text="classmethod" /> {#AxisMap-identity}

```python
identity(face_name: FaceName) -> AxisMap
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L118">source</a>

Same-axis, no-flip mapping for *face_name*.

Maps each of *face_name*'s along-axes to the partner's
same-named axis with no reversal — appropriate when both blocks
share node arrays in identical (i, j, k) orientation (the most
common case: same-axis SHARED seams).

### `AxisMap.partner_axes` <Badge type="info" text="method" /> {#AxisMap-partner_axes}

```python
partner_axes() -> tuple[AxisLabel, AxisLabel]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L133">source</a>

Return the two partner along-axis labels (without the flip flags).

### `AxisMap.flips` <Badge type="info" text="method" /> {#AxisMap-flips}

```python
flips() -> tuple[bool, bool]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L137">source</a>

Return the two flip flags (True = reversed indexing).

### `AxisMap.inverse` <Badge type="info" text="method" /> {#AxisMap-inverse}

```python
inverse(this_face: FaceName, partner_face: FaceName) -> AxisMap
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L141">source</a>

Return the AxisMap for the *partner* side of the same seam.

If `self` describes how *this_face*'s along-axes map onto
*partner_face*'s along-axes, the inverse describes how
*partner_face*'s along-axes map onto *this_face*'s along-axes.
Both sides agree on flip parity (a flip is symmetric).

### `AxisMap.validate_against_partner` <Badge type="info" text="method" /> {#AxisMap-validate_against_partner}

```python
validate_against_partner(this_face: FaceName, partner_face: FaceName) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L163">source</a>

Check that `self.partner_axes()` are exactly *partner_face*'s along-axes.

**Raises**

- `ValueError` — If any partner-axis label in *self* is not an along-axis of *partner_face*, or if the two label slots are not a permutation of *partner_face*'s along-axes.

## `Face` <Badge type="tip" text="dataclass" /> {#Face}

```python
Face(
    index: int,
    name: FaceName,
    bc: FaceBC,
    partner: Face | None = field(default=None, repr=False, compare=False),
    axis_map: AxisMap | None = field(default=None, repr=False),
    bc_code: int | None = None,
    association: Association | None = field(default=None, repr=False),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L189">source</a>

One face of a structured block.

**Parameters**

- **`index`** — 0-based index within the block's face list.
- **`name`** — Canonical face name: `'i_lo'`, `'i_hi'`, `'j_lo'`, `'j_hi'`, `'k_lo'` (wall), `'k_hi'` (outer).
- **`bc`** — Boundary condition type.
- **`partner`** — The partnered Face for SHARED / DIRICHLET / PERIODIC BCs.  None for WALL and FREE.
- **`axis_map`** — Mapping from this face's along-axes onto the partner's along-axes.  Required when `partner` is set; `None` otherwise.  See :class:`AxisMap`.
- **`bc_code`** — Optional explicit downstream natural-BC code for a `FREE` face, overriding the per-block default the cc.par writer would otherwise apply (`BlockMeta.free_face_bc` / `--free-bc`). `None` means "use the per-block default".  Only meaningful on `FREE` faces; must be a positive Xall natural-BC code in `[1, nccnat=19]` (e.g. `9` free/far-field, `3`/`4` inflow/outflow).  Lets a single block carry heterogeneous outer BCs — inflow on one face, extrapolation on the others.
- **`association`** — Geometry this face's nodes are projected onto (ADR 0003), or `None`.

### `Face.body_projected` <Badge type="info" text="property" /> {#Face-body_projected}

```python
body_projected: bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L234">source</a>

True when associated with a whole body (legacy flag, read-only).

### `Face.nodes` <Badge type="info" text="method" /> {#Face-nodes}

```python
nodes(block: Block) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/face.py#L238">source</a>

Return a view of this face's node coordinates in block.nodes.
