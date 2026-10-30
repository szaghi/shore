---
title: shore.split
outline: [2, 3]
---

# `shore.split`

Split structured blocks along ijk vertex indices, preserving adjacency.

Why splitting matters
---------------------

CFD solvers and parallel decomposers often need *smaller* blocks than
the topology builder naturally produces — for load balancing across
ranks, GPU streaming-multiprocessor occupancy, halo-exchange
granularity, or just to fit a single block within an L2 cache.  SHORE
emits its blocks at the natural topological scale (one cubed-sphere
sub-block, one flat-caps sub-block, etc.); this module is the
post-processing pass that subdivides those blocks into chunks while
recording every internal and external SHARED adjacency that results.

What this is, what it isn't
---------------------------

This is a *vertex-index* splitter.  Given a list of
:class:`shore.mesh.block.HexBlock` instances and a list of
:class:`BlockSplit` directives — each saying "split block X along
axis A at vertex indices [i1, i2, ...]" — return a fresh list of new
blocks plus all the face BC plumbing.  The original block list is
untouched (deep-copied for safety).

The original starting topology is preserved in spirit: if you had a
6-block cubed sphere, splitting `sub0` into 2 chunks gives you 7
blocks (2 sub0 chunks + sub1 + sub2 + sub3 + cap_north + cap_south).
What changes is only the internal block count and the adjacency
graph; world coordinates of every node are bit-exactly preserved.

Capabilities:

- **Multi-axis in one call.**  A block may be cut on i, j and k at once
  (one :class:`BlockSplit` per axis, same label); it becomes an
  `ni x nj x nk` chunk grid.
- **Automatic paired co-splits.**  A cut that subdivides a SHARED /
  DIRICHLET seam is propagated through the seam's AxisMap onto the
  partner block (to a fixpoint), so the caller never authors paired
  directives and the two sides stay conformal.  Cut *indices* on the
  two sides need not be equal — they are mapped (with flip) through the
  AxisMap.

Out of scope (deliberate):

- Non-vertex-aligned splits (interpolation between existing nodes).
- Splitting periodic faces (cubed-sphere sub-block ↔ sub-block
  j-seams).  The sub-blocks are aliased numpy views into a contiguous
  ring; splitting fights with that aliasing and is rejected.

Adjacency rules
---------------

A block split into a chunk grid exposes two kinds of seams:

- **Internal seams**: the cut planes themselves.  Adjacent chunks are
  wired SHARED with identity AxisMaps and bit-exact node coincidence
  (the chunks alias the same column from the source via numpy views).
- **External seams**: each original face is inherited only by the chunks
  that still expose it; an original inter-block SHARED/DIRICHLET seam is
  re-paired chunk-sheet to chunk-sheet onto the matching partner chunks
  (each side keeps its own BC — DIRICHLET stays DIRICHLET).

## `BlockSplit` <Badge type="tip" text="dataclass" /> {#BlockSplit}

```python
BlockSplit(label: str, axis: str, at)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/split.py#L144">source</a>

One split directive: which block, which axis, at which indices.

**Parameters**

- **`label`** — Label of the block to split (must match a block's `label` in the input list).
- **`axis`** — Axis to split along: `"i"`, `"j"`, or `"k"`.
- **`at`** — Sequence of *vertex* indices at which to cut.  Each index must satisfy `1 <= idx <= n - 1` along the chosen axis (so each resulting chunk has at least 2 vertices in that direction). Indices are sorted and deduplicated; `[32]` and `[16, 32]` produce 2 and 3 chunks respectively.

**Examples**

```python
>>> BlockSplit(label="sub0", axis="i", at=[32])
# BlockSplit(label='sub0', axis='i', at=(32,))
>>> BlockSplit(label="background_center", axis="k", at=[10, 20, 30])
# BlockSplit(label='background_center', axis='k', at=(10, 20, 30))
```

## `split_blocks` <Badge type="info" text="function" /> {#split_blocks}

```python
split_blocks(
    blocks: list[HexBlock],
    splits: list[BlockSplit],
    *,
    propagate_seams: bool = True,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/split.py#L191">source</a>

Apply split directives and return a fresh, fully-wired block list.

A block may be cut on **several axes in one call**, and cuts that
subdivide a SHARED / DIRICHLET seam are **auto-propagated to the
partner block** so the caller never has to author paired co-splits.
The result is a conformal multi-block decomposition with every
internal cut-plane and every original inter-block seam wired (partner
pointers + AxisMaps) onto the right chunks.

**Parameters**

- **`blocks`** — Input block list (e.g. `Mesh.blocks` or a primitive's `.build`).  Not mutated (chunks alias node memory read-only).
- **`splits`** — :class:`BlockSplit` directives.  Multiple directives may name the same block (one per axis); their cuts are unioned per axis.
- **`propagate_seams`** — When True (default), a cut on an in-plane axis of a SHARED / DIRICHLET face is mirrored through the seam's AxisMap onto the partner block (to a fixpoint), so both sides of every seam are cut consistently.  Set False only to reproduce the legacy "caller must pre-pair every co-split" behaviour.

**Returns**

- `list[HexBlock]` — Fresh blocks.  An original block split into an `ni x nj x nk` chunk grid is replaced (in row-major i-major-then-j-then-k order) by its chunks, labelled `<label>__i<a>__j<b>__k<c>` (only the cut axes appear in the suffix); unsplit blocks keep their label and original position.

**Raises**

- `ParameterError` — Unknown label, out-of-range cut, or a seam whose two sides cannot be made conformal (e.g. a partner with an incompatible extent).

## `propagate_cuts` <Badge type="info" text="function" /> {#propagate_cuts}

```python
propagate_cuts(
    blocks: list[HexBlock],
    requested: dict[str, dict[str, set[int]]],
) -> dict[str, dict[str, set[int]]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/split.py#L339">source</a>

Return *requested* cuts after seam-propagation to a fixpoint.

Pure helper for the balance planner: given per-block per-axis cut sets
(vertex indices), return a deep copy extended with every cut induced
across SHARED/DIRICHLET seams — exactly the set :func:`split_blocks`
would execute.  Lets the planner reason about the *actual* chunking
(and enforce the multigrid cell floor) on the propagated result, so its
plan matches the split output block-for-block.

## `chunk_axis_sizes` <Badge type="info" text="function" /> {#chunk_axis_sizes}

```python
chunk_axis_sizes(n_vertices: int, at: set[int]) -> list[int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/split.py#L361">source</a>

Cell counts of each chunk produced by cutting an axis at *at*.

## `original_blocks_face` <Badge type="info" text="function" /> {#original_blocks_face}

```python
original_blocks_face(blocks: list[HexBlock], label: str, face: FaceName)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/split.py#L550">source</a>

Yield the Face object for *label*.*face* among *blocks* (0 or 1 items).
