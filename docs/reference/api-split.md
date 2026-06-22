# shore.split

Subdivide a list of [`HexBlock`](./api-topology#hexblock) instances
along ijk vertex indices, preserving the adjacency graph. Used by
the [`shore split`](./split) CLI; also reusable from Python when
the parallel-decomposition step lives inside a larger workflow.

## `BlockSplit`

```python
@dataclass
class BlockSplit:
    label: str
    axis:  Literal["i", "j", "k"]
    at:    tuple[int, ...]
```

One split directive: which block to subdivide, along which axis,
at which vertex indices. Cuts are sorted and deduplicated at
construction.

| Parameter | Description |
|---|---|
| `label` | Label of the block to split. Validated against the actual block list at `split_blocks` time. |
| `axis` | `"i"` / `"j"` / `"k"`. |
| `at` | Sequence of vertex indices. Each index must satisfy `1 <= idx <= n - 1` so every resulting chunk has at least 2 vertices on the split axis. `N` indices → `N + 1` chunks. |

```python
from shore.split import BlockSplit

# 2 chunks: sub0__i0, sub0__i1
BlockSplit("sub0", "i", at=[32])

# 4 chunks: background_center__k0, ..., __k3
BlockSplit("background_center", "k", at=[10, 20, 30])
```

## `split_blocks`

```python
def split_blocks(
    blocks: list[HexBlock],
    splits: list[BlockSplit],
    *,
    propagate_seams: bool = True,
) -> list[HexBlock]
```

Apply the split directives and return a fresh, fully-wired list of new blocks.
Blocks not named in `splits` pass through unchanged. **Multiple directives may
name the same block** (one per axis); their cuts are unioned per axis, so a
block becomes an `ni × nj × nk` chunk grid in a single pass.

| Parameter | Description |
|---|---|
| `blocks` | Input block list (not mutated; chunks alias node memory read-only). |
| `splits` | `BlockSplit` directives; multiple per block allowed (multi-axis). |
| `propagate_seams` | When `True` (default), a cut on an in-plane axis of a SHARED / DIRICHLET face is mirrored through the seam's `AxisMap` onto the partner block (to a fixpoint), so both sides stay conformal — callers need not author paired co-splits. Set `False` for the legacy "caller pre-pairs every co-split" behaviour. |

**Adjacency rules:**

- **Internal seams.** Adjacent chunks share a `FaceBC.SHARED` face with
  **bit-exact node coincidence** (numpy-view aliased, no duplication) and an
  identity `AxisMap`.
- **External seams.** Each original face is inherited only by the chunks that
  still expose it. An original inter-block SHARED / DIRICHLET seam is re-paired
  chunk-sheet to chunk-sheet onto the matching partner chunks — each side keeps
  its own BC (a DIRICHLET cap↔sub seam stays DIRICHLET). Per-face `bc_code`
  overrides are carried through.
- **Auto seam co-split.** With `propagate_seams=True`, a cut that subdivides a
  seam induces the corresponding cut on the partner (mapped through the
  `AxisMap`, with flip and partner extent); cut indices on the two sides need
  not be equal. See
  [Design decisions](/guide/design-decisions#splitting-auto-propagates-cuts-across-seams).

**Block label format:** `<original_label>__<axis><k>` per cut axis, e.g.
`sub0__i0__k1`. Only the cut axes appear in the suffix.

**Raises**

- `ParameterError` — unknown block label, out-of-range cut index, or a seam
  whose two sides cannot be made conformal.

```python
from shore.split import BlockSplit, split_blocks

# Cut sub0 on k only — auto-propagates to the j-ring sub-blocks and the
# DIRICHLET caps, so every cubed-sphere block ends up k-split conformally.
new_blocks = split_blocks(mesh.blocks, [BlockSplit("sub0", "k", at=[10])])
```

## `propagate_cuts`

```python
def propagate_cuts(
    blocks: list[HexBlock],
    requested: dict[str, dict[str, set[int]]],
) -> dict[str, dict[str, set[int]]]
```

Return `requested` (per-block per-axis vertex-cut sets) extended with every cut
induced across SHARED / DIRICHLET seams — exactly the set `split_blocks` would
execute. Pure helper for the balance planner, which uses it to floor-check the
*propagated* chunking so its plan matches the split output block-for-block.

## `chunk_axis_sizes`

```python
def chunk_axis_sizes(n_vertices: int, at: set[int]) -> list[int]
```

Cell counts of each chunk produced by cutting an axis of `n_vertices` vertices
at the vertex indices `at`. Used with `propagate_cuts` to enforce the multigrid
cell floor.

## See also

- [`shore split` CLI](./split) — the user-facing wrapper.
- [`shore.io.adjacency`](./api-adjacency) — JSON sidecar I/O.
- [`shore.mesh.topology`](./api-topology) — `HexBlock`, `Face`,
  `AxisMap`.
- [`shore cc-par`](./cc-par) — the downstream consumer of the
  adjacency sidecar.
