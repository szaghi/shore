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
) -> list[HexBlock]
```

Apply every split directive and return a fresh list of new blocks.
Blocks not named in `splits` pass through unchanged (deep-copied).

**Adjacency rules:**

- **Internal seams.** Each split block is replaced by `len(at) + 1`
  chunks; consecutive chunks share a `FaceBC.SHARED` face along the
  split axis with **bit-exact node coincidence** — the chunks alias
  the same column from the source block via numpy views, no node
  duplication.
- **External faces** *not* perpendicular to the split axis appear on
  every chunk, sliced into pieces along the perpendicular axis. The
  original block's BC is inherited.
- **Faces perpendicular** to the split axis appear only on the first
  / last chunk respectively.
- **Cross-block partner rewiring.** When the original block had a
  `SHARED` face whose partner was *another* block in `blocks`, the
  partner's face BC is updated to point at the right new chunk.
  Cross-block splits along compatible axes must use **identical
  indices**; perpendicular cross-block splits raise
  `ParameterError`.

**Block label format:** `<original_label>__<axis><k>` where `k`
counts chunks in axis-ascending order. E.g. `sub0__k0`, `sub0__k1`.

**Raises**

- `ParameterError` — unknown block label, out-of-range cut index,
  duplicate split of the same block in one call, or incompatible
  cross-block splits across a SHARED seam.

```python
from shore.split import BlockSplit, split_blocks
from shore.io.adjacency import write_adjacency_json

new_blocks = split_blocks(
    mesh.blocks,
    [
        BlockSplit("sub0", "k", at=[10]),
        BlockSplit("sub1", "k", at=[10]),
        BlockSplit("sub2", "k", at=[10]),
        BlockSplit("sub3", "k", at=[10]),
    ],
)

# Write the per-block .geo + the adjacency JSON sidecar
for blk in new_blocks:
    ...
write_adjacency_json("mesh.adjacency.json", new_blocks)
```

The CLI `shore split` writes the same adjacency sidecar that
[`shore cc-par`](./cc-par) consumes; see the [`shore split` CLI](./split)
for the TOML schema and the JSON sidecar format.

## See also

- [`shore split` CLI](./split) — the user-facing wrapper.
- [`shore.io.adjacency`](./api-adjacency) — JSON sidecar I/O.
- [`shore.mesh.topology`](./api-topology) — `HexBlock`, `Face`,
  `AxisMap`.
- [`shore cc-par`](./cc-par) — the downstream consumer of the
  adjacency sidecar.
