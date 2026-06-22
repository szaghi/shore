# shore split

Split structured blocks along ijk vertex indices, preserving the
adjacency graph.

CFD solvers and parallel decomposers often need *smaller* blocks than
the topology builder naturally produces — for load balancing across
ranks, GPU streaming-multiprocessor occupancy, halo-exchange
granularity, or cache fit. SHORE emits its blocks at the natural
topological scale (one cubed-sphere sub-block per quadrant, one
flat-caps sub-block per cardinal direction). `shore split` is the
post-processing pass that subdivides them while recording every
internal and external SHARED adjacency that results.

## Synopsis

```
shore split [OPTIONS] -c CONFIG INPUTS...
```

## Arguments

| Argument | Description |
|---|---|
| `INPUTS...` | One or more input `.geo` files. Block labels are derived from each file's stem; if all inputs share a common prefix (e.g. `sphere_sub0.geo`, `sphere_sub1.geo`), the prefix is stripped. |

## Options

| Option | Description |
|---|---|
| `-c` / `--config PATH` | TOML split config (`kind = "split"`) — required. |
| `-o` / `--output PATH` | Output stem. Writes `<stem>_<label>.geo` per chunk plus `<stem>.adjacency.json`. Default: parent dir of the first input + the inferred base prefix. |
| `--adjacency PATH` | Input `.adjacency.json` sidecar(s) describing the pre-existing inter-block seams of the input `.geo` files. **Repeatable** (one per component). Without it, `.geo` files carry only geometry, so every face defaults to FREE and the topology's SHARED/DIRICHLET seams are lost through the split. Pass each component's sidecar to seed and preserve its seams. |
| `-v` / `--verbose` | Print every input block and every output file. |

::: tip Always pass `--adjacency` for topology meshes
A `.geo` file is geometry only. To split a cubed-sphere, C-grid, or any
multi-block mesh and keep its SHARED/DIRICHLET seams, pass the sidecar that
`shore mesh --adjacency-json` (or a previous `shore split`) wrote. The splitter
seeds the input blocks' faces from it, then carries the seams onto the chunks.
:::

## Output

For each chunk produced by the splits:

- One `.geo` file at `<output_stem>_<chunk_label>.geo`.

Plus exactly one JSON sidecar at `<output_stem>.adjacency.json` describing
the BC of every face on every block, including partner references for
`SHARED` / `DIRICHLET` faces. See [Adjacency JSON](#adjacency-json) below.

Block labels of the output chunks follow `<original_label>__<axis><k>`
where `axis ∈ {i, j, k}` and `k` counts chunks in axis-ascending
order, e.g. `sub0__k0`, `sub0__k1`.

## TOML schema (`kind = "split"`)

```toml
version = 1
kind    = "split"

# Optional output stem; CLI -o overrides this.
# output = "split_out/mesh"

# One [[splits]] entry per block to split.  Multiple cuts on the same
# block produce N+1 chunks (N = number of cut indices).
[[splits]]
label = "sub0"
axis  = "k"        # 'i' | 'j' | 'k'
at    = [4]        # vertex indices on the chosen axis; 1 <= idx <= n-1

[[splits]]
label = "background_center"
axis  = "i"
at    = [8, 16, 24]   # 4 chunks
```

Unknown keys raise `ParameterError`. Cut indices are deduplicated and
sorted automatically.

## Adjacency rules

### Internal seams

Each split block is replaced by `len(at) + 1` chunks. Consecutive
chunks share a `SHARED` face along the split axis with **bit-exact
node coincidence** (the chunks alias the same column from the source
block via numpy views — no node duplication).

Example: split `sub0` along `k` at `[4]` yields chunks `sub0__k0` and
`sub0__k1`; `sub0__k0.k_hi` is `SHARED` with `sub0__k1.k_lo`.

### External faces

Faces *not* perpendicular to the split axis (i.e. faces whose name
starts with a different axis letter) appear on **every chunk** —
they're sliced into pieces along the perpendicular axis. The original
block's BC for that face is inherited by every chunk.

Faces perpendicular to the split axis (`<axis>_lo` and `<axis>_hi`)
appear only on the first / last chunk respectively.

### Cross-block partner rewiring (automatic seam co-split)

When a cut subdivides a `SHARED` / `DIRICHLET` seam, the cut is
**automatically propagated to the partner block** (to a fixpoint),
and the partner's face is re-paired to the matching chunk. You no
longer author paired co-splits by hand.

- A cut on a face *perpendicular* to a shared seam induces the
  corresponding cut on the partner's mapped axis (through the seam's
  `AxisMap`, with the flip and the partner's extent). Cut **indices
  on the two sides need not be equal** — they are mapped.
- Propagation runs to a fixpoint, so a cut can induce further cuts on
  a third block, all kept conformal.
- Each side keeps its own BC — a `DIRICHLET` cap↔sub seam stays
  `DIRICHLET` after splitting.

Example: split `sub0` along `k` only. The k-cut auto-propagates to
its j-ring partners (`sub1`, `sub3`) and to both polar caps
(`DIRICHLET`), so every cubed-sphere block ends up k-split conformally
— a single `[[splits]]` entry suffices.

::: tip This replaces a former limitation
Earlier versions rejected an unmatched perpendicular split on a shared
seam ("split all four sub-blocks at the same indices"). That manual
pairing is gone — see
[Design decisions → auto-propagate cuts](/guide/design-decisions#splitting-auto-propagates-cuts-across-seams).
Pass `propagate_seams=False` to the Python API to restore the legacy
behaviour.
:::

## Adjacency JSON

The sidecar at `<stem>.adjacency.json` is a single object with two
top-level keys: `version` (currently 2) and `blocks` (a list).
Each block carries:

```json
{
  "label": "sub0__k0",
  "shape": [20, 9, 5],
  "faces": {
    "i_lo": {
      "bc": "DIRICHLET",
      "partner": {"label": "cap_south", "face": "j_lo"},
      "axis_map": [["i", false], ["k", false]]
    },
    "i_hi": {
      "bc": "DIRICHLET",
      "partner": {"label": "cap_north", "face": "i_lo"},
      "axis_map": [["j", false], ["k", false]]
    },
    "j_lo": {
      "bc": "SHARED",
      "partner": {"label": "sub3__k0", "face": "j_hi"},
      "axis_map": [["i", false], ["k", false]]
    },
    "j_hi": {
      "bc": "SHARED",
      "partner": {"label": "sub1__k0", "face": "j_lo"},
      "axis_map": [["i", false], ["k", false]]
    },
    "k_lo": {"bc": "WALL", "partner": null, "axis_map": null},
    "k_hi": {
      "bc": "SHARED",
      "partner": {"label": "sub0__k1", "face": "k_lo"},
      "axis_map": [["i", false], ["j", false]]
    }
  }
}
```

`bc` is one of `WALL`, `FREE`, `SHARED`, `DIRICHLET`, `PERIODIC`.
`partner` is `null` for `WALL` / `FREE`, otherwise an object with
`label` and `face`.

A `FREE` face may additionally carry an optional `"bc_code": <int>` — an
explicit Xall natural-BC code in `[1, 19]` (e.g. `3` inflow, `9`
extrapolation) that overrides the per-block `free_face_bc` default in
[cc.par](/reference/cc-par). The key is emitted only when set, so older
sidecars round-trip unchanged. See
[Design decisions → per-face BC codes](/guide/design-decisions#per-face-natural-bc-codes).

`axis_map` is `null` for `WALL` / `FREE`, otherwise a 2-element list
of `[partner_axis, flipped]` pairs. The order of the pairs matches
this face's *along-seam axes* — the two axes lying in the face plane
in canonical order: `("j", "k")` for `i_lo`/`i_hi`, `("i", "k")` for
`j_lo`/`j_hi`, `("i", "j")` for `k_lo`/`k_hi`. Each pair says which of
the partner face's axes the corresponding along-axis maps onto, and
whether the indexing direction is reversed (`true`) or preserved
(`false`). See `shore.mesh.face.AxisMap` for the in-Python representation.

> **Versioning.** `version: 1` files (which lacked `axis_map`) are
> rejected by `read_adjacency_json` with a regenerate-the-sidecar
> message. Sidecars are derived artefacts; re-run the topology /
> split that produced them.

Read with the standard library:

```python
import json
adj = json.loads(open("split.adjacency.json").read())
for blk in adj["blocks"]:
    print(blk["label"], blk["shape"])
```

Or use the typed helper:

```python
from shore.io.adjacency import read_adjacency_json
adj = read_adjacency_json("split.adjacency.json")  # validates version
```

## Examples

Split a single primitive box along i at index 6:

```bash
shore primitive box --min 0 0 0 --max 10 5 3 --ni 12 --nj 8 --nk 6 -o box

cat > split.toml <<EOF
version = 1
kind    = "split"
[[splits]]
label = "box"
axis  = "i"
at    = [6]
EOF

shore split box.geo -c split.toml -o box_split
# → box_split_box__i0.geo (7×8×6)
#   box_split_box__i1.geo (6×8×6)
#   box_split.adjacency.json
```

Split a cubed-sphere mesh along k. A single `[[splits]]` entry on `sub0`
is enough — the k-cut auto-propagates across the j-ring to the other
sub-blocks and the DIRICHLET caps:

```bash
shore mesh sphere.stl --ni 30 --nj 40 --nk 13 -o sphere --adjacency-json sphere.adjacency.json

cat > split.toml <<EOF
version = 1
kind    = "split"

[[splits]]
label = "sub0"
axis  = "k"
at    = [6]
EOF

shore split sphere_sub0.geo sphere_sub1.geo sphere_sub2.geo sphere_sub3.geo \
            sphere_cap_north.geo sphere_cap_south.geo \
            -c split.toml -o sphere_split --adjacency sphere.adjacency.json -v
# → all 6 blocks k-split (12 chunks) + adjacency JSON, fully conformal
```

::: warning Pass `--adjacency`
Without the sidecar the input `.geo` blocks are seam-free, so the k-cut
has no seams to propagate across and the caps/sub-blocks are NOT
co-split. Always pass the mesh's `--adjacency` for topology meshes.
:::

## Capabilities & limitations

- **Multi-axis per block in one call.** List one `[[splits]]` entry per
  axis for the same block; it becomes an `ni × nj × nk` chunk grid in a
  single pass. (No need to chain calls as in earlier versions.)
- **Automatic seam co-split.** Cuts propagate across SHARED / DIRICHLET
  seams to the partner block — no manual paired entries.
- **Vertex-aligned cuts only.** No interpolation between existing nodes.
- **Periodic faces are not split.** The cubed-sphere sub-block ↔
  sub-block j-ring is aliased numpy views; splitting it is rejected.

## See also

- [`shore.split`](/reference/api-split) — the Python API
  (`from shore.split import BlockSplit, split_blocks`).
- [`shore.io.adjacency`](/reference/api-adjacency) — reading / writing the
  JSON sidecar (`wire_blocks_from_adjacency`, `propagate_cuts`).
- [`shore balance`](balance) — plans the splits this command applies.
- [Design decisions](/guide/design-decisions) — why splitting works this way.
- [`shore export`](export) — bundling `.geo` files into a `.vtm`
  MultiBlock for ParaView; complementary to splitting.
