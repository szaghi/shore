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
| `-v` / `--verbose` | Print every input block and every output file. |

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

### Cross-block partner rewiring

When the original block had a `SHARED` face whose partner was
*another* block in the input, the partner's face BC is updated to
point at the right new chunk. Three cases:

1. **Partner not split** + face along the split axis → only one chunk
   (first or last) exposes that face; the partner's pointer is
   updated to reference that chunk.
2. **Partner not split** + face perpendicular to the split axis →
   every chunk now exposes the face, but the partner has only one
   face to share with. **Not reconstructible — the splitter raises.**
   To proceed, also split the partner block at matching indices on
   the same axis.
3. **Both split on the same perpendicular axis with identical
   cuts** → chunks pair up 1:1; partner pointers are wired between
   matching chunks (e.g. `sub0__k0.j_hi` ↔ `sub1__k0.j_lo`).

Case (2) is what catches you when splitting the cubed-sphere
sub-blocks: the j-periodic ring shares each `sub_q.j_lo` with
`sub_(q-1).j_hi`. To split any sub-block along i or k you must split
*all four* sub-blocks at the same indices on that axis.

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

Split all cubed-sphere sub-blocks of a sphere mesh along k at the
same index (required because of the j-periodic ring):

```bash
shore mesh sphere.stl --ni 30 --nj 40 --nk 13 -o sphere

cat > split.toml <<EOF
version = 1
kind    = "split"

[[splits]]
label = "sub0"
axis  = "k"
at    = [6]

[[splits]]
label = "sub1"
axis  = "k"
at    = [6]

[[splits]]
label = "sub2"
axis  = "k"
at    = [6]

[[splits]]
label = "sub3"
axis  = "k"
at    = [6]
EOF

shore split sphere_sub0.geo sphere_sub1.geo sphere_sub2.geo sphere_sub3.geo \
            sphere_cap_north.geo sphere_cap_south.geo \
            -c split.toml -o sphere_split -v
# → 8 sub-block chunks + 2 unchanged caps + adjacency JSON
```

## Limitations (v1)

- **Single-axis-per-block per call.** Compose by calling `shore split`
  twice on the result of the first call (e.g. split along i, then
  split each chunk along j).
- **Vertex-aligned cuts only.** No interpolation between existing
  nodes.
- **Cross-block matching is required for perpendicular splits on
  SHARED seams.** The splitter raises a clear error when this is
  violated; add the corresponding `[[splits]]` entry to the partner
  block.

## See also

- [`shore.split.split_blocks`](/reference/api-vtk) — the Python API
  (`from shore.split import BlockSplit, split_blocks`).
- [`shore.io.adjacency`](/reference/api-vtk) — reading / writing the
  JSON sidecar.
- [`shore export`](export) — bundling `.geo` files into a `.vtm`
  MultiBlock for ParaView; complementary to splitting.
