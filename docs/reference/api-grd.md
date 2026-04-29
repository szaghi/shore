# shore.io.grd

Pack and merge Xall binary `.grd` mesh files. SHORE produces one
`.geo` per block; `write_grd` concatenates them into the single
Fortran-unformatted binary that the Xall solver consumes, and
`merge_grd` concatenates multiple `.grd` files for multi-component
Xall runs.

For the binary format and CLI, see [`shore grd`](./grd) and
[`shore grd-merge`](./grd-merge).

## `write_grd`

```python
def write_grd(
    path: str | Path,
    geo_paths: list[str | Path],
) -> Path
```

Pack one or more `.geo` files into a single Xall `.grd`. The input
list determines both the block count and the block ordering — which
must match the ordering used to generate the companion
[`cc.par`](./api-cc-par).

| Parameter | Description |
|---|---|
| `path` | Output file path. `.grd` extension appended if absent. |
| `geo_paths` | Ordered list of `.geo` paths, one per block. Must not be empty. |

**Returns** the actual `Path` written.

The dimension records store **cell counts** (`ni-1`, `nj-1`, `nk-1`),
not node counts. Coordinates are written as flat `float64` arrays
in **Fortran column-major** order (i-fastest), matching the format
the Xall reader expects. Records use Fortran unformatted sequential
binary format (4-byte length markers, native endianness).

## `merge_grd`

```python
def merge_grd(inputs: list[str | Path], out: str | Path) -> Path
```

Concatenate two-or-more `.grd` files into one. The merged file
lists every block from `inputs[0]` first, then `inputs[1]`, etc.;
the header block count is rewritten to the sum. Shape and
coordinate records are passed through byte-for-byte — no node
values are touched.

| Parameter | Description |
|---|---|
| `inputs` | Two-or-more input `.grd` paths, in the desired output order. |
| `out` | Output path. `.grd` extension appended if absent. |

**Returns** the actual `Path` written.

**Raises** `ValueError` for fewer than 2 inputs or a malformed
input file.

### Coordinated invocation

Invoke `merge_grd` and [`shore.io.cc_par.merge_cc_par`](./api-cc-par)
**as a pair, with the same input ordering**. The companion
`merge_cc_par` shifts patch block-indices and `family` pointers
(and box `associated_block` indices) to match the merged block
list; without the matching `cc.par` the merged `.grd` is
unusable.

```python
from shore.io.grd import merge_grd
from shore.io.cc_par import merge_cc_par

merge_grd(
    ["background.grd", "wall.grd"],
    "assembly.grd",
)
merge_cc_par(
    ["background.cc.par", "wall.cc.par"],
    "assembly.cc.par",
    grd_filename="assembly.grd",
)
```

## Example

```python
from shore.io.grd import write_grd

write_grd(
    "wall.grd",
    [
        "wall_sub0.geo", "wall_sub1.geo",
        "wall_sub2.geo", "wall_sub3.geo",
        "wall_cap_north.geo", "wall_cap_south.geo",
    ],
)
```

## See also

- [`shore grd` CLI](./grd) — the command-line wrapper for
  `write_grd`.
- [`shore grd-merge` CLI](./grd-merge) — wrapper for `merge_grd`.
- [`shore.io.cc_par`](./api-cc-par) — the companion ASCII chimera
  descriptor.
- [`shore.io.geo`](./api-geo) — `.geo` reader / writer.
