# `shore grd`

Pack one or more `.geo` files into a single Xall binary `.grd` mesh file.

```bash
shore grd GEO_FILES... --output MESH.grd
```

## Motivation

Xall consumes a **single** binary `.grd` file per simulation
component, containing all blocks in one Fortran unformatted stream. SHORE
produces one `.geo` file per block (one file per topology sub-block, or one
file per chunk after `shore split`). `shore grd` is the bridge: it reads the
per-block `.geo` files in the order you specify and packs them into the binary
format Xall expects, ready to pair with a `cc.par`.

The block ordering in the `.grd` must match the ordering used when generating
the companion `cc.par` — both derive their block numbering from the same
ordered list.

## Binary format reference

`.grd` is **Fortran unformatted sequential binary**: every `WRITE` statement
in the Fortran writer wraps its payload between two identical 4-byte
record-length markers (native-endian `int32`). Xall reads it with matching
`READ` calls that rely on the markers to skip fields.

Record layout (block count `nb`, block index `n = 1..nb`):

```
Record 1         : nb                    int32     number of blocks
Records 2..nb+1  : ni_n, nj_n, nk_n     3x int32  per-block cell counts (= node count - 1)
For each block n:
  Record nb+1+3(n-1)+1  : X flat         nb_n float64   X coordinates, i-fastest (Fortran column-major)
  Record nb+1+3(n-1)+2  : Y flat         nb_n float64
  Record nb+1+3(n-1)+3  : Z flat         nb_n float64
```

where `nb_n = (ni_n+1) * (nj_n+1) * (nk_n+1)` nodes.

### Cell counts vs. node counts

The dimension records store **cell counts** (`ni - 1`, `nj - 1`, `nk - 1`),
not node counts. This matches the convention in the reference Fortran
pre-processor (`geogrd.f90`): it reads `(ni, nj, nk)` from a `.geo` file and
immediately subtracts 1 before writing. SHORE `.geo` files store node counts;
`shore grd` performs the subtraction automatically.

### Coordinate ordering

Coordinates are written as flat `float64` arrays in **Fortran column-major**
order (i-index varies fastest). This matches the `order="F"` used in `write_geo`
and the layout Xall's Fortran reader expects.

## Arguments and options

| Argument / Option | Description |
|---|---|
| `GEO_FILES...` | One or more `.geo` input files. Block order in the `.grd` matches the order given here. |
| `-o` / `--output PATH` | Output `.grd` file. The `.grd` extension is added automatically if absent. |

## Examples

```bash
# Single block
shore grd sphere.geo --output sphere.grd

# 6-block cubed-sphere near-body mesh
shore grd wall_sub0.geo wall_sub1.geo wall_sub2.geo wall_sub3.geo \
          wall_cap_north.geo wall_cap_south.geo \
          --output wall.grd

# 5-block background mesh after a split
shore grd background_center.geo \
          background_sub_e.geo background_sub_n.geo \
          background_sub_w.geo background_sub_s.geo \
          --output background.grd
```

The block order must agree with the `cc.par` generated from the same
`.adjacency.json` sidecar. Both `shore grd` and `shore cc-par` derive their
block numbering from the ordered input; keep the lists consistent.

## Python API

```python
from shore.io.grd import write_grd

# Pass an ordered list of .geo file paths
out = write_grd("wall.grd", [
    "wall_sub0.geo",
    "wall_sub1.geo",
    "wall_sub2.geo",
    "wall_sub3.geo",
    "wall_cap_north.geo",
    "wall_cap_south.geo",
])
```

`write_grd` accepts `str` or `pathlib.Path` for both arguments. It raises
`ValueError` if the path list is empty.

## Typical workflow

```bash
# 1. Build the near-body mesh
shore mesh sphere.stl --ni 30 --nj 40 --nk 13 -o wall

# 2. (Optional) split blocks for parallel decomposition
shore split wall_sub0.geo ... -c split.toml -o wall_split

# 3. Write the binary grid
shore grd wall_sub0.geo wall_sub1.geo wall_sub2.geo wall_sub3.geo \
          wall_cap_north.geo wall_cap_south.geo \
          --output wall.grd

# 4. Write the chimera descriptor (same block order)
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd
```

## Merging multiple `.grd` files (`shore grd-merge`)

When an Xall job combines several independent Chimera components — for
example a near-body wall mesh plus a far-field background — each
component is generated as its own `.grd` + `cc.par` pair. Running them
through Xall as one job needs both pairs concatenated into a single
`assembly.grd` and a single `assembly.cc.par` whose patch indices are
shifted to refer to the merged block list.

`shore grd-merge` does the file half:

```bash
shore grd-merge background.grd wall.grd -o assembly.grd
```

The merged file lists every block from `inputs[0]` first, then
`inputs[1]`, etc. The header block count is rewritten to the sum;
shape and coordinate records are passed through byte-for-byte (no
node values are touched). Order is the caller's choice — it matters
because Xall's overset solver uses block ordering to break level /
priority ties — but it **must** match the order used for the
companion [`shore cc-par-merge`](./cc-par#merging-multiple-cc-par-files-shore-cc-par-merge).

### Python API

```python
from shore.io.grd import merge_grd

merge_grd(
    ["background.grd", "wall.grd"],
    "assembly.grd",
)
```

Raises `ValueError` for fewer than 2 inputs or a malformed `.grd` file.

## See also

- [`shore cc-par`](./cc-par) — generate the companion Xall chimera descriptor.
- [`shore cc-par-merge`](./cc-par#merging-multiple-cc-par-files-shore-cc-par-merge) —
  the companion merge for chimera descriptors. Always invoke `grd-merge`
  and `cc-par-merge` as a pair, with the same input ordering.
- [`shore split`](./split) — subdivide blocks and produce the `.adjacency.json` sidecar.
- [`shore export`](./export) — export `.geo` files to VTK XML (ParaView), not Xall format.
