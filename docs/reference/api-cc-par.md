# shore.io.cc_par

Python API for the Xall `cc.par` writer / merger and the
`BlockMeta` per-block metadata dataclass. The CLI wrappers
([`shore cc-par`](./cc-par), [`shore cc-par-merge`](./cc-par-merge))
are thin adapters around these functions.

The file format itself is documented at
[`shore cc-par` — File format reference](./cc-par#file-format-reference).

## Public exports

| Name | Role |
|---|---|
| `write_cc_par` | Translate a SHORE adjacency payload into a `cc.par` file. |
| `merge_cc_par` | Concatenate multiple `cc.par` files for a multi-component Xall run. |
| `BlockMeta` | Per-block metadata dataclass (level / group / priority / comment / free_face_bc). |

## `write_cc_par`

```python
def write_cc_par(
    adjacency: dict[str, Any] | str | Path,
    path: str | Path,
    *,
    grd_filename: str,
    output_basename: str = "cc.",
    save_ghost_cells: bool = True,
    increase_overlap: bool = False,
    extend_internal_wall: bool = False,
    multigrid_levels: tuple[int, int] = (4, 1),
    boundary_layer_thickness: float = -0.1,
    numerical_beach: float = -1.0,
    block_metadata: dict[str, BlockMeta] | None = None,
    boxes: list[Box] | None = None,
) -> Path
```

Translate a SHORE [adjacency payload](./api-adjacency) into an Xall
`cc.par` file. Symmetric pairing of every connection patch is
validated before writing — an asymmetric pairing is fatal.

| Parameter | Description |
|---|---|
| `adjacency` | Either a parsed adjacency dict (as returned by [`read_adjacency_json`](./api-adjacency)) or a path to a `.adjacency.json` file. |
| `path` | Output `.cc.par` path. |
| `grd_filename` | Header line 1: name of the binary `.grd` Xall pairs with this `cc.par`. |
| `output_basename` | Header line 2 (default `"cc."`). |
| `save_ghost_cells`, `increase_overlap`, `extend_internal_wall` | Header booleans. See [Header fields](./cc-par#header-fields). |
| `multigrid_levels` | `(ngr, fgr)` — total levels and finest active. Default `(4, 1)`. |
| `boundary_layer_thickness`, `numerical_beach` | Header reals. Negative disables. |
| `block_metadata` | Optional per-block override `{label: BlockMeta}`. Missing labels get the default `BlockMeta()`. |
| `boxes` | Optional list of [`Box`](./api-boxes-io) instances; emitted as the boxes section (type-9 hexahedra). Each box's `associated_block` label is resolved to a 1-based index into the blocks list. |

**Returns** the actual `Path` written.

**Raises** `MeshInputError` for missing files, malformed adjacency,
unknown block label in `boxes`, or a symmetry-pairing violation.

## `BlockMeta`

```python
@dataclass
class BlockMeta:
    level:        int = 1
    group:        int = 0
    priority:     int = 1
    comment:      str = ""        # defaults to block label
    free_face_bc: int = 40        # offgen
```

One blocks-section line's metadata, plus the BC code that the
writer emits for that block's `FREE` faces.

| Field | Description |
|---|---|
| `level` | Overset priority among overlapping bodies. Higher overrides lower when chimera holes overlap. |
| `group` | Group index — bundle blocks belonging to the same body. |
| `priority` | Within-group priority for donor search tie-breaking. |
| `comment` | Trailing `! comment` text. Defaults to the block label. |
| `free_face_bc` | BC code emitted for `FREE` faces of this block. **Default `40`** (`offgen`, chimera cell-center) — appropriate for body-fitted / buffer meshes. Use a negative natural-BC code (`-9` extrapolation, `-3` inflow, `-4` inoutflow) on far-field background blocks. |

```python
from shore.io.cc_par import BlockMeta, write_cc_par

meta = {
    "sub0": BlockMeta(free_face_bc=40),
    "sub1": BlockMeta(free_face_bc=40),
    "sub2": BlockMeta(free_face_bc=40),
    "sub3": BlockMeta(free_face_bc=40),
    "cap_north": BlockMeta(level=2, free_face_bc=-9),
    "cap_south": BlockMeta(level=2, free_face_bc=-9),
}
write_cc_par(
    "wall.adjacency.json",
    "wall.cc.par",
    grd_filename="wall.grd",
    block_metadata=meta,
)
```

## `merge_cc_par`

```python
def merge_cc_par(
    inputs: list[str | Path],
    out: str | Path,
    *,
    grd_filename: str,
) -> Path
```

Concatenate two-or-more `cc.par` files into one. The merged file's
blocks list is the concatenation of the inputs; patch `block`
fields, `family` pointers (partner patch indices), and box
`associated_block` indices are shifted by the cumulative block /
patch counts of the preceding inputs so all references stay
internally consistent.

| Parameter | Description |
|---|---|
| `inputs` | Two-or-more input `cc.par` paths in the desired output order. |
| `out` | Output path. |
| `grd_filename` | Header line 1 of the merged file (overrides the per-input value). |

**Returns** the actual `Path` written.

**Raises** `ValueError` for fewer than 2 inputs; `FileNotFoundError`
for a missing input. Header field disagreement between inputs is
emitted as `UserWarning` but does not block the merge.

::: warning Pair with grd-merge
The merged `cc.par`'s patch `block` indices and box
`associated_block` indices are 1-based positions into the merged
`.grd`'s block list. Always invoke `merge_cc_par` and
[`merge_grd`](./api-grd) **as a pair, with the same input ordering** —
otherwise the merged files don't agree.
:::

```python
from shore.io.cc_par import merge_cc_par
from shore.io.grd import merge_grd

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

## See also

- [`shore cc-par` CLI](./cc-par) — wrapper for `write_cc_par` and
  full file format reference.
- [`shore cc-par-merge` CLI](./cc-par-merge) — wrapper for
  `merge_cc_par`.
- [`shore.io.adjacency`](./api-adjacency) — read / write the input
  payload.
- [`shore.io.grd`](./api-grd) — companion binary mesh writer /
  merger.
- [Marker boxes](./boxes) — the `boxes` argument is documented in
  detail there.
- [`shore.io.boxes`](./api-boxes-io) — the `Box` dataclass.
