---
title: shore.io.cc_par
outline: [2, 3]
---

# `shore.io.cc_par`

Xall `cc.par` writer.

`cc.par` is the chimera pre-processor input file consumed by the
Xall / Xall tooling.  It pairs with a `.grd` mesh file to declare
per-block multigrid level / group / priority and to enumerate every
patch (face boundary) on every block.  Some patches are simple BCs
(walls, inflow); others are *connections* between adjacent blocks.

This module is a one-way bridge: it consumes the SHORE
`.adjacency.json` sidecar (see :mod:`shore.io.adjacency`) and emits a
fully populated `cc.par` ASCII file ready for Xall ingestion.  It
does **not** read `cc.par` files — the Xall side has its own
`cc_par_object` reader.

Orientation codes are derived **from the adjacency JSON's
`axis_map` field** — there are no per-topology lookup tables in
this writer.  Any topology that produces a fully-wired adjacency
(every SHARED/DIRICHLET face carries `axis_map` and `partner`)
generates correct `cc.par` output.

Format reference
----------------
`cc.par` has 5 sections:

1. **Header** (10 lines): grd filename, output basename, 3 booleans,
   blank, 2 multigrid levels, blank, boundary-layer thickness,
   numerical-beach width, blank.
2. **Blocks**: `N` followed by `N` lines of
   `level group priority ! comment`.
3. **Patches**: `N` followed by `N` lines of
   `block face BC family Imin Imax Jmin Jmax Kmin Kmax ! patch_index`.
   `face` is 1..6 for Imin/Imax/Jmin/Jmax/Kmin/Kmax.  Extents are in
   **cell counts**, i.e. `n_nodes - 1` per axis.  `BC > 99` with
   `family > 0` denotes a connection; `family` is the partner
   patch's 1-based index, and the BC value's three digits encode the
   partner-axis orientation (see :func:`_orientation_code`).
4. **Edges**: always `0` for SHORE-generated meshes.
5. **Boxes**: hexahedral marker volumes that flag cells as
   "internal wall" for the chimera hole-cutter.  Optional; pass a
   list of :class:`shore.io.boxes.Box` via `boxes=...` to emit
   them.  Defaults to `0 ! boxes number`.

Reference implementation: `overset_exploded_cc_par_object.F90` in
`Xnavis-projects/overset-exploded` (the Fortran writer this Python
output is meant to round-trip with).

BC code table
-------------
Mapping from SHORE's `FaceBC` enum to `cc.par` BC codes:

- `WALL`       → `1` (`parete`, viscous wall) + `family = 0`.
- `FREE`       → caller-chosen via :class:`BlockMeta.free_face_bc`,
  defaulting to `40` (`offgen`, chimera face cell-center).  Body-
  fitted and buffer meshes use `40`; far-field background patches
  use a natural BC, which MUST be a positive code in `[1, nccnat=19]`
  (e.g. `9` for free/far-field extrapolation, `3`/`4` for
  inflow/outflow).  `calface` stores the code verbatim into `icc`
  and `calspig2` / `calvert` reject any face `icc < 1` ("Errore
  in calspig: icc sulle facce deve essere &gt; 0"), so a *negative* code
  such as `-9` is invalid and fails where two natural faces meet.
- `SHARED` / `DIRICHLET` → connection orientation code (`>= 100`)
  + `family = partner patch index`.  The 3-digit orientation code
  `lᵢ lⱼ lₖ` is derived from this face's `axis_map` plus the
  partner face's name (see :func:`_orientation_code_from_axis_map`).
  The normal digit is odd when the two faces are on opposite sides of
  their blocks and even when on the same side, as the overset
  pre-processor reads it (`metrix.f90` `normind`, ADR 0036).  A face
  whose sidecar entry carries `patches` (shared with several partners)
  gets one patch per partner over that partner's sub-range.
- `PERIODIC`   → `cc.par` has no first-class periodic BC; emits
  `40` with a `UserWarning`.
- `SYMMETRY`   → `2` (natural-BC index of the symmetry plane; the
  preprocessor maps it to its internal `BC_NATURAL_SIMMETRY = -2`)
  + `family = 0`.

Limitations
-----------
- Boxes: type-9 hexahedron only on emit.  The Xall reader also
  accepts type-1 (parallelepiped), type-2 (cylinder, deprecated),
  and type-3 (coaxial cylinders) shorthand forms; the merger
  passes those through verbatim from input files but the writer
  never produces them.
- No `edges` section — unsupported by SHORE-style adjacency.

## `BlockMeta` <Badge type="tip" text="dataclass" /> {#BlockMeta}

```python
BlockMeta(
    level: int = 1,
    group: int = 0,
    priority: int = 1,
    comment: str = '',
    free_face_bc: int = _BC_OFFGEN,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/cc_par.py#L165">source</a>

Per-block metadata for a `cc.par` blocks-section line.

**Parameters**

- **`level`** — Overset level (priority among overlapping bodies).  Default 1.
- **`group`** — Group index.  Default 0.
- **`priority`** — Within-group priority.  Default 1.
- **`comment`** — Trailing `! comment` for the line.  Defaults to the block label.
- **`free_face_bc`** — BC code emitted for `FREE` faces of this block.  Default `40` (`offgen`, chimera-cell-center) — appropriate for body-fitted / buffer meshes.  On far-field background blocks use a positive natural-BC code in `[1, nccnat=19]` (e.g. `9` for free/far-field extrapolation, `3`/`4` for inflow/outflow).  Negative codes (`-9`, `-3`) are invalid: overset stores the code verbatim and rejects any face `icc < 1` at edge/vertex cells.

## `write_cc_par` <Badge type="info" text="function" /> {#write_cc_par}

```python
write_cc_par(
    adjacency: dict[str, Any] | str | Path,
    path: str | Path,
    *,
    grd_filename: str,
    output_basename: str = 'cc.',
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
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/cc_par.py#L237">source</a>

Write a `cc.par` file from a SHORE adjacency payload.

**Parameters**

- **`adjacency`** — Either a parsed adjacency dict (as produced by :func:`shore.io.adjacency.write_adjacency_json` and read back via :func:`shore.io.adjacency.read_adjacency_json`) or a path to a `.adjacency.json` file on disk.
- **`path`** — Output `cc.par` path.
- **`grd_filename`** — Name of the `.grd` mesh file Xall should pair with this `cc.par` (header line 1).
- **`output_basename`** — Output basename Xall will use for its derived files (header line 2).  Default `"cc."`.
- **`save_ghost_cells, increase_overlap, extend_internal_wall`** — Header booleans (lines 3-5).
- **`multigrid_levels`** — Two-element tuple `(N_levels, finest_level_active)` (line 7).
- **`boundary_layer_thickness, numerical_beach`** — Header reals (lines 9-10).
- **`block_metadata`** — Optional `{label: BlockMeta}` dict.  Any block not listed gets a default `BlockMeta()`.
- **`boxes`** — Optional list of :class:`shore.io.boxes.Box` instances written as the `cc.par` boxes section (after the edges section). Each box is emitted as type-9 (hexahedron); its `associated_block` label is resolved to a 1-based block index against the adjacency block list.  Defaults to `None`, which writes `0 ! boxes number`.

**Returns**

- `Path` — The actual path written.

**Raises**

- `ParameterError` — If a box references a block label not present in the adjacency.

## `merge_cc_par` <Badge type="info" text="function" /> {#merge_cc_par}

```python
merge_cc_par(inputs: list[str | Path], out: str | Path, *, grd_filename: str) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/cc_par.py#L801">source</a>

Concatenate multiple `cc.par` files into one for an Xall multi-mesh run.

Used together with :func:`shore.io.grd.merge_grd` when an Xall job
combines several independent Chimera components (e.g. wall +
background): the merged `cc.par` lists every block and patch from
*inputs[0]* first, then *inputs[1]*, etc.  All block indices and
patch `family` references in inputs[1:] are shifted by the
cumulative block / patch counts of the preceding inputs so the
merged file is internally consistent.

**Parameters**

- **`inputs`** — Two-or-more input `cc.par` paths, in the desired output order. The same order **must** be used when calling :func:`shore.io.grd.merge_grd` — `cc.par` patch block-indices are 1-based positions into the `.grd` block list.
- **`out`** — Output path; `.par` extension is **not** auto-appended (cc.par filenames are conventional, not extension-driven).
- **`grd_filename`** — Name of the merged `.grd` file, written into header line 1 (overriding whatever filename inputs[0]'s header carried).

**Returns**

- `Path` — The path actually written.

**Notes**

Header reconciliation: every header field other than `grd_filename`
is taken from inputs[0] verbatim (output basename, ghost-cells flag,
multigrid levels, BL thickness, numerical beach).  When inputs disagree
on those fields a `UserWarning` is emitted and inputs[0]'s values
are used; surface drift loudly without blocking.

Block-comment prefixing: each block's comment is prefixed with the
input file's stem (e.g. `wall.sub0`, `background.center`) so
blocks from different topologies are distinguishable in the merged
file.  The cc.par format treats comments as advisory only — Xall
indexes blocks by their 1-based position.

**Raises**

- `ValueError` — If `inputs` has fewer than 2 entries or any input is malformed.
