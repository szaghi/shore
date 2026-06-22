# `shore cc-par`

Convert a SHORE `.adjacency.json` sidecar to an Xall `cc.par` chimera
pre-processor input file.

```bash
shore cc-par INPUT.adjacency.json OUTPUT.par --grd MESH.grd [OPTIONS]
```

## Motivation

Xall is a Fortran CFD solver that consumes **two**
files per simulation: a binary mesh (`.grd`) and an ASCII chimera
descriptor (`cc.par`). The `cc.par` file lists every block's overset
metadata (level, group, priority) and every face *patch* — boundary
condition, IJK extents, and partner reference for inter-block
connections.

SHORE produces the structured grid (packed into a binary `.grd` via
`shore grd`) and a JSON adjacency sidecar that captures the multi-block
connectivity. `shore cc-par` is the bridge: it translates the SHORE
adjacency JSON into the Xall ASCII format, ready to feed into the
Xall pre-processor.

## File format reference

`cc.par` is plain ASCII, free-format, parsed by Fortran list-directed
`read(npar,*)` calls. **Order matters** — the parser reads sections
sequentially. **Comments after `!`** are ignored by the reader but
emitted by SHORE / overset-exploded for human readability.

The canonical readers / writers are:

- `Xnavis-projects/overset/src/readpar.f90` — Xall reader (the
  authoritative format definition).
- `Xnavis-projects/overset-exploded/src/lib/overset_exploded_cc_par_object.F90`
  — modular reader/writer used by `overset-exploded` for
  load-balancing.
- `Xnavis-projects/overset/src/calcorn.f90` lines 119–125, 392–411 —
  the orientation-digit decoder (`tadia` / `nadia` subroutines).
- `Xnavis-projects/overset/src/mod/modpar.f90` — BC code constants.

The file has **5 sections**, in this order:

```
┌─ Header ──────────────────────┐  10 lines
│ '<grd-filename>'              │  line 1
│ '<output-basename>'           │  line 2
│ <save-ghost-cells>            │  line 3 (.true. / .false.)
│ <increase-overlap>            │  line 4 (.true. / .false.)
│ <extend-internal-wall>        │  line 5 (.true. / .false.)
│ <blank>                       │  line 6
│ <ngr> <fgr>                   │  line 7 (multigrid levels)
│ <blank>                       │  line 8
│ <boundary-layer-thickness>    │  line 9 (real)
│ <numerical-beach>             │  line 10 (real)
│ <blank>                       │  line 11
└───────────────────────────────┘
┌─ Blocks ──────────────────────┐
│ <N>                           │  number of blocks
│ <blank>                       │
│ <level> <group> <priority>  ! <comment>     ×N
│ <blank>                       │
└───────────────────────────────┘
┌─ Patches ─────────────────────┐
│ <M>                           │  number of patches
│ <blank>                       │
│ <block> <face> <BC> <family> <Imin> <Imax> <Jmin> <Jmax> <Kmin> <Kmax> ! <patch_index>     ×M
│ <blank>                       │
└───────────────────────────────┘
┌─ Edges ───────────────────────┐
│ <K>                           │  number of wall-edge entries (always 0 for SHORE)
│ <blank>                       │
│ ... (K lines if K > 0)        │
└───────────────────────────────┘
┌─ Boxes ───────────────────────┐
│ <B>                           │  number of hole-cutting boxes (always 0 for SHORE)
│ <blank>                       │
│ ... (B records of 1+8 lines)  │
└───────────────────────────────┘
```

### Header fields

| Line | Field | Type | Meaning |
|---|---|---|---|
| 1 | grd filename | quoted string | Path/name of the binary `.grd` mesh file Xall pairs with this `cc.par`. |
| 2 | output basename | quoted string | Basename Xall uses for derived output files (e.g. `cc.`). |
| 3 | save ghost cells | `.true.` / `.false.` | Write block files with ghost-cell frames. |
| 4 | increase overlap | `.true.` / `.false.` | Enlarge inter-block overlap during chimera donor search. |
| 5 | extend internal wall | `.true.` / `.false.` | Extend wall BC into chimera-boundary cells. |
| 7 | `ngr fgr` | two integers | Multigrid: `ngr` = total levels, `fgr` = finest active level (typically `4 1`). |
| 9 | boundary layer thickness | real | Estimate of physical boundary-layer thickness; affects internal-wall and overlap heuristics. Negative disables. |
| 10 | numerical beach | real | Width of numerical-beach damping zone. Negative disables. |

Lines 6, 8, 11 are intentionally blank — parser advances past them
with bare `read(npar,*)` calls.

### Blocks section

```
N ! blocks number
                  ← blank line
<level> <group> <priority> ! <comment>     ← N lines, one per block
                  ← blank line
```

| Field | Meaning |
|---|---|
| `level` | Overset *priority among bodies*. Higher `level` overrides lower when chimera holes overlap. |
| `group` | Group index — used to bundle blocks belonging to the same body. |
| `priority` | Within-group priority for donor search tie-breaking. |
| `comment` | Free-form text. Conventionally the block name or a short label. |

Block index is **implicit**: the first listed block is block 1, etc.
This index is referenced by patches in the next section.

### Patches section

The heart of the format. **Every face of every block is one patch.**
A block with no split patches has exactly 6 lines (one per face);
faces split into multiple sub-patches contribute more.

```
M ! patches number
                              ← blank line
<block> <face> <BC> <family> <Imin> <Imax> <Jmin> <Jmax> <Kmin> <Kmax> ! <patch_index>     ← M lines
                              ← blank line
```

**10 numeric columns + 1 trailing comment**:

| Col | Field | Type | Meaning |
|---|---|---|---|
| 1 | `block` | int | 1-based index into the blocks section. |
| 2 | `face` | int 1–6 | Which face: `1`=Imin, `2`=Imax, `3`=Jmin, `4`=Jmax, `5`=Kmin, `6`=Kmax. |
| 3 | `BC` | int | Boundary-condition code (see below). |
| 4 | `family` | int | For walls: family index for grouping (priority). For connections (`BC ≥ 100`): partner patch's 1-based `patch_index`. `0` when neither applies. |
| 5–6 | `Imin Imax` | int | Patch i-extent in **cell** counts (not nodes). |
| 7–8 | `Jmin Jmax` | int | Patch j-extent in cells. |
| 9–10 | `Kmin Kmax` | int | Patch k-extent in cells. |
| trailing | `! patch_index` | int comment | 1-based patch index (must equal the line's position in the section; SHORE writes it for readability). |

For a face perpendicular to axis A, the A-extents collapse (`Amin == Amax = 0` for the lo face, `n_A_cells` for the hi face) and the other two axes span the full block face.

### BC code constants

From `Xnavis-projects/overset/src/mod/modpar.f90` and
`overset-exploded_block_object.F90`:

| Code | Symbol | Meaning |
|---|---|---|
| `0` | `stdflag` | Active interior cell (not a BC). |
| `1` | `parete` | Adiabatic viscous wall. |
| `11` | `inawal` | Inactive wall, adiabatic. |
| `12` | `isowal` | Isothermal wall. |
| `13` | `isoina` | Inactive wall, isothermal. |
| `-1` | `BC_NATURAL_WALL` | Natural-BC viscous wall. |
| `-2` | symmetry | Symmetry plane. |
| `-3` | inflow | Inflow. |
| `-4` | inoutflow | In/outflow. |
| `-5` | assigned inflow | Assigned-inflow. |
| `-6` | assigned pressure | Assigned-pressure. |
| `-7` | assigned normal velocity | |
| `-8` | assigned Riemann | |
| `-9` | extrapolated | Extrapolation (typical far-field). |
| `-10` | moving wall | |
| `-11` | inactive wall | |
| `20` | `offchi` | Offset for chimera **face-center** recipients. `21..26` = chimera face at I0/IN/J0/JN/K0/KN. |
| `27` | `innerx` | Chimera interior point. |
| `28` | `innerw` | Chimera internal wall. |
| `29` | `innerb` | Numerical-beach (spiaggia) cell. |
| `40` | `offgen` | Offset for chimera **cell-center** recipients. `41..46` = chimera face at I0/IN/J0/JN/K0/KN. |
| `60` | `offbiu` | Offset for one-to-one adjacency (recipient side). `61..66` = adjacent at I0/IN/J0/JN/K0/KN. |
| `70` | `estrap` | Extrapolation for moving bodies. |
| `80` | `xedge` | Edge cell. |
| **`≥ 100`** | **connection** | **Resolved one-to-one connection: orientation-coded, requires `family > 0`.** |

Codes ≥ 100 are not constants — they are *computed* from the
orientation between two adjacent blocks (see next section). The
parser distinguishes connections from BCs by `BC > 99 AND family > 0`
(`readpar.f90:147`). Patches with `BC ≥ 100` must come in **symmetric
pairs**: if patch *p* names *q* as its family, then patch *q* must
name *p* as its family (`readpar.f90:219-228`); the symmetry check is
fatal.

### Edges section

```
0 ! edges number
                ← blank line
```

Wall-edge entries flag concave wall corners for boundary-layer
treatment. SHORE produces no concave wall geometry, so this section
is always `0` followed by a blank line.

### Boxes section

```
N ! boxes number
                       ← blank line
        9       <ibl>  ! box k (type 9, hexahedron)
        v1.x  v1.y  v1.z
        v2.x  v2.y  v2.z
        ...
        v8.x  v8.y  v8.z
                       ← blank line
```

Boxes are **hole-cutting volumes** — closed regions that mark cells
inside them as "internal wall". Conceptually each box contains a
solid body and itself fits inside the body's near-mesh: the union
of one body's boxes lies in the thin shell between the STL surface
and the outer extent of the body's near-mesh, so that a *background*
mesh's cells overlapping the body land inside at least one box and
get hole-cut.

SHORE emits **type 9 (hexahedron)** only — eight vertex records
per box, ordered by the canonical (i, j, k) bit layout described in
`shore.io.boxes` (and reproduced in [Marker boxes](#marker-boxes)
below). Other types (parallelepiped / cylinder shorthands the Xall
reader supports) are passed through unchanged by the merger but
never emitted by the writer.

`<ibl>` is the 1-based block index this box is associated with;
the box inherits that block's `level / group / priority` for the
chimera precedence rules in `check_parete_interna`
(`overset/src/scatole.f90`). See the [Marker boxes](#marker-boxes)
section below for the SHORE Python API and the `--boxes` CLI flag.

## Adjacency map: how connections work

Every pair of conformal adjacent blocks contributes **two** patch
lines — one on each side of the seam. Both reference each other via
the `family` field, and both carry orientation codes that describe
how their `(i, j, k)` index spaces map across the seam.

### Worked example: same-axis seam

Two cubes `A` and `B` sharing a face: `A`'s `j_hi` (face 4) is glued
to `B`'s `j_lo` (face 3). Both cubes are `8×8×8` cells, oriented the
same way (no flip).

Patches:

```
1   4   135   2    0   8   8   8   0   8  ! 1     A.j_hi -> B.j_lo
2   3   145   1    0   8   0   0   0   8  ! 2     B.j_lo -> A.j_hi
```

- `135` (A side): A's i → B's I-min increasing (along, l=1); A's j
  (perpendicular at Jmax) → B's J-min increasing inward (l=3); A's k
  → B's K-min increasing (l=5).
- `145` (B side): B's i → A's I-min increasing (along, l=1); B's j
  (perpendicular at Jmin) → A's J-max decreasing inward (l=4); B's k
  → A's K-min increasing (l=5).
- The two codes differ even though the seam is geometrically symmetric:
  each side reports its perpendicular partner from its own frame, and
  the partner is at Jmin on one side and Jmax on the other. Both codes
  describe the same seam — see the [orientation digits](#orientation-digits-in-detail)
  table for the digit semantics.
- `family = 2` and `family = 1` cross-reference the two patches.
- A's j-extents are `8 8` (collapsed at j-max); B's j-extents are
  `0 0` (collapsed at j-min).

### Worked example: cross-axis seam (cubed-sphere cap-equator)

Sub-block `sub0` is one of the 4 equatorial blocks in the cubed-sphere
topology; its north pole row connects to `cap_north`. Geometrically,
**sub0's j-axis is preserved into cap_north's j-axis**, but **sub0's i
maps to cap_north's i** (perpendicular axis swap from Imax→Imin):

```
1   2   135  25   8   8   0   8   0   4  ! 2     sub0.i_hi -> cap_north.i_lo
5   1   235   2   0   0   0   8   0   4  ! 25    cap_north.i_lo -> sub0.i_hi
```

- `135` (sub side): sub's i → cap's I-min increasing (perpendicular,
  l=1); sub's j → cap's J-min increasing (preserved, l=3); sub's k →
  cap's K-min increasing (preserved, l=5).
- `235` (cap side): cap's i → sub's I-max decreasing (perpendicular
  flipped, l=2); cap's j → sub's J-min increasing (preserved, l=3);
  cap's k → sub's K-min increasing (preserved, l=5).

The two BCs differ even though the geometric seam is the same — each
side describes the mapping *from its own axes to its partner's*. The
asymmetry comes from "Imax meets Imin": the perpendicular-axis digit
on the sub side records "partner Imin going inward", while on the
cap side it records "partner Imax going inward".

### Orientation digits in detail

A connection BC code has three digits `li lj lk`:

```
BC = 100·li + 10·lj + lk
```

Each digit `l ∈ {1..6}` says how *this side's* axis maps onto a
*partner* axis:

| `l` | Partner axis & start | Direction |
|---|---|---|
| `1` | I, starting at Imin | increasing inward |
| `2` | I, starting at Imax | decreasing inward (flipped) |
| `3` | J, starting at Jmin | increasing inward |
| `4` | J, starting at Jmax | decreasing inward (flipped) |
| `5` | K, starting at Kmin | increasing inward |
| `6` | K, starting at Kmax | decreasing inward (flipped) |

Even-valued digits (`2, 4, 6`) flag a **flipped** axis correspondence;
odd-valued digits (`1, 3, 5`) preserve direction. This is implemented
in `tadia` / `nadia` subroutines (`calcorn.f90:392-411`):

```fortran-free-form
! tadia: along-seam axis index
if      (l.eq.1) then; ia = nla(1) + n          ! Imin, increasing
else if (l.eq.2) then; ia = nla(2) + 1 - n      ! Imax, decreasing
else if (l.eq.3) then; ja = nla(3) + n          ! Jmin, increasing
else if (l.eq.4) then; ja = nla(4) + 1 - n      ! Jmax, decreasing
else if (l.eq.5) then; ka = nla(5) + n          ! Kmin, increasing
else if (l.eq.6) then; ka = nla(6) + 1 - n      ! Kmax, decreasing
```

For a same-axis seam with no flips, the side whose perpendicular face
is at min (e.g. `j_lo`) emits `135` and the side at max (e.g. `j_hi`)
emits `145`. The codes are *not* symmetric: each side describes the
mapping in its own frame, and the perpendicular digit reflects which
boundary of the partner the seam lies on.

Decoding examples:

| Code | li | lj | lk | Meaning |
|---|---|---|---|---|
| `135` | 1 | 3 | 5 | i→Imin↑, j→Jmin↑, k→Kmin↑ — *no flips, perpendicular goes Imin* |
| `145` | 1 | 4 | 5 | i→Imin↑, j→Jmax↓, k→Kmin↑ — *perpendicular Jmax-flipped (mirror of 135)* |
| `235` | 2 | 3 | 5 | i→Imax↓, j→Jmin↑, k→Kmin↑ — *perpendicular Imax-flipped* |
| `415` | 4 | 1 | 5 | i→Jmax↓, j→Imin↑, k→Kmin↑ — *axis swap (i↔j) + flip on i* |
| `245` | 2 | 4 | 5 | i→Imax↓, j→Jmax↓, k→Kmin↑ — *both i and j flipped* |

### Cell-count vs. node-count extents

Patch IJK extents are in **cells**, not nodes. A face perpendicular
to i on a block of `Ni × Nj × Nk` *nodes* (so `Ni-1 × Nj-1 × Nk-1`
cells) has:

```
Imin = 0,  Imax = 0       (i_lo face)        Imin = Ni-1,  Imax = Ni-1   (i_hi face)
Jmin = 0,  Jmax = Nj-1                      ↳ same J extents
Kmin = 0,  Kmax = Nk-1                      ↳ same K extents
```

SHORE's adjacency JSON records node-shape `(Ni, Nj, Nk)`; the writer
subtracts 1 to produce the cell-count extents the reader expects.

## SHORE → `cc.par` BC mapping

How `shore cc-par` translates each SHORE `FaceBC` enum value into the
Xall BC code. See the [BC code constants](#bc-code-constants) table
above for the full Xall BC vocabulary.

| SHORE `FaceBC` | `cc.par` BC code | `family` field | Notes |
|---|---|---|---|
| `WALL` | `1` (`parete`, viscous wall) | `0` | Body wall. |
| `FREE` | `--free-bc` (default `40`), or per-block `free_face_bc` via `--meta`, or per-**face** `bc_code` | `0` | `40` = chimera cell-center recipient (body-fitted / buffer meshes). For far-field background blocks, override with a **positive** natural BC in `[1, 19]`: `3` inflow, `4` outflow, `9` extrapolation (see warning below). Per-block `--meta` overrides `--free-bc`; a per-face `bc_code` (see [per-face BC](#per-face-bc-codes)) overrides both. |
| `SHARED` | orientation code (≥ 100) | partner patch index | Code derived from the face's `axis_map` plus the partner face name — see [orientation digits](#orientation-digits-in-detail). For a same-axis seam with no flips the codes are `135` / `145` per side. |
| `DIRICHLET` | orientation code | partner patch index | Same derivation as `SHARED`; cross-axis seams (e.g. cubed-sphere cap-equator) are produced naturally from the `axis_map` without per-topology lookup tables. |
| `PERIODIC` | `40` + `UserWarning` | `0` | `cc.par` has no first-class periodic BC; use Xall chimera boxes for periodic stitching. |

::: warning Natural-BC codes for FREE faces must be POSITIVE
The overset preprocessor stores a natural-BC code **verbatim** into its `icc`
matrix and then rejects any face with `icc < 1`
(`Errore in calspig: icc sulle facce deve essere > 0`). Use a **positive**
code in `[1, 19]` — `9` for free/far-field extrapolation, `3`/`4` for
inflow/outflow. A **negative** code such as `-9` is invalid and fails where two
natural faces meet at an edge. (Earlier docs recommended `-9`; that was wrong.)
:::

### Per-face BC codes {#per-face-bc-codes}

A single block sometimes needs *heterogeneous* outer BCs — e.g. inflow on one
face and extrapolation on the others. A `FREE` [`Face`](/reference/api-mesh)
can carry an explicit `bc_code` (positive, `[1, 19]`) that overrides the
per-block `free_face_bc` for that face only:

```python
from shore.mesh.face import FaceBC
# inflow (3) on the -X face; this block's other FREE faces stay default (9)
blk._set_face_bc("i_lo", FaceBC.FREE, bc_code=3)
```

The code is carried through the adjacency sidecar and emitted by `cc-par`. See
[Design decisions → per-face BC](/guide/design-decisions#per-face-natural-bc-codes).

## Cubed-sphere cap-equator codes

For the cubed-sphere topology, each sub-block's polar face
(`sub_q.i_lo` for south, `sub_q.i_hi` for north) connects to a
*lateral* face of the matching cap block (`cap_*.i_lo` / `i_hi` /
`j_lo` / `j_hi`), with a per-`q` axis remapping that mirrors the cap
boundary overwrite in `topology/sphere_cap.py`. The eight orientation
codes that fall out of the axis_map derivation are:

| Sub | Pole | Cap face | Sub-side code | Cap-side code |
|---|---|---|---|---|
| sub0 | north | `cap_north.i_lo` | `135` | `235` |
| sub1 | north | `cap_north.j_hi` | `415` | `325` |
| sub2 | north | `cap_north.i_hi` | `245` | `245` |
| sub3 | north | `cap_north.j_lo` | `325` | `425` |
| sub0 | south | `cap_south.j_lo` | `315` | `315` |
| sub1 | south | `cap_south.i_hi` | `235` | `135` |
| sub2 | south | `cap_south.j_hi` | `425` | `415` |
| sub3 | south | `cap_south.i_lo` | `145` | `145` |

These values are not hard-coded in the writer: they are computed
from the `axis_map` field that `CubedSphereTopology` stamps onto each
seam at build time, then serialised into the adjacency JSON. The
table is shown here as a reference for cross-checking the writer's
output, not as a writer input.

## CLI options

| Flag | Default | Description |
|---|---|---|
| `--grd TEXT` | *required* | Name of the `.grd` mesh file Xall pairs with this `cc.par` (header line 1). |
| `--output-basename TEXT` | `"cc."` | Output basename Xall uses for derived files (header line 2). |
| `--free-bc INT` | `40` | BC code for `FREE` faces on every block that has no entry in `--meta`. Use a **positive** natural BC (`3`/`4`/`9`) for far-field background blocks — see the warning above. |
| `--meta TEXT` | `null` | Per-block `BlockMeta` overrides. Accepts an inline JSON object or a path to a JSON file (inline JSON is tried first). Keys are block labels; values are objects with any subset of `{level, group, priority, comment, free_face_bc}`. Missing keys inherit their defaults. Any block absent from the map uses `--free-bc` for its FREE faces. |
| `--mgl INT INT` | `4 1` | Multigrid levels: total, finest active. |
| `--boundary-layer-thickness FLOAT` | `-0.1` | Header line 9. |
| `--numerical-beach FLOAT` | `-1.0` | Header line 10. |

### `--meta` format

`--meta` accepts either an inline JSON string or a path to a JSON file.
Resolution order: the value is first parsed as JSON; if that fails it is
treated as a file path (the file's content is then parsed as JSON). This
means inline JSON is always unambiguous regardless of the working
directory.

The JSON value must be an **object** mapping block labels to partial
`BlockMeta` fields:

```json
{
  "sub0":      { "level": 1, "group": 0, "priority": 1, "free_face_bc": 40 },
  "cap_north": { "level": 2, "free_face_bc": 9 }
}
```

Any field can be omitted; omitted fields fall back to their defaults
(`level=1`, `group=0`, `priority=1`, `free_face_bc=--free-bc`). Blocks
absent from the map entirely get all defaults plus `free_face_bc` from
`--free-bc`.

Error cases:
- Value is not valid JSON **and** is not a path to an existing file → exit 1.
- Value is a valid file path but the file content is not valid JSON → exit 1.
- Value parses as JSON but is not an object (e.g. a JSON array) → exit 1.

## Examples

```bash
# Body-fitted near-body mesh (FREE faces are chimera fringes)
shore cc-par sphere.adjacency.json sphere.cc.par --grd sphere.grd

# Far-field background mesh — all FREE faces are extrapolation
shore cc-par background.adjacency.json background.cc.par \
  --grd background.grd --free-bc -9

# Mixed near-body: caps get a different overset level and far-field BC
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
  --meta '{"cap_north": {"level": 2, "free_face_bc": -9},
           "cap_south": {"level": 2, "free_face_bc": -9}}'

# Same, loading the metadata from a JSON file
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
  --meta wall_meta.json
```

## Python API

```python
from shore.io.cc_par import BlockMeta, write_cc_par

write_cc_par(
    "wall.adjacency.json",     # path or parsed dict
    "wall.cc.par",
    grd_filename="wall.grd",
    output_basename="cc.",
    multigrid_levels=(4, 1),
    boundary_layer_thickness=-0.1,
    numerical_beach=-1.0,
    block_metadata={
        # Equatorial sub-blocks — chimera cell-centre handoff
        "sub0": BlockMeta(free_face_bc=40),
        "sub1": BlockMeta(free_face_bc=40),
        "sub2": BlockMeta(free_face_bc=40),
        "sub3": BlockMeta(free_face_bc=40),
        # Polar caps — different overset level, far-field BC
        "cap_north": BlockMeta(level=2, free_face_bc=-9),
        "cap_south": BlockMeta(level=2, free_face_bc=-9),
    },
)
```

The `block_metadata` dict is optional; any block not listed gets a
default `BlockMeta()` (level=1, group=0, priority=1, free_face_bc=40).
The `--meta` CLI flag exposes the same control without writing Python.

## Multi-component runs

Independent Chimera components (e.g. a near-body wall + a far-field
background) are produced as separate `.grd` + `cc.par` pairs and
combined for a single Xall job by:

- [`shore cc-par-merge`](./cc-par-merge) — concatenate the
  `cc.par` files with the right block-index and patch-`family`
  shifts.
- [`shore grd-merge`](./grd-merge) — concatenate the binary `.grd`
  files in the same order.

The two commands must be invoked with the same input ordering; see
each page for details.

## Marker boxes

The boxes section of `cc.par` (cells flagged as internal-wall by
the chimera hole-cutter) has its own page covering the data model,
the `--boxes` flag, the `shore boxes-from-stl` generator, and the
`shore boxes-export` ParaView visualiser:

- [Marker boxes](./boxes) — full reference.
- [Algorithm — marker boxes](/algorithm/boxes) — voxel-fill +
  greedy-merge derivation.

## Limitations

- **No `edges` section.** Unsupported by SHORE-style adjacency.
- **One-way bridge.** No `cc.par` reader is provided; the Xall
  Fortran tooling owns the canonical reader.
- **Topology-agnostic orientation.** Cross-axis seams are no longer
  hard-coded per topology: every connection face in the adjacency
  JSON carries an `axis_map` that the writer uses to derive the
  orientation code. New topologies need only stamp valid `axis_map`
  values on their seams (see `shore.mesh.face.AxisMap`); no edits
  to this writer are required.

## See also

- [`shore grd`](./grd) — pack per-block `.geo` files into the binary
  `.grd` that pairs with this `cc.par`.
- [`shore cc-par-merge`](./cc-par-merge) — multi-component merge.
- [`shore grd-merge`](./grd-merge) — companion `.grd` merger.
- [`shore split`](./split) — produces the `.adjacency.json` SHORE
  consumes.
- [Adjacency sidecar format](./split#adjacency-sidecar) — JSON schema.
- [Marker boxes](./boxes) — `cc.par` boxes section reference.
- [overset-exploded](https://github.com/szaghi/overset-exploded) —
  the Fortran reference implementation of `cc.par` I/O.

