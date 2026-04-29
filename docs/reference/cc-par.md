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

```fortran
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
| `FREE` | `--free-bc` (default `40`), or per-block `free_face_bc` via `--meta` | `0` | `40` = chimera cell-center recipient (body-fitted / buffer meshes). For far-field background blocks, override with a natural BC: `-3` inflow, `-4` inoutflow, `-9` extrapolation. Per-block overrides via `--meta` take precedence over `--free-bc`. |
| `SHARED` | orientation code (≥ 100) | partner patch index | Code derived from the face's `axis_map` plus the partner face name — see [orientation digits](#orientation-digits-in-detail). For a same-axis seam with no flips the codes are `135` / `145` per side. |
| `DIRICHLET` | orientation code | partner patch index | Same derivation as `SHARED`; cross-axis seams (e.g. cubed-sphere cap-equator) are produced naturally from the `axis_map` without per-topology lookup tables. |
| `PERIODIC` | `40` + `UserWarning` | `0` | `cc.par` has no first-class periodic BC; use Xall chimera boxes for periodic stitching. |

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
| `--free-bc INT` | `40` | BC code for `FREE` faces on every block that has no entry in `--meta`. Use `-3`/`-4`/`-9` etc. for far-field background blocks. |
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
  "cap_north": { "level": 2, "free_face_bc": -9 }
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

## Merging multiple `cc.par` files (`shore cc-par-merge`)

Several independent Chimera components in one Xall run (e.g. a
near-body wall + a far-field background) are generated as separate
`.grd` + `cc.par` pairs. To submit them as a single Xall job, the two
`cc.par` files have to be combined into one whose blocks list is the
concatenation of the inputs, with patch block-indices and `family`
pointers shifted so they still refer to the right targets.

`shore cc-par-merge` does the merge:

```bash
shore cc-par-merge background.cc.par wall.cc.par \
                   -o assembly.cc.par --grd assembly.grd
```

### What gets shifted, what doesn't

For inputs in caller-specified order `[A, B, C, ...]`:

- **Block list**: simple concatenation. Block 1 of `B` becomes
  block `len(A.blocks) + 1` in the merged file, and so on.
- **Patch `block` field**: each patch in `B` has its 1-based block
  index shifted by `len(A.blocks)`; same recipe per subsequent input.
- **Patch `family` field**: `cc.par` connection patches reference
  their partner by 1-based patch index. Since each input is
  internally complete (no SHARED across the merge boundary), every
  family pointer stays *inside* its own input — but the absolute
  index changes. The merger shifts `family` by `len(A.patches)` for
  patches in `B`, etc.
- **Patch extents** (Imin/Imax/Jmin/Jmax/Kmin/Kmax): unchanged.

Symmetric pairing survives the merge: every connection patch with
`bc ≥ 100` continues to point at a partner patch whose `family`
points back at it.

### Header reconciliation

Header fields other than the `.grd` filename are taken from `inputs[0]`
verbatim (output basename, ghost-cells flag, multigrid levels,
boundary-layer thickness, numerical-beach width). The `.grd` filename
on header line 1 is overridden by the `--grd` flag — that's how the
merged file knows which `.grd` it pairs with.

If subsequent inputs disagree on any header field, a `UserWarning` is
emitted but the merge proceeds with `inputs[0]`'s value. Surface drift
loudly without blocking.

### Block-comment auto-prefixing

Each block line in the merged `cc.par` has its trailing comment
prefixed with the input file's stem so blocks from different
topologies are distinguishable when reading the merged file:

```
        1         0         1 ! 1 : background.center
        2         0         1 ! 6 : wall.sub0
        2         0         1 ! 10 : wall.cap_north
```

`cc.par` treats comments as advisory only — Xall indexes blocks by
their 1-based position — so this is purely for human readability.
A trailing `.cc` in the input filename (typical for `wall.cc.par`)
is stripped from the prefix.

### Python API

```python
from shore.io.cc_par import merge_cc_par

merge_cc_par(
    ["background.cc.par", "wall.cc.par"],
    "assembly.cc.par",
    grd_filename="assembly.grd",
)
```

Raises `ValueError` for fewer than 2 inputs, `FileNotFoundError` for
a missing input.

### Coordinated invocation

`cc-par-merge` and `grd-merge` must be called with the **same** input
ordering — `cc.par` patch block-indices are 1-based positions into
the `.grd` block list. A typical end-to-end pair:

```bash
shore grd-merge    background.grd    wall.grd    -o assembly.grd
shore cc-par-merge background.cc.par wall.cc.par -o assembly.cc.par \
                   --grd assembly.grd
```

The example at `examples/cc_par_assembly.py` builds the wall and
background meshes from scratch, runs both merges, and prints the
equivalent CLI commands.

## Marker boxes

Marker boxes are the last section of `cc.par`. Each box defines a
closed hexahedral region of the simulation domain; cells inside the
box are flagged as **internal wall** by the chimera hole-cutter
(`overset/src/scatole.f90`'s `check_parete_interna`).

**Use case.** A near-body mesh (e.g. the cubed sphere wall) contains
a solid body explicitly; the background mesh that overlaps the same
volume must mark its body-interior cells as donors-only. The
canonical way is to define one or more boxes per body whose union:

- contains the STL surface of the body, and
- itself fits inside the body's near-mesh (i.e. between the STL and
  the outer extent of the near-mesh).

A simple body like a sphere can use a single AABB centred on the
body. Slender or branched bodies usually need several boxes.

### Vertex ordering (type 9)

SHORE emits boxes in type-9 (hexahedron) form: 8 vertices per box in
the canonical "i-fastest" Fortran-column-major order, where vertex
index `1 + i + 2*j + 4*k` corresponds to corner `(i, j, k) ∈
{0, 1}³`:

```
1 = (i=lo, j=lo, k=lo)
2 = (i=hi, j=lo, k=lo)
3 = (i=lo, j=hi, k=lo)
4 = (i=hi, j=hi, k=lo)
5 = (i=lo, j=lo, k=hi)
6 = (i=hi, j=lo, k=hi)
7 = (i=lo, j=hi, k=hi)
8 = (i=hi, j=hi, k=hi)
```

The `(i, j, k)` directions need not be world-aligned: any rotated or
sheared parallelepiped is valid as long as the connectivity above is
preserved (faces `{1,2,3,4}` and `{5,6,7,8}`, opposite face pairs
`{1,3,5,7}` ↔ `{2,4,6,8}` and `{1,2,5,6}` ↔ `{3,4,7,8}`).

`shore.io.boxes.Box.__post_init__` rejects degenerate, coplanar,
self-intersecting (bow-tie), and inverted hexahedra at construction.

### Block association

Each box references one SHORE block via `Box.associated_block` (the
block's label). At write time the label is resolved to a 1-based
index into the `cc.par` blocks list; the box inherits the block's
`level / group / priority` for the chimera precedence rules.

Choose the block whose mesh you want the box to override: for a
near-body wall, that's typically one of the wall blocks themselves
(level=2 say); the box then hole-cuts a coarser background block
(level=1).

### Python API

```python
from shore.io.boxes import Box
from shore.io.cc_par import write_cc_par

# AABB convenience: corners + label.
box1 = Box.from_aabb(
    vmin=(-0.6, -0.6, -0.6),
    vmax=(+0.6, +0.6, +0.6),
    associated_block="sub0",
)

# Or explicit 8 vertices for a skewed / rotated box.
import numpy as np
verts = np.array([...])  # (8, 3) in canonical order
box2 = Box(vertices=verts, associated_block="cap_north")

write_cc_par(
    "wall.adjacency.json",
    "wall.cc.par",
    grd_filename="wall.grd",
    boxes=[box1, box2],
)
```

### CLI: `--boxes`

`shore cc-par` accepts a `--boxes` flag that mirrors `--meta`:
either inline JSON or a path to a JSON file. The payload is a list
of `{"associated_block": <label>, "vertices": [[x,y,z], ...]}`
entries:

```bash
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
             --boxes '[
                {"associated_block": "sub0",
                 "vertices": [[-0.6,-0.6,-0.6], [0.6,-0.6,-0.6],
                              [-0.6,0.6,-0.6], [0.6,0.6,-0.6],
                              [-0.6,-0.6,0.6], [0.6,-0.6,0.6],
                              [-0.6,0.6,0.6], [0.6,0.6,0.6]]}
             ]'

shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
             --boxes boxes.json
```

A box referencing a block label not present in the adjacency JSON
exits with a clear "unknown block" error. Geometric validation
errors (degenerate hexahedron, coplanar vertices, etc.) surface as
non-zero exits with a message identifying the failing check.

### Merge behaviour

`shore cc-par-merge` (and `merge_cc_par` in Python) concatenate boxes
across inputs. Each box's associated-block index is shifted by the
cumulative block count of the preceding inputs — the same offset
that applies to patch `block` fields — so a box originally pointing
at "block 1" of the second input lands at "block 1 + len(A.blocks)"
in the merged file. Box vertex coordinates pass through unchanged.

Boxes of types other than 9 (parallelepiped / cylinder shorthands)
are tolerated by the merger and round-tripped verbatim, even though
SHORE itself never emits them.

### Generating boxes from an STL

Hand-writing 8-vertex hexahedra for every body is tedious and
error-prone. `shore.boxes.boxes_from_stl` (and the
`shore boxes-from-stl` CLI) take a body STL and return a list of
`Box` instances ready to feed into `write_cc_par(boxes=...)` or
`shore cc-par --boxes`.

Three algorithms, picked via the `method` argument:

| Method | Description | When to use |
|---|---|---|
| `voxel-fill` (default) | Voxelise the body's AABB at side `2 * max_outward_gap / sqrt(3)`. Test each voxel-centre with `trimesh.contains`; emit each strictly-inside-body voxel as an axis-aligned `Box`. Box vertex outward gap to the body surface is bounded by `max_outward_gap` by construction. Guarantees body-volume coverage. | Always — the only method that covers the body interior. |
| `obb` | Single oriented bounding box from area-weighted PCA. Tight. Best effort: emits a `UserWarning` when the measured outward gap exceeds `max_outward_gap`. | Quick previews / sanity checks. |
| `aabb` | Single world-axis-aligned bounding box. Tight. Same best-effort warning. | Quick previews on body-fitted-cube bodies. |

#### Two constraints, one algorithm

The boxes' role in chimera hole-cutting imposes two simultaneous
geometric constraints:

1. **Coverage of body volume.** Every body-interior cell-centre of
   the background mesh must lie inside at least one box, so the
   chimera marks it as donors-only.
2. **Confinement to the near-mesh.** Every box vertex must lie
   inside the body's near-mesh (the wall mesh). Box vertices that
   sit past the near-mesh outer extent spuriously hole-cut
   legitimate flow cells.

A single AABB or OBB satisfies (1) trivially (it bounds the body)
but typically fails (2) because the bounding cube's corners stick
out by `r * (sqrt(3) - 1)` past a circumscribed sphere of radius
`r` — much further than typical near-mesh thicknesses.

`voxel-fill` satisfies both by construction:

- The voxel side is chosen so the half-diagonal does not exceed
  `max_outward_gap`. A voxel with centre strictly inside the body
  therefore has every vertex at most `max_outward_gap` past the
  body surface in the worst case (corner outside, on a body-tangent
  voxel). Constraint (2) holds.
- The voxel set tiles the body's deep interior. Coverage of (1) is
  exact in the discretisation: every cell-centre at `signed_distance
  >= side / 2` from the surface is guaranteed to lie inside some
  voxel, and shallower body-interior cells are typically already
  inside the wall mesh and don't need box marking.

#### Picking `max_outward_gap`

Supply a value smaller than the body's near-mesh wall-normal
thickness. For the `cc_par_assembly` example (cubed sphere, wall
thickness ~0.214), `max_outward_gap=0.10` is conservative and
produces ~389 boxes on the icosphere.

Box count scales with body-volume / voxel-volume = body-volume *
`(sqrt(3) / (2 * gap))^3`. Halving the gap budget produces ~8x more
boxes (3D voxel scaling). Xall reads the boxes once at startup and
the per-cell box-test loop is O(N_boxes) — a few thousand boxes per
body is fine.

#### Watertightness

`voxel-fill` calls `trimesh.contains()` for the inside/outside
test, which assumes a watertight body STL. Non-watertight inputs
are rejected with a clear `MeshInputError` recommending repair.
The single-box methods (`aabb`, `obb`) work on non-watertight STLs
too — they don't need inside/outside testing.

#### Why single-box methods are best-effort

For most bodies, a single OBB has corner outward-gap dominated by
the body's bounding-cube geometry: `r * (sqrt(3) - 1)` past a
sphere of radius `r`. On a unit-radius sphere that's ~0.73 — much
larger than typical near-mesh thickness budgets. The single-box
methods emit one box and warn when the measured outward gap exceeds
the budget, recommending the user switch to `voxel-fill`. They
remain useful as quick previews and for very loose budgets on
near-cubic bodies.

#### Python API

```python
from shore.boxes import boxes_from_stl

boxes = boxes_from_stl(
    "sphere.stl",
    associated_block="sub0",
    max_outward_gap=0.10,
    method="voxel-fill",         # default
)
write_cc_par(
    "wall.adjacency.json",
    "wall.cc.par",
    grd_filename="wall.grd",
    boxes=boxes,
)
```

For best-effort `aabb` / `obb`, an additional `face_samples_per_axis`
kwarg (default 4) controls the box-surface sampling density used
for the outward-gap post-check. Each face is sampled on a `k × k`
grid (8 + 6 * k * k samples per box). Raise to 6-8 for sharply
concave geometries.

#### CLI

```bash
shore boxes-from-stl sphere.stl --block sub0 --max-outward-gap 0.10 \
                     -o sphere_boxes.json

shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
             --boxes sphere_boxes.json
```

The output is a JSON list in the exact schema `shore cc-par --boxes`
consumes — the two commands compose without intermediate massaging.

### Visualising boxes in ParaView

Marker boxes are usually wrong the first time you place them — too
small misses body cells, too large hole-cuts legitimate near-mesh
cells. Eyeballing the boxes against the STL surface and the wall
mesh is the fastest way to validate them.

Two paths produce ParaView-ready files:

**Standalone** — boxes only, useful for opening in a separate
ParaView pane next to the existing `assembly.vtm`:

```bash
shore boxes-export sphere_boxes.json -o sphere_boxes.vtm
```

```python
from shore.io.vtk import write_boxes_vtm

write_boxes_vtm(boxes, "sphere_boxes.vtm")
```

**Layered** — boxes appended to the existing per-topology `.vtm` so
the body wall, near-mesh, and marker boxes show in one MultiBlock:

```python
from shore.io.vtk import write_caps_vtm

write_caps_vtm(
    "wall.vtm",
    grids=[wall_block_grids...],
    labels=[wall_block_labels...],
    boxes=body_boxes,
)
```

Both forms emit one `vtkUnstructuredGrid` cell per box, of cell type
`VTK_HEXAHEDRON`. The vertex order is remapped from SHORE's
canonical `(i, j, k)` bit layout to VTK's hexahedron winding so
ParaView renders the box as a translucent solid by default — switch
to wireframe via the GUI for a clearer view of the box vs. the
mesh.

Each cell carries two cell-data fields:

- `associated_block` (string): the SHORE block label the box was
  attached to. Useful for colouring boxes by topology in mixed
  scenes.
- `box_index` (int, 1-based): the index in the original `Box` list,
  matching the `MultiBlock` child key (`box_<n>`). Useful for
  picking a specific box or filtering by selection.

## Limitations

- **Boxes: type 9 only on emit.** SHORE emits hexahedral (type-9)
  boxes only. Other types in input files are passed through by the
  merger but never produced by the writer.
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
- [`shore grd-merge`](./grd#merging-multiple-grd-files-shore-grd-merge) —
  the companion merge for `.grd` binaries. Always invoke
  `cc-par-merge` and `grd-merge` as a pair, with the same input
  ordering.
- [`shore split`](./split) — produces the `.adjacency.json` SHORE
  consumes.
- [Adjacency sidecar format](./split#adjacency-sidecar) — JSON schema.
- [overset-exploded](https://github.com/szaghi/overset-exploded) —
  the Fortran reference implementation of `cc.par` I/O.
