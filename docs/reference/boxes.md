# Marker boxes (`shore boxes-from-stl`, `shore boxes-export`)

Marker boxes are the last section of `cc.par` (see
[cc.par format reference](./cc-par#boxes-section)). Each box defines
a closed hexahedral region of the simulation domain; cells inside
the box are flagged as **internal wall** by the chimera hole-cutter
(`overset/src/scatole.f90`'s `check_parete_interna`).

This page covers the box data model, the two CLI commands that
generate (`boxes-from-stl`) and visualise (`boxes-export`) marker
boxes, the `--boxes` flag on `shore cc-par`, and the underlying
algorithm. The Python API surface lives in
[`shore.io.boxes`](./api-boxes-io) (the `Box` dataclass) and
[`shore.boxes`](./api-boxes) (`boxes_from_stl`); see
[Algorithm — marker boxes](/algorithm/boxes) for the
voxel-fill + greedy-merge derivation.

## Use case

A near-body mesh (e.g. the cubed sphere wall) contains a solid body
explicitly; the background mesh that overlaps the same volume must
mark its body-interior cells as donors-only. The canonical way is
to define one or more boxes per body whose union:

- contains the STL surface of the body, and
- itself fits inside the body's near-mesh (i.e. between the STL and
  the outer extent of the near-mesh).

A simple body like a sphere can be covered by a small handful of
slab-shaped boxes. Slender or branched bodies usually need more.

## Vertex ordering (type 9)

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

## Block association

Each box references one SHORE block via `Box.associated_block` (the
block's label). At write time the label is resolved to a 1-based
index into the `cc.par` blocks list; the box inherits the block's
`level / group / priority` for the chimera precedence rules.

Choose the block whose mesh you want the box to override: for a
near-body wall, that's typically one of the wall blocks themselves
(level=2 say); the box then hole-cuts a coarser background block
(level=1).

## Generating boxes from an STL: `shore boxes-from-stl`

Hand-writing 8-vertex hexahedra for every body is tedious and
error-prone. `shore boxes-from-stl` (and the underlying
[`shore.boxes.boxes_from_stl`](./api-boxes)) take a body STL and
return a list of `Box` instances ready to feed into
`write_cc_par(boxes=...)` or `shore cc-par --boxes`.

```bash
shore boxes-from-stl sphere.stl --block sub0 --max-outward-gap 0.30 \
                     -o sphere_boxes.json
```

### Methods

| Method | Description | When to use |
|---|---|---|
| `voxel-fill` (default) | Voxelise the body's AABB at side `max_outward_gap / sqrt(3)`. Include every voxel that overlaps the body (`signed_distance(centre) > -half_diag`) — both centre-inside and body-shell voxels. Greedy-merge face-adjacent grid blocks into slab-shaped AABBs whenever the merged block's 8 corners still satisfy `outward_gap <= max_outward_gap`. The output covers the entire body (interior + surface) with a small slab cover; box count tracks body surface area, not volume, after merging. | Always — the only method that produces a watertight cover of the body. |
| `obb` | Single oriented bounding box from area-weighted PCA. Tight. Best effort: emits a `UserWarning` when the measured outward gap exceeds `max_outward_gap`. | Quick previews / sanity checks. |
| `aabb` | Single world-axis-aligned bounding box. Tight. Same best-effort warning. | Quick previews on body-fitted-cube bodies. |

### Two constraints, one algorithm

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

`voxel-fill` satisfies both by construction. See
[Algorithm — marker boxes](/algorithm/boxes) for the full
derivation; the short version:

- **Coverage (1).** The seed grid admits every voxel whose centre is
  inside the body (`signed > 0`) plus every body-shell voxel whose
  centre is at most one half-diagonal outside (`-half_diag < signed
  <= 0`). Together these tile the entire body — interior and surface
  shell — so any background cell-centre inside the body lands in at
  least one box.
- **Confinement (2).** Voxel side is `max_outward_gap / sqrt(3)`,
  giving `2 * half_diag = side * sqrt(3) = max_outward_gap`: the
  worst-case outward gap from a body-shell voxel's corner to the
  body surface is exactly the budget. The greedy merge pass
  rechecks the 8 corners of every candidate merged AABB before
  accepting it, so merging never violates the budget.

### Picking `max_outward_gap`

Supply a value smaller than the body's near-mesh wall-normal
thickness. For the `chimera_assembly` example (unit-radius sphere,
wall thickness ~0.214), `max_outward_gap=0.30` is conservative and
produces a few dozen merged boxes on the icosphere.

The seed grid scales as body-volume / voxel-volume =
`body-volume * (sqrt(3) / gap)^3`, so halving the gap roughly 8×
the *seed* count. After greedy merging, the final box count tracks
body surface area rather than volume, so halving the gap roughly
4× the *merged* count. Xall reads the boxes once at startup and
the per-cell box-test loop is `O(N_boxes)` — a few thousand boxes
per body is fine.

### Watertightness

`voxel-fill` calls `trimesh.contains()` for the inside/outside
test, which assumes a watertight body STL. Non-watertight inputs
are rejected with a clear `MeshInputError` recommending repair.
The single-box methods (`aabb`, `obb`) work on non-watertight STLs
too — they don't need inside/outside testing.

### Why single-box methods are best-effort

For most bodies, a single OBB has corner outward-gap dominated by
the body's bounding-cube geometry: `r * (sqrt(3) - 1)` past a
sphere of radius `r`. On a unit-radius sphere that's ~0.73 — much
larger than typical near-mesh thickness budgets. The single-box
methods emit one box and warn when the measured outward gap exceeds
the budget, recommending the user switch to `voxel-fill`. They
remain useful as quick previews and for very loose budgets on
near-cubic bodies.

### Python API

```python
from shore.boxes import boxes_from_stl

boxes = boxes_from_stl(
    "sphere.stl",
    associated_block="sub0",
    max_outward_gap=0.30,
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

## Composing with `shore cc-par`: `--boxes`

`shore cc-par` accepts a `--boxes` flag that mirrors `--meta`:
either inline JSON or a path to a JSON file. The payload is a list
of `{"associated_block": <label>, "vertices": [[x,y,z], ...]}`
entries:

```bash
# Generate, then bake into cc.par
shore boxes-from-stl sphere.stl --block sub0 --max-outward-gap 0.30 \
                     -o sphere_boxes.json
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
             --boxes sphere_boxes.json

# Or hand-write the JSON inline
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
             --boxes '[
                {"associated_block": "sub0",
                 "vertices": [[-0.6,-0.6,-0.6], [0.6,-0.6,-0.6],
                              [-0.6,0.6,-0.6], [0.6,0.6,-0.6],
                              [-0.6,-0.6,0.6], [0.6,-0.6,0.6],
                              [-0.6,0.6,0.6], [0.6,0.6,0.6]]}
             ]'
```

A box referencing a block label not present in the adjacency JSON
exits with a clear "unknown block" error. Geometric validation
errors (degenerate hexahedron, coplanar vertices, etc.) surface as
non-zero exits with a message identifying the failing check.

The Python equivalent:

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

## Visualising boxes in ParaView: `shore boxes-export`

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

## Merge behaviour

[`shore cc-par-merge`](./cc-par-merge) (and `merge_cc_par` in
Python) concatenate boxes across inputs. Each box's
associated-block index is shifted by the cumulative block count of
the preceding inputs — the same offset that applies to patch
`block` fields — so a box originally pointing at "block 1" of the
second input lands at "block 1 + len(A.blocks)" in the merged file.
Box vertex coordinates pass through unchanged.

Boxes of types other than 9 (parallelepiped / cylinder shorthands)
are tolerated by the merger and round-tripped verbatim, even though
SHORE itself never emits them.

## Limitations

- **Type 9 only on emit.** SHORE emits hexahedral (type-9) boxes
  only. Other types in input files are passed through by the merger
  but never produced by the writer.
- **`voxel-fill` requires a watertight STL.** Repair the STL before
  calling, or fall back to `aabb` / `obb` for a best-effort cover.

## See also

- [`shore cc-par`](./cc-par) — the `cc.par` writer that consumes
  boxes via `--boxes`.
- [`shore cc-par-merge`](./cc-par-merge) — multi-component merge.
- [Algorithm — marker boxes](/algorithm/boxes) — voxel-fill +
  greedy-merge derivation.
- [`shore.io.boxes`](./api-boxes-io) — the `Box` dataclass.
- [`shore.boxes`](./api-boxes) — `boxes_from_stl`.
- [cc.par + boxes pipeline](/guide/cc-par-pipeline) — end-to-end
  walkthrough.
