# cc.par + boxes pipeline

This walkthrough shows the full chimera-component pipeline:
generate two independent meshes (a body-fitted cubed-sphere wall +
a flat-caps cylinder background), produce one `.grd` and one
`cc.par` per component, generate the marker boxes that hole-cut
the body, then **merge** the two pairs into a single Xall-ready
assembly.

The reference example is at `examples/05-chimera-assembly/cc_par_assembly.py` — every
step in this guide maps to a function call in that script, which
runs end-to-end without external data.

## What you'll build

Two meshes covering the same domain `[-10, 10]² × [-10, 50]`:

| Mesh | Topology | Blocks | Role |
|---|---|---|---|
| Near-body wall | [Cubed sphere](/topologies/cubed-sphere) | 6 | Wraps a unit-radius (diameter 1.0) sphere body at the origin, marched outward to wall-normal extent ~0.2. |
| Background | [Flat-caps cylinder](/topologies/primitives) | 5 | Cylindrical domain with no central hole — central square + 4 trapezoidal sub-blocks. |

Plus, packaged for the solver:

- 2 binary `.grd` files (one per mesh) → 1 merged `assembly.grd`.
- 2 ASCII `cc.par` files (one per mesh, the wall one with marker
  boxes) → 1 merged `assembly.cc.par`.
- 1 ParaView `.vtm` per mesh (and optionally a standalone boxes `.vtm`) for visual inspection.

## Step 1 — Generate the meshes

The two meshes are exactly what
[`assembly_sphere_in_cylinder.py`](/guide/assembly) builds — see
that walkthrough for the topology, parameter choices, and ParaView
inspection. The short form:

```python
from shore.mesh import CubedSphereTopology, Mesh
from shore.primitive import FlatCapsCylinder

# Near-body wall: 6 cubed-sphere blocks
wall = Mesh.from_stl(
    "sphere.stl",
    topology=CubedSphereTopology(),
    ni=30, nj=40, nk=13,
    theta_cap_deg=30.0,
)
wall.march(ds=0.01, growth=1.1)
wall.write_geo("wall")          # wall_sub0.geo .. wall_cap_south.geo

# Background: 5 flat-caps blocks
bg = FlatCapsCylinder.build(
    r_out=10.0, z0=-10.0, z1=50.0,
    ni=16, nc=24, nk=60,
    inner_frac=0.4,
)
bg.write_geo("background")      # background_center.geo .. background_sub_s.geo
```

## Step 2 — Per-component .grd

Pack each mesh's `.geo` files into a binary `.grd`. The block
ordering you choose here is the block ordering the companion
`cc.par` will reference, so be consistent.

```python
from shore.io.grd import write_grd

write_grd("wall.grd", [
    "wall_sub0.geo", "wall_sub1.geo",
    "wall_sub2.geo", "wall_sub3.geo",
    "wall_cap_north.geo", "wall_cap_south.geo",
])

write_grd("background.grd", [
    "background_center.geo",
    "background_sub_e.geo", "background_sub_n.geo",
    "background_sub_w.geo", "background_sub_s.geo",
])
```

CLI equivalent:

```bash
shore grd wall_sub0.geo wall_sub1.geo wall_sub2.geo wall_sub3.geo \
          wall_cap_north.geo wall_cap_south.geo -o wall.grd
shore grd background_center.geo background_sub_e.geo background_sub_n.geo \
          background_sub_w.geo background_sub_s.geo -o background.grd
```

## Step 3 — Per-component adjacency sidecars

`shore cc-par` needs one `.adjacency.json` per mesh — the JSON
sidecar that captures every face's BC, partner, and `axis_map`.
Both `Mesh` (cubed sphere) and `FlatCapsCylinder` already wire
their internal seams; emit the sidecar with `write_adjacency_json`:

```python
from shore.io.adjacency import write_adjacency_json

write_adjacency_json(wall.blocks, "wall.adjacency.json")
write_adjacency_json(bg.blocks,   "background.adjacency.json")
```

If you need to subdivide blocks (parallel decomposition), call
[`shore split`](/reference/split) before this step — it produces
the same JSON sidecar on the split block list.

## Step 4 — Generate marker boxes from the body STL

Boxes hole-cut the background mesh's body-interior cells. Generate
them from the STL with [`shore boxes-from-stl`](/reference/boxes):

```python
from shore.boxes import boxes_from_stl

boxes = boxes_from_stl(
    "sphere.stl",
    associated_block="sub0",        # any wall block; level/group inherited
    max_outward_gap=0.10,           # < wall-normal extent (~0.2)
    method="voxel-fill",            # default; covers body interior + surface
)
```

`max_outward_gap` must be smaller than the wall mesh's wall-normal
thickness (so the boxes fit inside the wall). The voxel-fill
algorithm + greedy merge produces a small number of slab-shaped
boxes that cover the entire body. See
[Algorithm — marker boxes](/algorithm/boxes) for the derivation.

CLI equivalent:

```bash
shore boxes-from-stl sphere.stl --block sub0 --max-outward-gap 0.10 \
                     -o sphere_boxes.json
```

## Step 5 — Per-component cc.par

For the **wall** mesh (body-fitted), `FREE` faces are chimera
fringes → BC code `40` (`offgen`, the default). The cubed-sphere
caps can use a different overset level if you want them to take
priority over the equatorial sub-blocks. Bake the boxes into the
wall's `cc.par`:

```python
from shore.io.cc_par import BlockMeta, write_cc_par

write_cc_par(
    "wall.adjacency.json",
    "wall.cc.par",
    grd_filename="wall.grd",
    block_metadata={
        "sub0": BlockMeta(level=1, group=0, priority=1, free_face_bc=40),
        "sub1": BlockMeta(level=1, group=0, priority=1, free_face_bc=40),
        "sub2": BlockMeta(level=1, group=0, priority=1, free_face_bc=40),
        "sub3": BlockMeta(level=1, group=0, priority=1, free_face_bc=40),
        "cap_north": BlockMeta(level=2, group=0, priority=1, free_face_bc=40),
        "cap_south": BlockMeta(level=2, group=0, priority=1, free_face_bc=40),
    },
    boxes=boxes,
)
```

For the **background** mesh (far-field), `FREE` faces are
extrapolation → BC code `-9`. No boxes needed (boxes live with the
body):

```python
write_cc_par(
    "background.adjacency.json",
    "background.cc.par",
    grd_filename="background.grd",
    block_metadata={
        "center": BlockMeta(level=0, group=1, priority=1, free_face_bc=-9),
        "sub_e":  BlockMeta(level=0, group=1, priority=1, free_face_bc=-9),
        "sub_n":  BlockMeta(level=0, group=1, priority=1, free_face_bc=-9),
        "sub_w":  BlockMeta(level=0, group=1, priority=1, free_face_bc=-9),
        "sub_s":  BlockMeta(level=0, group=1, priority=1, free_face_bc=-9),
    },
)
```

CLI equivalent:

```bash
shore cc-par wall.adjacency.json wall.cc.par --grd wall.grd \
             --boxes sphere_boxes.json \
             --meta '{"cap_north": {"level": 2}, "cap_south": {"level": 2}}'

shore cc-par background.adjacency.json background.cc.par \
             --grd background.grd --free-bc -9
```

## Step 6 — Merge into the assembly

Two coordinated calls — `merge_grd` and `merge_cc_par` **must use
the same input ordering**. The patch `block` field, `family`
pointers, and box `associated_block` indices are 1-based positions
into the merged blocks list, so the two files only agree if they
were merged with the same ordering.

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

CLI equivalent:

```bash
shore grd-merge    background.grd    wall.grd    -o assembly.grd
shore cc-par-merge background.cc.par wall.cc.par -o assembly.cc.par \
                   --grd assembly.grd
```

The merged files describe one 11-block Xall job: 5 background
blocks (indices 1..5) + 6 wall blocks (indices 6..11).

## Step 7 — Visualise

ParaView is the cheapest sanity check — opening the merged `.vtm`
shows that the wall sits where you think it does and that the
boxes actually cover the body:

```python
from shore.io.vtk import write_caps_vtm

# Wall + boxes in one MultiBlock
write_caps_vtm(
    "wall.vtm",
    [b.nodes for b in wall.blocks],
    labels=[b.label for b in wall.blocks],
    boxes=boxes,
)

# Background alone
write_caps_vtm(
    "background.vtm",
    [b.nodes for b in bg.blocks],
    labels=[b.label for b in bg.blocks],
)
```

Open both in ParaView. Switch the boxes to wireframe (the default
is translucent solid) for a clearer view of how the boxes
intersect the wall and the body STL.

## Reference

- [`examples/05-chimera-assembly/cc_par_assembly.py`](https://github.com/szaghi/shore/blob/main/examples/05-chimera-assembly/cc_par_assembly.py) — the full script.
- [`shore cc-par`](/reference/cc-par) — file format reference.
- [`shore cc-par-merge`](/reference/cc-par-merge) — multi-component merge.
- [`shore grd`](/reference/grd) and [`shore grd-merge`](/reference/grd-merge) — binary mesh writers.
- [Marker boxes](/reference/boxes) — full reference.
- [Algorithm — marker boxes](/algorithm/boxes) — voxel-fill + greedy-merge derivation.
- [Assembly walkthrough](/guide/assembly) — the underlying mesh build, without `cc.par`.
