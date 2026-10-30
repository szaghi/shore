# Project scripts

Everything SHORE does to a grid is a **command** applied to a `Project`
(`shore.project`, [ADR 0004](/adr/0004-command-layer)). The CLI grid
commands, the GUI and your own scripts all go through it, and the
project records each command as one line of a **journal**: a Python
script that rebuilds the same grids bit for bit. This page is the
reference of that layer: the `shore run` command, every `Project`
method, and the journal and project files. For the ideas behind it see
[Projects and journals](/concepts/projects); for a walk-through,
[tutorial 6](/tutorials/projects).

## `shore run`

Replay a journal: a Python script of `Project` method calls, written by
`--journal` (on `mesh`, `primitive *`, `cap-tip`, `split`, `grd`,
`cc-par`, `proc-input`) or `Project.save_journal`.

```
shore run [OPTIONS] JOURNAL
```

| Argument / option | Default | Description |
|---|---|---|
| `JOURNAL` | required | The journal script. Paths in it are used as written, relative to the directory you run from. |
| `--save-project DIR` | — | Also save the resulting project directory (`project.toml`, `journal.py`, `inputs/`, `adjacency.json`, `grids/<component>/*.npz`; see [File formats](formats#project-directory)). |

```bash
shore body sphere -o sphere.stl
shore mesh sphere.stl --ni 33 --nj 64 --nk 25 --ds 0.01 -o out/sphere --journal sphere.py
shore run sphere.py                          # rewrites out/sphere_*.geo identically
shore run sphere.py --save-project sphere.shore
```

The journal of that run:

```python
# SHORE journal: replay with `shore run <this file>`.
from shore.project import Project

p = Project()
p.load_geometry('sphere.stl')
p.build_topology('cubed_sphere', nk=25, ds=0.01, ni=33, nj=64, theta_cap=0.15707963267948966, i_spacing='uniform', i_beta=3.0, j_spacing='uniform', j_beta=3.0)
p.march()
p.write_geo('out/sphere')
```

Edit it like any script: change a parameter, add a `p.snap(exact=True)`,
point `load_geometry` at an edited CAD file, and replay.

::: warning A journal is Python
`shore run`, `Project.load` and `run_journal` execute the file. Run only
journals you trust.
:::

## The `Project` API

```python
from shore.project import CancelToken, Project, run_journal
```

Each method below builds one command (a frozen dataclass in
`shore.project.commands`) and applies it; the command's `to_source()` is
its journal line. Methods act on the **current component**
(`use_component`). State commands can be undone; output commands
(`write_*`) are recorded but never undone and are not saved with the
project. A command that raises leaves the project as it was.

### Geometry and components

| Method | Command | What it does | ADR |
|---|---|---|---|
| `load_geometry(path, anchor_strategy="auto", skip_hygiene=False, *, name="body")` | `LoadGeometry` | Load an STL, or STEP/IGES/BREP (`cad` extra), under `name`. Reloading a name replaces it and discards the meshes of the components built on it. | [0002](/adr/0002-geometry-backend), [0016](/adr/0016-multi-component-projects) |
| `use_component(name)` | `UseComponent` | Make component `name` current, creating it empty if new. The first component is `"main"`. | [0016](/adr/0016-multi-component-projects) |

### Topologies and marching

| Method | Command | What it does | ADR |
|---|---|---|---|
| `build_topology(topology, *, nk, ds=1e-3, growth=1.15, k_law="geometric", k_beta=3.0, k_thickness=None, options=None, geometry="body", **params)` | `BuildTopology` | The k=0 wall and the blocks of `cubed_sphere`, `ogrid`, `cgrid`, `ch` or `oh` on `geometry`. `params` are the topology's (`ni`, `nj`, `theta_cap`, … or `ni_body`, `ni_wake`, `n_stations`, `wake_length`, …), checked by name; `options` the `ch` / `oh` frame. See [`shore mesh`](mesh). | [0004](/adr/0004-command-layer) |
| `march(smooth_strength=0.2, smooth_iters=2, blend_normals_k=None, smooth_mode="limited", march_scheme="algebraic", *, progress=None, cancel=None)` | `March` | Extrude the blocks from the wall. `progress(fraction, message)` and a `CancelToken` are hooks; a failed or cancelled march restores the grid. `blend_normals_k` is removed: `0` warns, anything else is refused. | [0005](/adr/0005-marching-schemes), [0034](/adr/0034-remove-normal-blending) |
| `build_primitive(kind, **params)` | `BuildPrimitive` | `box`, `annulus`, `flat_caps`, `h_grid_rect` or `h_grid_tube` with the primitive's `build` parameters (`vmin`, `vmax`, `ni`, …). See [`shore primitive`](primitive). | |
| `build_tip_cap(body_geo, **options)` | `BuildTipCap` | Wing-tip cap blocks from a C-grid body `.geo`; options as [`shore cap-tip`](cap-tip) (`tip`, `le_frac`, `te_frac`, …). Replaces the mesh. | |
| `load_blocks(geo_paths, adjacency=())` | `LoadBlocks` | Replace the mesh with blocks read from `.geo` files, wired by `.adjacency.json` sidecars. | |
| `split(splits)` | `Split` | Split blocks by `(label, axis, at)` directives, as [`shore split`](split). | |

### Association

| Method | Command | What it does | ADR |
|---|---|---|---|
| `associate(block, entity, geometry="body", kind="body", indices=())` | `Associate` | Link face, edge or vertex `entity` of a block with entities of `geometry` (`body`, `surface`, `curve`, `point`), or clear it with `kind="none"`. | [0003](/adr/0003-association), [0039](/adr/0039-association-editing) |
| `snap(exact=False)` | `Snap` | Project the associated nodes onto their geometry (vertex wins over edge over face). `exact=True` projects onto the exact CAD surfaces. Seams stay closed. | [0003](/adr/0003-association), [0040](/adr/0040-seams-kept-closed) |

### Spacing and wall lines

| Method | Command | What it does | ADR |
|---|---|---|---|
| `redistribute(block, edge, law="uniform", *, beta=3.0, cluster_at=None, first_step=None, last_step=None, interior="uniform")` | `Redistribute` | Apply a spacing law along an edge to every line of its axis, through the seams. | [0041](/adr/0041-redistribute-edge) |
| `orthogonalize(strength=1.0)` | `Orthogonalize` | Reshape the lines across the walls so they leave the wall along its normal. | [0048](/adr/0048-orthogonal-wall-lines) |
| `wall_layer(first_step)` | `WallLayer` | Cluster the nodes on those lines geometrically from the wall, the first cell exactly `first_step`. | [0050](/adr/0050-wall-layer) |

### Blocking

| Method | Command | What it does | ADR |
|---|---|---|---|
| `create_blocking(ni, nj, nk, *, geometry="body", margin=0.0, vmin=None, vmax=None)` | `CreateBlocking` | Start the component's blocking with one block `b0`. | [0042](/adr/0042-blocking-skeleton-premesh) |
| `split_blocking(block, axis, fraction)` | `SplitBlocking` | Cut a block, and its sheet, at a node. | [0042](/adr/0042-blocking-skeleton-premesh) |
| `delete_blocks(blocks)` | `DeleteBlocks` | Remove blocking blocks (a hole). | [0042](/adr/0042-blocking-skeleton-premesh) |
| `ogrid(blocks=(), *, faces=(), exclude=(), factor=0.5, nodes=5)` | `Ogrid` | Insert an O-grid layer on faces of a block set. | [0043](/adr/0043-ogrid-insertion) |
| `move_vertex(block, vertex, xyz)` | `MoveVertex` | Move a blocking vertex, onto its association if it has one. | [0045](/adr/0045-move-blocking-vertex) |
| `associate_blocking(block, entity, geometry="body", kind="body", indices=())` | `AssociateBlocking` | `associate` for blocking faces, edges and vertices. | [0042](/adr/0042-blocking-skeleton-premesh), [0049](/adr/0049-stl-feature-curves) |
| `premesh(exact=False)` | `Premesh` | Make the component's mesh from its blocking. | [0042](/adr/0042-blocking-skeleton-premesh) |

### Collars and surface grids

| Method | Command | What it does | ADR |
|---|---|---|---|
| `build_collar(other, *, seed, n_loop, ds_surface, growth_surface, m_a, m_b, ds_wall, growth_wall, nk, label="collar", plane=None, geometry="body", loop_law="uniform", loop_beta=3.0, loop_cluster_at=None, loop_curvature=0.0, split=False, wake_length=None, volume=None)` | `BuildCollar` | The collar block at the junction of `geometry` and `other` (a geometry name or a file path); the half-collar with `plane=(point, normal)`; `split=True` cuts it at sharp corners. Replaces the current component's mesh. Options as [`shore collar`](collar). | [0010](/adr/0010-collar-block), [0013](/adr/0013-half-collar-block), [0020](/adr/0020-junction-curve-spacing), [0024](/adr/0024-split-collar-at-sharp-corners), [0026](/adr/0026-smoothed-collar-volume-by-default) |
| `surface_grid(block, face, *, surface, smoothing="elliptic", geometry="body")` | `FillFace` | Rebuild a block face on CAD SURFACE `surface` from its boundary: TFI in (u, v), then elliptic smoothing (or `"none"`). | [0007](/adr/0007-surface-grids) |

### Outputs

Output commands write files and change no state.

| Method | Command | What it writes | ADR |
|---|---|---|---|
| `write_geo(stem, per_block=False)` | `WriteGeo` | The `.geo` files: a primitive's own naming, else one `<stem>.geo` for a single non-cubed-sphere block and `<stem>_<label>.geo` per block otherwise (always the latter with `per_block`). | |
| `write_adjacency(path)` | `WriteAdjacency` | The `.adjacency.json` sidecar. | [0003](/adr/0003-association) |
| `write_grd(path, geo_paths)` | `WriteGrd` | A `.grd` packed from `.geo` files, as `shore grd`. | |
| `write_cc_par(adjacency, path, *, grd_filename, output_basename="cc.", free_face_bc=40, multigrid_levels=(4, 1), boundary_layer_thickness=-0.1, numerical_beach=-1.0, meta={}, boxes=())` | `WriteCcPar` | A `cc.par` from a sidecar, as `shore cc-par`. | [0036](/adr/0036-cc-par-as-xall-reads-it) |
| `write_proc_input(grd, path, *, np=1, meta={}, header_comment="generated by SHORE")` | `WriteProcInput` | A `proc.input`, as `shore proc-input`. | |
| `write_xall(directory, *, grd="mesh.grd", np=1, components=None, output_basename="cc.", multigrid_levels=(4, 1), boundary_layer_thickness=-0.1, numerical_beach=-1.0)` | `WriteXall` | The Xall run of every component: per component files, merged `.grd` / `cc.par`, one `proc.input`; `components` maps names to `free_face_bc`, `meta`, `boxes`, `body`, `proc_meta`, as [`shore xall`](xall). | [0016](/adr/0016-multi-component-projects), [0017](/adr/0017-xall-options-per-component) |
| `write_cgns(path, components=())` | `WriteCgns` | A CGNS file (`cgns` extra), all components with a mesh by default. | [0035](/adr/0035-cgns-plot3d-export) |
| `write_plot3d(path, components=())` | `WritePlot3d` | A multi-block Plot3D grid file. | [0035](/adr/0035-cgns-plot3d-export) |

### History, journal and files

| Member | What it does | ADR |
|---|---|---|
| `undo()` | Revert the last state-changing command (and drop the outputs recorded after it); returns it. Snapshots live in an on-disk `UndoCache`. Refused for history loaded without replay. | [0004](/adr/0004-command-layer) |
| `redo()` | Re-apply the last undone command. A new command clears the redo stack. | [0004](/adr/0004-command-layer) |
| `history` | The applied commands, oldest first. | |
| `journal()` / `save_journal(path)` | The history as a replayable script / written to a file. | [0004](/adr/0004-command-layer) |
| `save(directory)` | The project directory: inputs copied into `inputs/` with their SHA-256, the state-changing journal, the grids of every component. | [0016](/adr/0016-multi-component-projects), [0033](/adr/0033-load-without-replay) |
| `Project.load(directory, *, replay=True)` | Check the inputs, replay the journal and require bit-identical grids; with `replay=False`, read the saved grids instead (no march, no undo, no blocking). | [0016](/adr/0016-multi-component-projects), [0033](/adr/0033-load-without-replay) |
| `run_journal(path, input_dir=None)` | Execute a journal and return its `p`. | |
| `close()` | Delete the on-disk undo snapshots. | |
| `mesh`, `blocking`, `blockings`, `geometries`, `components`, `current_component` | The current component's mesh and blocking, every component's blocking, the loaded geometries (by name), every component's mesh, the current component's name. | |

## A whole Xall component in one script

```python
from pathlib import Path
from shore.project import Project

p = Project()
p.load_geometry("sphere.stl")
p.build_topology("cubed_sphere", nk=13, ni=30, nj=40, ds=0.01, growth=1.1)
p.march()
p.split([(lbl, "k", (6,)) for lbl in ("sub0", "sub1", "sub2", "sub3", "cap_north", "cap_south")])
p.write_geo("out/wall", per_block=True)
p.write_adjacency("out/wall.adjacency.json")
p.write_grd("out/wall.grd", sorted(Path("out").glob("wall_*.geo")))
p.write_cc_par("out/wall.adjacency.json", "out/wall.cc.par", grd_filename="wall.grd")
p.write_proc_input("out/wall.grd", "out/wall.proc.input", np=4)
```

`write_grd` packs from the `.geo` files, as `shore grd` does, so the
`.grd` carries the `.geo` print precision. The `*-merge` commands,
`balance`, `boxes-*`, the `body`/`profile` generators and the inspection
commands are file tools outside the project. Progress and cancel:

```python
from shore.project import CancelToken

token = CancelToken()              # token.cancel() from another thread stops the march
p.march(progress=lambda f, msg: print(f"{f:5.0%} {msg}"), cancel=token)
p.undo()                           # back to the un-marched grid, bit for bit
p.redo()
p.save("sphere.shore")
```

## Associations

A block face (`"k_lo"`), edge (`"i_lo_k_lo"`) or vertex (`"v0"` .. `"v7"`)
is linked to geometry with `associate`; `snap` then projects its nodes
onto it, vertices before edges before faces winning at shared nodes
([ADR 0003](/adr/0003-association),
[ADR 0039](/adr/0039-association-editing)):

```python
p.load_geometry("part.step", name="part")
p.build_primitive("box", vmin=(-0.8, -0.4, 0.55), vmax=(0.8, 0.4, 0.9), ni=5, nj=4, nk=3)
p.associate("box", "k_lo", "part", "surface", (2,))    # a face with CAD surfaces
p.associate("box", "j_lo_k_lo", "part", "curve", (5,))  # an edge with a curve
p.associate("box", "v0", "part", "point", (1,))         # a vertex with a point
p.snap()
p.associate("box", "k_lo", kind="none")                 # clear
```

An entity takes geometry of its own dimension or higher (a face: surfaces
or the body; an edge: also curves; a vertex: also points). Associations are
stored as entity fingerprints, so they are found again after the CAD file
is re-imported, and they are saved with the project. `snap` resolves every
geometry of the project by name; `"body"` is the mesh's own body.

## Redistributing an edge

`redistribute` applies a spacing law along a block edge's axis
([ADR 0041](/adr/0041-redistribute-edge)):

```python
p.redistribute("sub0", "j_lo_k_lo", "tanh2", beta=2.0)          # cluster at both ends
p.redistribute("box", "j_lo_k_lo", "pinned", first_step=0.01, last_step=0.2)
```

The law (`uniform`, `tanh`, `tanh2` with `beta` and `cluster_at`, `pinned`
with `first_step` / `last_step` lengths and an `interior` law) is taken over
the edge's arc length. Every grid line of the block along that axis is
resampled at the same arc-length fractions of its own length, and so is
every block that shares a face containing the axis through a SHARED or
DIRICHLET seam (reversed across a flipped seam). Nodes stay on their old
lines, line ends do not move, and seams stay closed. A face shared with
several blocks (CH's body) or a periodic face on the way is refused.
Associated nodes are not re-snapped: run `snap` after.

## Blocking

A blocking is a block topology made by hand, ICEM-style: a skeleton of hex
blocks with a node count on every edge, edited by commands and premeshed
into ordinary blocks ([ADR 0042](/adr/0042-blocking-skeleton-premesh)).
Twenty-six blocks around a body:

```python
p = Project()
p.load_geometry("sphere.stl")
p.create_blocking(13, 13, 13, vmin=(-3, -3, -3), vmax=(3, 3, 3))  # or margin=1.0
p.split_blocking("b0", "i", 1 / 3)       # cut at node round(fraction * (n - 1))
p.split_blocking("b1", "i", 0.5)         # three columns along i: b0, b1, b2
p.split_blocking("b0", "j", 1 / 3)       # cuts b0, b1, b2: high parts b3, b4, b5
p.split_blocking("b3", "j", 0.5)
p.split_blocking("b0", "k", 1 / 3)
p.split_blocking("b9", "k", 0.5)         # 27 blocks; b13 is the centre one
for face in ("i_lo", "i_hi", "j_lo", "j_hi", "k_lo", "k_hi"):
    p.associate_blocking("b13", face)    # the centre block's faces on the body
p.delete_blocks(["b13"])                 # the hole
p.premesh()
```

- `create_blocking(ni, nj, nk, *, geometry="body", margin=0.0, vmin=None,
  vmax=None)`: one block `b0`, the geometry's bounding box grown by
  `margin` times its extent on each side, or the box between the corners.
- `split_blocking(block, axis, fraction)`: cuts at the nearest node, and
  every block sharing a cut edge with it. The low part keeps its label; the
  high part of each cut block is `b<n>`, counting up in block order.
- `delete_blocks(labels)`: the faces left alone become boundary faces.
- `associate_blocking(block, entity, geometry, kind, indices)`: as
  `associate`, on the blocking. An entity without an association of its
  own takes the one its faces agree on (an edge), or its edges or faces
  (a vertex). An edge between faces associated with different geometry
  stays straight until it gets one of its own.
- `premesh(exact=False)`: replaces the component's mesh. Edges are equally
  spaced and projected; faces and blocks are filled by transfinite
  interpolation, faces projected. Shared faces become SHARED seams,
  associated boundary faces WALL, the others FREE, and the blocks keep the
  associations, so `snap` and `redistribute` work on them.

Replaying the journal with an edited geometry re-projects the blocking. A
project loaded without replay has its mesh but no blocking.

### O-grids

`ogrid` inserts a layer of blocks on faces of a block set
([ADR 0043](/adr/0043-ogrid-insertion)):

```python
p.ogrid(["b13"], factor=0.5, nodes=5)                 # a core and a ring of 6 blocks
p.ogrid(["b0"], exclude=[("b0", "k_lo")])             # half-O: open on k_lo
p.ogrid(faces=[("b4", "k_hi"), ("b22", "k_lo"), ...])  # a layer on the faces of a hole
```

- The block set is `blocks`, or every block when empty. The layered faces
  are `faces` (`(block, face)` pairs), or every face on the set's boundary
  less `exclude`.
- Every vertex of a layered face gets an inner copy, moved into the set by
  `(1 - factor)` times half the sum of its edges that leave the layered
  faces into the set. Inside one block, the core is the block scaled by
  `factor` about its centre.
- The set's blocks take the copies and keep their labels. Each layered face
  gets a block between it and its copy, with `nodes` nodes across; these are
  labelled `b<n>` in block and face order.
- The set's other boundary faces stay open: the O-grid passes through them,
  and copies on them stay in their plane. An open face may not be shared
  with a block outside the set (add that block to the set, or layer the
  face).
- Associations:
  - a layered face keeps its own, on the layer's outer face;
  - a replaced face or edge passes its association on;
  - an open face also passes it to the layer faces lying on it.

Around a body, layer the faces of the hole (or O-grid the centre block
before deleting its core): the grid lines then leave the wall nearly normal
to it (within 6 degrees around a sphere, against 41 for the plain lattice).

### Curves and points of an STL

An STL has CURVE and POINT entities, like a CAD model, for associations
([ADR 0049](/adr/0049-stl-feature-curves)).
- **Curves:** chains of its feature edges, the edges where the facets turn
  by more than `crease_deg` (30°), and its open boundaries. A wing's TE
  corner lines are curves, and so are its root and tip outlines.
- **Points:** the vertices where chains meet or end.

`p.geometries["body"].geometry.features()` lists them. Then
`p.associate("box", "i_lo_j_lo", "body", "curve", [k])` or
`p.associate_blocking(...)` puts an edge on curve `k`, and the GUI picks
curves as it does CAD edges. A smooth body (a sphere) has none.

### Moving vertices

`p.move_vertex("b0", "v7", (1.2, 1.3, 1.1))` moves vertex `v7` of block
`b0`, and every block sharing it
([ADR 0045](/adr/0045-move-blocking-vertex)). A vertex with an
association (its own, or inherited from its edges or faces) goes to the
projection of the point onto it, so it slides on its geometry. After a
geometry edit, premesh re-projects it.

### Wall-normal grid lines

After `premesh`, the lines across a wall run along the transfinite
interpolation and leave a curved wall obliquely. `p.orthogonalize(strength=1.0)`
reshapes them ([ADR 0048](/adr/0048-orthogonal-wall-lines)).

- In every block with one WALL face, each line across the wall becomes a
  curve that leaves the wall along its normal (tangent `strength` times its
  chord) and meets the opposite face as before.
- The wall and the opposite face stay put, and each line keeps its nodes'
  arc-length fractions.
- Blocks without a wall whose faces move are re-filled.
- A strength that would invert a valid block is refused.

Around a wing on a plane, the median angle between the first grid segment
and the wall normal goes from 30 to 8 degrees. Then `redistribute`
(clustering towards the wall) moves the nodes along the new lines.

### Wall layers

`p.wall_layer(1e-4)` clusters the nodes towards the walls
([ADR 0050](/adr/0050-wall-layer)).
- On every line across the wall of a block with one WALL face, the cells
  grow geometrically from the wall, the first exactly `first_step` long.
  The ratio comes from each line's length.
- The nodes stay on the line, and the wall node and the far node stay put.

Run it after `orthogonalize`, so that the clustered lines leave the wall
along its normal: `premesh → orthogonalize → wall_layer(ds)`. It needs a
layer topology (an O-grid around the body). On a plain lattice it would
invert the corner blocks and is refused, as is a first step that a line
cannot hold.
