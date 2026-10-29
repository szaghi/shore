# shore run

Replay a SHORE journal: a Python script of `Project` method calls, written by
`shore mesh --journal` or `Project.save_journal`. Replaying re-runs the same
commands and rewrites the same outputs, bit for bit (ADR 0004).

```bash
shore mesh sphere.stl --ni 40 --nj 64 --nk 30 -o out/sphere --journal sphere.py
shore run sphere.py                          # rewrites out/sphere_*.geo identically
shore run sphere.py --save-project sphere.shore
```

A journal reads like this:

```python
# SHORE journal: replay with `shore run <this file>`.
from shore.project import Project

p = Project()
p.load_geometry('sphere.stl')
p.build_topology('cubed_sphere', nk=30, ni=40, nj=64, theta_cap=0.15707963267948966)
p.march()
p.write_geo('out/sphere')
```

Edit it like any script: change a parameter, add a `p.snap(exact=True)`,
point `load_geometry` at an edited CAD file, and replay. Paths are used as
written, relative to the directory you run from.

::: warning A journal is Python
`shore run` executes the file. Run only journals you trust.
:::

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `JOURNAL` | — | The journal script. |
| `--save-project DIR` | — | Also save the resulting project directory: `journal.py`, `project.toml`, `adjacency.json`, `grids/*.npz`. `Project.load(DIR)` replays the journal and requires every grid to match the saved one bit for bit. |

## Python API

```python
from shore.project import Project, CancelToken

p = Project()
p.load_geometry("wing.step")
p.build_topology("cgrid", nk=40, ni_body=120, ni_wake=30, n_stations=12, wake_length=10.0)
token = CancelToken()
p.march(progress=lambda f, msg: print(f"{f:5.0%} {msg}"), cancel=token)
p.undo()          # back to the un-marched grid, bit for bit
p.redo()
p.save("wing.shore")
```

Every change goes through a command:

| Group | Commands |
|---|---|
| Build a mesh | `LoadGeometry`, `BuildTopology`, `March`, `Snap`, `BuildPrimitive`, `BuildTipCap`, `LoadBlocks` |
| Associate with geometry | `Associate` (`p.associate(block, entity, geometry, kind, indices)`, ADR 0039), then `Snap` |
| Change the spacing | `Redistribute` (`p.redistribute(block, edge, law, ...)`, ADR 0041) |
| Transform it | `Split` |
| Write outputs | `WriteGeo`, `WriteAdjacency`, `WriteGrd`, `WriteCcPar`, `WriteProcInput` |

`shore mesh`, `shore primitive *`, `shore cap-tip`, `shore split`,
`shore grd`, `shore cc-par` and `shore proc-input` all run over a `Project`
and accept `--journal`. One script can carry a whole Xall component:

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
commands are file tools outside the project. `undo` reverts the last
state-changing command and drops any outputs written after it. A failed or
cancelled `march` leaves the grid as it was before the call.

## Associations

A block face (`"k_lo"`), edge (`"i_lo_k_lo"`) or vertex (`"v0"` .. `"v7"`)
is linked to geometry with `associate`; `snap` then projects its nodes
onto it, vertices before edges before faces winning at shared nodes
([ADR 0003](../adr/0003-association.md),
[ADR 0039](../adr/0039-association-editing.md)):

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
([ADR 0041](../adr/0041-redistribute-edge.md)):

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
into ordinary blocks ([ADR 0042](../adr/0042-blocking-skeleton-premesh.md)).
Twenty-six blocks around a body:

```python
p = Project()
p.load_geometry("sphere.stl")
p.create_blocking(13, 13, 13, vmin=(-3, -3, -3), vmax=(3, 3, 3))  # or margin=1.0
p.split_blocking("b0", "i", 1 / 3)       # cut at node round(fraction * (n - 1))
p.split_blocking("b1", "i", 0.5)
...                                      # the same along j and k
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
