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
