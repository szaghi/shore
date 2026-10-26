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

Every change goes through a command (`LoadGeometry`, `BuildTopology`,
`March`, `Snap`, `WriteGeo`, `WriteAdjacency`). `undo` reverts the last
state-changing command and drops any outputs written after it. A failed or
cancelled `march` leaves the grid as it was before the call.
