---
title: Projects and replay
---

# Projects and replay

Every grid SHORE builds through a `Project` is the result of a list of
commands, and that list is a Python script: the journal. In this tutorial
you build a project with two components, undo and redo, save it, load it
with and without replay, write the Xall run and a CGNS file, and then
replay the journal on an edited input to get a new grid with no other
change. The last step does the same from the command line.

::: info You need
- SHORE installed, plus the `cgns` extra for the CGNS step
  (`pip install 'shore-mesh[cgns]'`).
- About 2 s for the script
  [`examples/tutorials/projects.py`](https://github.com/szaghi/shore/blob/main/examples/tutorials/projects.py).
:::

## 1. Build two components

A project holds named geometries and named components
([ADR 0016](/adr/0016-multi-component-projects)). Each mesh command acts on
the current component, `"main"` by default; `use_component` switches to
another one, creating it empty if it is new. A geometry is loaded once and
any component can build on it.

<<< @/../examples/tutorials/projects.py#components{python}

```text
components: {'main': 1, 'background': 1}
```

`main` holds an O-grid marched from the sphere; `background` a box of
13 × 13 × 13 nodes round it. In a Chimera run they are two independent
component grids.

## 2. Undo and redo

Every state-changing command records how to revert itself
([ADR 0004](/adr/0004-command-layer)). `undo` reverts the last one (output
commands such as `write_xall` are skipped: they are never undone), `redo`
applies it again.

<<< @/../examples/tutorials/projects.py#undo{python}

```text
after undo: None
after redo: ['box']
```

## 3. Read the journal

`journal()` renders the history as a script of the same method calls.
`save_journal` writes it; `shore run journal.py` (or `run_journal`) replays
it and gives the same grids bit for bit.

<<< @/../examples/tutorials/projects.py#journal{python}

```python
# SHORE journal: replay with `shore run <this file>`.
from shore.project import Project

p = Project()
p.load_geometry('sphere.stl', name='sphere')
p.build_topology('ogrid', nk=17, ds=0.002, growth=1.2, geometry='sphere', ni=17, nj=32)
p.march()
p.use_component('background')
p.build_primitive('box', vmin=(-3.0, -3.0, -3.0), vmax=(3.0, 3.0, 3.0), ni=13, nj=13, nk=13)
```

The undone and redone box appears once: the journal is the history as it
stands, not a log of what you typed. A journal is Python, so run only
journals you trust.

## 4. Save and load

`save(dir)` writes a portable project directory: `project.toml` (format,
components and their block labels, input hashes), `journal.py`, `inputs/`
(a copy of every file the journal reads) and `grids/<component>/` (the node
arrays and each component's face graph).

<<< @/../examples/tutorials/projects.py#save{python}

```text
replay == saved grids: True
```

The two loads differ:

| | `Project.load(dir)` | `Project.load(dir, replay=False)` |
|---|---|---|
| Inputs | checked against their SHA-256 | checked against their SHA-256 |
| Grids | rebuilt by replaying the journal, then required to equal the saved ones bit for bit | read from `grids/`, rewired from the saved face graphs |
| Cost | the whole build again (marches, premesh) | reading the files |
| Afterwards | undo and redo work as in the original session | the loaded history cannot be undone, and a loaded mesh cannot be marched |

A replaying load is also a determinism check: it fails if anything gives a
different grid. Load without replay ([ADR 0033](/adr/0033-load-without-replay))
when you only want to inspect, export or add components.

## 5. Write the solver files

`write_xall` writes each component's `.geo`, adjacency sidecar, `.grd` and
`cc.par` under `run/<component>/`, then the merged `.grd` and `cc.par` and
one `proc.input` balanced over all blocks of the run. `write_cgns` writes one
CGNS base per component.

<<< @/../examples/tutorials/projects.py#outputs{python}

```text
['background', 'cc.par', 'main', 'mesh.grd', 'proc.input']
```

Per-component options (marker boxes, `free_face_bc`, body numbers, per-block
metadata) are in [Multi-component Xall run](/cookbook/multi-component).

## 6. Replay on an edited input

The journal names its inputs by path. `run_journal(path, input_dir=d)`
resolves relative paths against `d`, so the same commands run on a
different file of the same name. Here the sphere grows from radius 0.5 to
0.6:

<<< @/../examples/tutorials/projects.py#replay{python}

```text
edited wall radius: 0.598
```

The wall now lies on the larger sphere (0.598: the mean over the STL
facets), with the same node counts, spacing and topology. This is how a
grid follows a design change. For a CAD model, an edit that keeps the
model's faces and edges recognisable replays the same way: the
associations re-find their entities by fingerprint
([ADR 0003](/adr/0003-association)). [Blocking a wing–body](/tutorials/blocking)
does it on a STEP file whose span and fuselage radius change.

## 7. The same from the command line

`--project DIR [--component C]` on `shore mesh`, `shore collar` and
`shore primitive *` loads the project in `DIR` (or creates it), builds
component `C` and saves the project back
([ADR 0018](/adr/0018-cli-collars-and-projects)). `shore xall` and
`shore cgns` write from a saved project. The commands below use
`sphere.stl` and `journal.py` from the steps above and `box.stl` /
`post.stl` from the [collar tutorial](/tutorials/collar).

```bash
shore mesh sphere.stl --topology ogrid --ni 17 --nj 32 --nk 17 \
    --ds 2e-3 --growth 1.2 -o out/sphere --project proj --journal sphere_journal.py
shore collar box.stl post.stl --seed 0.55 0 0.02 --n-loop 48 \
    --ds-surface 0.01 --growth-surface 1.12 --m-a 12 --m-b 12 \
    --ds-wall 1e-4 --growth-wall 1.2 --nk 25 --project proj
shore xall proj --out run --options xall.toml
shore cgns proj -o proj.cgns
shore cgns proj -o collar.cgns -c collar --no-replay
shore run journal.py --save-project replayed
```

with `xall.toml`:

```toml
[run]
grd = "run.grd"
np = 4

[components.collar]
body = 1
```

- `shore mesh` builds component `main` (its default), `shore collar`
  component `collar`. A geometry already loaded from a file with the same
  content is reused, so the earlier components survive; a file whose stem
  matches a geometry loaded from a different file is an error (name it
  with `--name-a`, `--name-b` or `--geometry-name`).
- `shore xall` writes `run/main/`, `run/collar/`, then `run/run.grd`,
  `run/cc.par` and `run/proc.input`.
- `--no-replay` on any of them opens the project from its saved grids.

## What you learned

- A project is geometries, components and the command history that built
  them; `use_component` selects where mesh commands act.
- `undo`/`redo` walk the history; the journal is that history as a script.
- `save` makes a portable directory; `load` replays and verifies it,
  `load(replay=False)` reads the saved grids.
- Replaying the journal on edited inputs rebuilds the grid with the same
  commands.
- `--project` and `--component` give the CLI the same model;
  `shore xall` and `shore cgns` write from it.

## Next

- [Blocking a wing–body](/tutorials/blocking): a 46-command journal
  replayed after a CAD edit.
- [Projects and journals](/concepts/projects): the command layer.
- [Save, load, replay](/cookbook/projects) and
  [Export CGNS and Plot3D](/cookbook/exports).
- [Project scripts](/reference/run), [`shore xall`](/reference/xall),
  [`shore cgns`](/reference/cgns).
