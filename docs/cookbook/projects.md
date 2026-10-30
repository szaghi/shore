# Save, load and replay a project

Use these recipes when a grid takes more than one command: a `Project` records
every command in a journal, so you can undo it, save the whole state to a
directory, load it back (replayed and checked, or fast without replay), and
edit the journal as a script.

The script runs in about 2 s.

## Undo and redo a command

<<< @/../examples/cookbook/projects.py#undo_redo{python}

```text
undo restores: True
['LoadGeometry', 'BuildTopology', 'March', 'Redistribute']
```

- Every `Project` method that changes state is one command: `undo()` restores
  the state before it bit for bit, `redo()` re-applies it. A new command after
  an undo clears the redo stack.
- Grid snapshots for undo live in an on-disk cache, not in memory.
- Output commands (`write_geo`, `write_adjacency`, `write_xall`, ...) are
  recorded in the journal but never undone: the files stay written.

## Save a project and load it back

<<< @/../examples/cookbook/projects.py#save_load{python}

```text
['adjacency.json', 'grids', 'inputs', 'journal.py', 'project.toml']
replayed load 0.17 s, load without replay 0.02 s
same blocks: True same journal: True
```

The project directory holds:

| Entry | Content |
|---|---|
| `project.toml` | format, SHORE version, command count, each component's block labels, input hashes |
| `journal.py` | the history as a replayable script, output commands dropped |
| `inputs/` | copies of the geometry files (SHA-256 checked on load) |
| `grids/<component>/` | each component's blocks (`.npz`) and face graph (`adjacency.json`) |
| `adjacency.json` | the face graph of the current component |

- `Project.load(dir)` checks the input hashes, replays the journal and
  requires every rebuilt grid to equal the saved one bit for bit; a changed
  grid raises `MeshInputError`.
- The directory is portable: inputs are copied in and the journal refers to
  them relative to the project.

## Load a project fast, without replay

`Project.load(dir, replay=False)` reads each component's blocks from its saved
grids and face graph and only records the journal: no march, no collar build.
Use it to open a large project to write outputs or make small edits. The
blocks, face BCs, journal and any Xall run written from it equal those of a
replayed load.

<<< @/../examples/cookbook/projects.py#no_replay_limits{python}

```text
refused: the history before this point was loaded without replay (Project.load(replay=False)): it cannot be undone.
refused: march: the current mesh (blocks) cannot be marched; ...
```

- The loaded history cannot be undone, and a loaded mesh cannot be marched
  again. Commands you run after the load can be undone as usual.
- Seams read back from grids are copies, not shared views; commands that edit
  nodes (`snap`, `redistribute`, `surface_grid`) write them through, so seams
  stay closed ([ADR 0040](/adr/0040-seams-kept-closed)).
- A project saved before the face graphs existed cannot be loaded without
  replay: load it once with replay and save it again.

## Edit a journal and replay it

<<< @/../examples/cookbook/projects.py#edit_journal{python}

```python
# SHORE journal: replay with `shore run <this file>`.
from shore.project import Project

p = Project()
p.load_geometry('sphere.stl')
p.build_topology('cubed_sphere', nk=9, ds=0.02, ni=33, nj=64, theta_cap=0.6108652381980153)
p.march()
p.redistribute('sub0', 'j_lo_k_lo', 'tanh2', beta=2.0)
p.write_geo('out/sphere')
```

- A journal is a plain Python script of `Project` calls, with defaults left
  out. Change a parameter, insert a command, point `load_geometry` at an
  edited file, and replay it with `run_journal` or `shore run`.
- Paths in a journal are used as written, relative to the directory you run
  from. The journal inside a saved project refers to `inputs/...`: run it from
  the project directory, or pass `run_journal(path, input_dir=project_dir)`.
- A journal is executed as Python: run only journals you trust.

```bash
shore run sphere_fine.py --save-project sphere_fine.shore
```

## Build a project from the CLI

`shore mesh`, `shore collar` and `shore primitive *` take `--project DIR
--component NAME`: they load the project (or create it), build that component
and save it. `--no-replay` opens the existing project without replay.

```bash
shore mesh sphere.stl --ni 33 --nj 64 --nk 9 --ds 0.02 --theta-cap-deg 35 \
    --project case.shore --component sphere
shore primitive box --min -4 -4 -4 --max 4 4 4 --ni 17 --nj 17 --nk 17 \
    --project case.shore --component background --no-replay
cd case.shore && shore run journal.py     # replays both components
```

Write the merged Xall files of all components with
`shore xall case.shore --out run`
([Multi-component Xall run](/cookbook/multi-component)).

## See also

- [Projects and journals](/concepts/projects)
- [`shore run`](/reference/run) and the command table there
- [`shore.project`](/reference/api/shore.project.project)
- [ADR 0004](/adr/0004-command-layer), [ADR 0016](/adr/0016-multi-component-projects),
  [ADR 0018](/adr/0018-cli-collars-and-projects), [ADR 0033](/adr/0033-load-without-replay)
