# Write a multi-component Xall run

An overset case is several independent grids, a near-body grid per body
and a background, that Xall reads as one block list. SHORE keeps each grid
as a component of one project and writes the merged run in one step; the
`*-merge` commands do the same from separate files.

The code is
[`examples/cookbook/multi_component.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/multi_component.py)
(about 10 s, most of it the marker boxes).

## Build the components in one project

<<< @/../examples/cookbook/multi_component.py#components{python}

```text
{'main': 6, 'background': 5}
```

- `use_component(name)` makes `name` current, creating it if new; every
  build command acts on the current component. The first one is `main`.
- From the CLI, `--project DIR --component NAME` on `shore mesh`,
  `shore collar` and `shore primitive ...` adds a component to a saved
  project (see [Save, load, replay](/cookbook/projects)).
- Each component keeps its own k-schedule and its own block labels; the
  labels only have to be unique within a component.

## Write the run with an options file

The options file sets run-wide `cc.par` header values under `[run]` and
per-component values under `[components.<name>]`:

<<< @/../examples/cookbook/multi_component.py#options{python}

```bash
shore xall case.shore --out run --options xall.toml
```

```text
✓ Wrote the Xall run of ['main', 'background'] to run
```

The Python call takes the same options as keyword arguments; it writes the
same files byte for byte:

<<< @/../examples/cookbook/multi_component.py#write{python}

```text
['run/background', 'run/cc.par', 'run/main', 'run/proc.input', 'run/run.grd']
```

| Key | Where | Default | Meaning |
|---|---|---|---|
| `grd` | `[run]` | `mesh.grd` | merged grid file name |
| `np` | `[run]` | 1 | MPI ranks of `proc.input`, balanced over the merged grid |
| `output_basename`, `multigrid_levels`, `boundary_layer_thickness`, `numerical_beach` | `[run]` | `cc.`, `[4, 1]`, `-0.1`, `-1.0` | `cc.par` header, one for the run |
| `free_face_bc` | component | 40 | `cc.par` code of its FREE faces: 40 chimera, a natural BC (1–19) for a far field, e.g. 9 extrapolation |
| `body` | component | 0 | Xall body number of all its blocks in `proc.input` |
| `boxes` | component | none | marker boxes, a JSON file from `shore boxes-from-stl` (relative to the options file) |
| `meta.<block>` | component | | per-block `level`, `group`, `priority`, `comment`, `free_face_bc` |
| `proc_meta.<block>` | component | | per-block `proc`, `group`, `body` in `proc.input` |

- `run/<component>/` holds each component's `.geo` files, adjacency
  sidecar, `.grd` and `cc.par`; `run/` holds the merged `run.grd`,
  `cc.par` and one `proc.input`.
- Components are merged in creation order. Unknown keys, unknown
  components or block labels, and components without a mesh are errors,
  reported before anything is written.
- `shore xall` requires `--out`.

Background primitives (`box`, `annulus`, `flat_caps`) have no wall: every
face that is not a seam is FREE, so `free_face_bc` sets all their boundary
codes ([ADR 0055](/adr/0055-background-primitives-have-no-wall)).

## Merge separate component files

When the components come from separate runs, merge their files with the
same input order in all three commands:

```bash
shore grd-merge main.grd background.grd -o run.grd
shore cc-par-merge main.cc.par background.cc.par -o cc.par --grd run.grd
shore proc-input run.grd proc.input --np 8
```

```text
✓ Wrote run.grd (2 inputs merged)
✓ Wrote cc.par (2 inputs merged)
✓ Wrote proc.input
...
  imbalance: 0.0%
```

- `grd-merge` concatenates the block payloads byte for byte and rewrites
  the block count. The result here is identical to `run/run.grd`.
- `cc-par-merge` shifts block and patch indices of the later inputs and
  prefixes each block comment with its input's stem (`main.sub0`,
  `background.center`). Header fields come from the first input; a
  disagreement is a warning.
- `shore proc-input` on the merged `.grd` balances all blocks together.
  `shore proc-input-merge` instead concatenates per-component
  `proc.input` files without re-balancing: two components each balanced
  on their own ranks pile onto the same low ranks (150 % imbalance for
  this case, 6 + 2 ranks merged). Use it only when you assigned the ranks
  yourself.

## See also

- [Chimera grids](/concepts/chimera) and [Marker boxes](/concepts/marker-boxes).
- [Chimera assembly](/tutorials/chimera-assembly) and
  [Xall bundle](/tutorials/xall-bundle): the same workflow told end to end.
- [Split and balance](/cookbook/split-balance): more blocks for more ranks.
- [ADR 0016](/adr/0016-multi-component-projects),
  [ADR 0017](/adr/0017-xall-options-per-component),
  [ADR 0018](/adr/0018-cli-collars-and-projects).
- CLI: [`shore xall`](/reference/xall), [`shore grd-merge`](/reference/grd-merge),
  [`shore cc-par-merge`](/reference/cc-par-merge),
  [`shore proc-input-merge`](/reference/proc-input#shore-proc-input-merge),
  [`shore boxes-from-stl`](/reference/boxes).
- API: [`Project.write_xall`](/reference/api/shore.project.project#Project-write_xall).
