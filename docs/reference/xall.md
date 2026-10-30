# shore xall

Write the Xall run of every component of a SHORE project: per component,
its `.geo` files, adjacency sidecar, `.grd` and `cc.par`; then the merged
`.grd` and `cc.par`, and one `proc.input` balanced over the whole run
([ADR 0016](/adr/0016-multi-component-projects), [ADR 0017](/adr/0017-xall-options-per-component), [ADR 0018](/adr/0018-cli-collars-and-projects)).

## Synopsis

```
shore xall [OPTIONS] PROJECT_DIR
```

## Arguments and options

| Argument / option | Description |
|---|---|
| `PROJECT_DIR` | A saved project (`Project.save`, or `--project` of `shore mesh` / `shore collar` / `shore primitive`). |
| `--out DIR` (required) | Directory for the run files: `DIR/<component>/...`, then `DIR/<grd>`, `DIR/cc.par`, `DIR/proc.input`. |
| `--options FILE` | TOML file with the run and per-component options (below). |
| `--no-replay` | Open the project from its saved grids without replaying its journal ([ADR 0033](/adr/0033-load-without-replay)). |

## Options file

```toml
[run]                       # one per run: merge_cc_par needs equal headers
grd = "run.grd"             # merged grid file name (default mesh.grd)
np = 64                     # MPI ranks of proc.input (default 1)
output_basename = "cc."
multigrid_levels = [4, 1]
boundary_layer_thickness = -0.1
numerical_beach = -1.0

[components.background]
free_face_bc = 9            # natural BC of its FREE faces (9: extrapolation)
body = 0                    # Xall body number of all its blocks

[components.main]
body = 1
boxes = "boxes.json"        # marker boxes (shore boxes-from-stl), relative to this file

[components.main.meta.ogrid]          # per-block cc.par metadata
level = 2
priority = 3
comment = "hull"

[components.collar.proc_meta.collar]  # per-block proc.input overrides
proc = 1
```

An unknown table or key, a component that does not exist or has no mesh,
and a block label the component does not have are errors (exit code 1),
reported before anything is written.
