# CLI commands

Every `shore` command and sub-command, grouped as in the sidebar. Each
page lists the options exactly as `shore <command> --help` prints them,
with their defaults. `shore --version` prints the installed version
(`shore 0.4.0`); `shore --help` and `shore <command> --help` list the
commands and options.

## Grids

| Command | Purpose | Extras |
|---|---|---|
| [`shore mesh`](mesh) | Full pipeline: STL or CAD → k=0 wall → marched volume → `.geo` files, for the `cubed_sphere`, `ogrid`, `cgrid`, `ch` and `oh` topologies. | `cad` for STEP/IGES/BREP |
| [`shore primitive box`](primitive#shore-primitive-box) | Cartesian box, one block. | |
| [`shore primitive annulus`](primitive#shore-primitive-annulus) | Annular cylinder, one block, optionally a j sector. | |
| [`shore primitive flat-caps`](primitive#shore-primitive-flat-caps) | 5-block cylinder without a central hole (Chimera background). | |
| [`shore primitive h-grid-rect`](primitive#shore-primitive-h-grid-rect) | Rectangular-duct H-grid, one block. | |
| [`shore primitive h-grid-tube`](primitive#shore-primitive-h-grid-tube) | Circular-bore H-grid, one block, j-periodic. | |
| [`shore collar`](collar) | Collar block at the junction of two bodies; half-collar on a symmetry plane; split at sharp corners. | `cad` for CAD bodies |
| [`shore cap-tip`](cap-tip) | Wing-tip cap blocks from a C-grid body `.geo`. | |
| [`shore project`](project) | Project a lat-lon sphere grid onto a body: the k=0 layer only, as `.npz`. | `cad` for CAD input |
| [`shore run`](run) | Replay a journal; optionally save the project it builds. | whatever the journal uses |
| [`shore gui`](gui) | The desktop GUI over a project (also `shore-gui`). | `gui` |

## Bodies

| Command | Purpose | Extras |
|---|---|---|
| [`shore profile naca`](profile) | NACA 4- or 5-digit section as a Selig `.dat`. | |
| [`shore body sphere`](body#shore-body-sphere) | Icosahedral sphere STL. | |
| [`shore body cylinder`](body#shore-body-cylinder) | Closed cylinder STL, axis z. | |
| [`shore body torus`](body#shore-body-torus) | Torus STL, axis z. | |
| [`shore body box`](body#shore-body-box) | Box STL. | |
| [`shore body capsule`](body#shore-body-capsule) | Capsule STL, axis z. | |
| [`shore body extrude`](body#shore-body-extrude) | Wing STL extruded from a Selig section (taper, sweep, twist). | |

## Xall and chimera

| Command | Purpose | Extras |
|---|---|---|
| [`shore xall`](xall) | The Xall run of every component of a project: per-component files, merged `.grd` and `cc.par`, one `proc.input`. | |
| [`shore grd`](grd) | Pack `.geo` files into a binary `.grd`. | |
| [`shore grd-merge`](grd-merge) | Concatenate `.grd` files of several components. | |
| [`shore cc-par`](cc-par) | `.adjacency.json` sidecar → `cc.par` (BCs, connection patches, marker boxes). | |
| [`shore cc-par-merge`](cc-par-merge) | Concatenate `cc.par` files, shifting block and patch indices. | |
| [`shore proc-input`](proc-input) | `proc.input` from a `.grd`, blocks balanced over MPI ranks. | |
| [`shore proc-input-merge`](proc-input#shore-proc-input-merge) | Concatenate `proc.input` files. | |
| [`shore boxes-from-stl`](boxes#generating-boxes-from-an-stl-shore-boxes-from-stl) | Hole-cutting marker boxes from a body STL, as JSON for `cc-par --boxes`. | |
| [`shore boxes-export`](boxes#visualising-boxes-in-paraview-shore-boxes-export) | Marker boxes → `.vtm` for ParaView. | `vis` |
| [`shore split`](split) | Split blocks along index planes, keeping the face graph. | |
| [`shore balance`](balance) | Plan splits for MPI load balance (`split.toml` or `bsplit.par`). | |

## Exports and inspection

| Command | Purpose | Extras |
|---|---|---|
| [`shore cgns`](cgns) | A project as one CGNS/HDF5 file. | `cgns` |
| [`shore plot3d`](plot3d) | A project's blocks as one multi-block Plot3D grid. | |
| [`shore export`](export) | `.geo` / `.npz` → VTK (`.vts`, `.vtm`) for ParaView. | `vis` |
| [`shore info`](info) | Shape and extent of a `.geo` or `.npz`. | |
| [`shore check`](check) | Cell Jacobians of a `.geo`, layer by layer. | |
| [`shore view`](view) | Interactive PyVista window on a `.geo` or `.npz`. | `vis` |

## Files

| Command | Purpose | Extras |
|---|---|---|
| [`shore config`](config#shore-config-subcommand) | `template`: print a default mesh config; `validate PATH`: check one. | |
| [File formats](formats) | `.geo`, `.grd`, `cc.par`, `proc.input`, `.adjacency.json`, `.npz`, project directory, CGNS, Plot3D. | |

## Options shared by the grid commands

| Option | Commands | Description |
|---|---|---|
| `--project DIR` | `mesh`, `collar`, `primitive *` | Add the result to the project in `DIR` (created if missing) and save it ([ADR 0018](/adr/0018-cli-collars-and-projects)). |
| `--component C` | `mesh`, `collar`, `primitive *` | The component of `--project` to build: default `main` (`collar` for `shore collar`). |
| `--no-replay` | `mesh`, `collar`, `primitive *`, `xall`, `cgns`, `plot3d`, `gui` | Open the project from its saved grids instead of replaying its journal ([ADR 0033](/adr/0033-load-without-replay)). |
| `--journal PATH` | `mesh`, `primitive *`, `cap-tip`, `split`, `grd`, `cc-par`, `proc-input` | Also write the run as a replayable journal for [`shore run`](run). |
| `--adjacency-json PATH` | `mesh`, `primitive *`, `collar`, `cap-tip` | Also write the `.adjacency.json` face graph, the input of `shore cc-par`. |
| `-c` / `--config PATH` | `mesh`, `primitive box/annulus/flat-caps`, `split` | Read parameters from a TOML file; explicit flags win ([Config files](config)). |

Parameter precedence everywhere: explicit flag → `--config` file →
built-in default. Errors print one line and exit with code 1.

## Extras

`pip install 'shore-mesh[extra]'`:

| Extra | Installs | Needed by |
|---|---|---|
| `cad` | `build123d` | STEP/IGES/BREP input to `mesh`, `collar`, `project`, `Project.load_geometry` |
| `gui` | PySide6-Essentials, pygfx, wgpu, rendercanvas | `shore gui`, `shore-gui` |
| `cgns` | `h5py` | `shore cgns`, `Project.write_cgns` |
| `vis` | `pyvista` | `export`, `view`, `boxes-export` |

## See also

- [Project scripts](run): the Python layer behind every grid command
- [Python API](/reference/api/)
- [Quick start](/guide/quickstart)
