# shore plot3d

Write the blocks of a saved SHORE project as one multi-block Plot3D grid
file ([ADR 0035](../adr/0035-cgns-plot3d-export.md)).

## Synopsis

```
shore plot3d [OPTIONS] PROJECT_DIR
```

## Arguments and options

| Argument / option | Description |
|---|---|
| `PROJECT_DIR` | A saved project. |
| `-o`, `--output FILE` | The Plot3D grid file to write (conventionally `.xyz` or `.x`). |
| `-c`, `--component NAME` | A component to write; repeatable. Default: every component with a mesh, in project order. |
| `--no-replay` | Load the project from its saved grids without replaying its journal ([ADR 0033](../adr/0033-load-without-replay.md)). |

## Format

Multi-block, whole (no IBLANK), 3-D grid; Fortran unformatted sequential
records, little endian, double precision:

| Record | Content |
|---|---|
| 1 | number of blocks (`int32`) |
| 2 | `ni, nj, nk` of every block (`int32`, node counts) |
| 2 + n | `X`, `Y`, `Z` of block n, each i-fastest (`float64`) |

Blocks follow the components in order, and the blocks of a component in
their order. Plot3D carries no names, connectivity or boundary conditions:
use [`shore cgns`](./cgns) for those.
