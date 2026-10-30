# shore cgns

Write a saved SHORE project as one CGNS/HDF5 file ([ADR 0035](/adr/0035-cgns-plot3d-export)):
one base per component, one structured zone per block, the 1-to-1
connections and boundary conditions of the face graph.

Needs `h5py`: `pip install 'shore-mesh[cgns]'`.

## Synopsis

```
shore cgns [OPTIONS] PROJECT_DIR
```

## Arguments and options

| Argument / option | Description |
|---|---|
| `PROJECT_DIR` | A saved project (`Project.save`, or `--project` of `shore mesh` / `shore collar` / `shore primitive`). |
| `-o`, `--output FILE` | The CGNS file to write. |
| `-c`, `--component NAME` | A component to write; repeatable. Default: every component with a mesh, in project order. |
| `--no-replay` | Load the project from its saved grids without replaying its journal ([ADR 0033](/adr/0033-load-without-replay)). |

## What the file holds

| SHORE | CGNS |
|---|---|
| component | `CGNSBase_t` named after it (cell and physical dimension 3) |
| block | structured `Zone_t` named by its label; `GridCoordinates` `CoordinateX/Y/Z`, double precision |
| SHARED / DIRICHLET / PERIODIC face | `GridConnectivity1to1_t` in both zones (`PointRange`, `PointRangeDonor`, `Transform` from the face's `AxisMap`) |
| WALL face | `BC_t` `BCWall` over the face |
| SYMMETRY face | `BC_t` `BCSymmetryPlane` |
| FREE face | `BC_t` `BCTypeUserDefined` with a `Descriptor_t` reading `shore:FREE` (an overset fringe or a natural boundary, decided by the solver setup) |

- 1-to-1 connections must be point-matched. A face shared with several
  smaller faces (CH's `body.k_hi` and its three channel blocks) gets one
  connection per partner, over the partner's sub-range.
- The O-grid's periodic j seam is not a face of the block graph (as in
  `cc.par`), so it is not written as a connection.
- Not written: overset hole cutting and donor data, flow solutions.

## Example

```bash
shore xall run/ --out xall/            # the Xall run
shore cgns run/ -o run.cgns            # the same grids for CGNS tools
shore cgns run/ -o collar.cgns -c collar --no-replay
```
