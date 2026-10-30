# Export CGNS, Plot3D and VTK

Xall reads `.grd` and `cc.par`; most other solvers and post-processors read
CGNS or Plot3D, and ParaView reads VTK. SHORE writes all three from a
project or a mesh.

The code is
[`examples/cookbook/exports.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/exports.py)
(CGNS and Plot3D, needs the `cgns` extra, about 2 s) and
[`examples/cookbook/exports_vtk.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/exports_vtk.py)
(VTK, needs the `vis` extra). Both start from a project with two
components, a cubed sphere and a background box:

<<< @/../examples/cookbook/exports.py#project{python}

## Write a CGNS file

::: code-group
```bash [CLI]
shore cgns case.shore -o case.cgns
shore cgns case.shore -c main -o sphere.cgns
```
<<< @/../examples/cookbook/exports.py#cgns{python} [Python]
:::

```text
base main: 6 zones, 24 1-to-1 connections, 12 BCs
base background: 1 zones, 0 1-to-1 connections, 6 BCs
```

| Block model | CGNS |
|---|---|
| component | `CGNSBase_t`, named after the component |
| block | structured `Zone_t`, named after the block label |
| SHARED / DIRICHLET / PERIODIC face | `GridConnectivity1to1_t` in both zones, with the transform from the face's `AxisMap` |
| WALL face | `BC_t` `BCWall` |
| SYMMETRY face | `BC_t` `BCSymmetryPlane` |
| FREE face | `BC_t` `BCTypeUserDefined`, tagged `shore:FREE` |

- Needs `h5py`: `pip install "shore-mesh[cgns]"`. The file is CGNS/HDF5,
  library version 3.4.
- `-c` / `components=` picks components (repeatable); the default is every
  component with a mesh.
- The cubed sphere has 4 seams per block, written once in each zone: 24
  connections. A face shared with several blocks (CH's `body.k_hi`) gets
  one connection per partner; nodes that match no partner are refused.
- FREE faces (outer boundaries, overlap boundaries) are user-defined:
  set the far-field or overset condition in your solver set-up.
- A box primitive keeps the block default `k_lo = WALL`, so the
  background box above has a `BCWall` at its `k_lo` face. Override it in
  the solver set-up when the box is a far field.

## Write a Plot3D grid

::: code-group
```bash [CLI]
shore plot3d case.shore -o case.xyz
```
<<< @/../examples/cookbook/exports.py#plot3d{python} [Python]
:::

```text
case.xyz: 7 blocks, dims [[33, 17, 17], [33, 17, 17], [33, 17, 17], [33, 17, 17], [17, 17, 17], [17, 17, 17], [33, 33, 33]]
```

- Multi-block whole 3-D grid: no IBLANK, Fortran unformatted (4-byte
  record markers), double precision, little endian; node counts in the
  header.
- Components are written in project order, blocks in component order.
  Plot3D carries no connectivity or boundary conditions; use CGNS when the
  reader needs them.

## Write VTK for ParaView

<<< @/../examples/cookbook/exports_vtk.py#vtk{python}

```text
[('sub0', (33, 17, 17)), ('sub1', (33, 17, 17)), ('sub2', (33, 17, 17)), ('sub3', (33, 17, 17)), ('cap_north', (17, 17, 17)), ('cap_south', (17, 17, 17))]
```

From `.geo` files:

```bash
shore export sphere_sub*.geo sphere_cap*.geo -o sphere.vtm
shore export sphere_sub0.geo -o sub0.vtm
```

- Needs `pyvista`: `pip install "shore-mesh[vis]"`.
- Cell data `Jacobian` is attached by default (`with_jacobian=False`,
  `--no-jacobian` to skip).
- A single `.geo` exported to a `.vts` is treated as a j-periodic O-grid
  block: a column is added that closes j. That is right for an `ogrid`
  block and wrong for any other; export other single blocks to `.vtm`
  (`-o name.vtm`), which keeps them open.

## See also

- [Block model](/concepts/block-model): faces, BCs and axis maps.
- [Projects and journals](/concepts/projects) and
  [Save, load, replay](/cookbook/projects).
- [ADR 0035](/adr/0035-cgns-plot3d-export): CGNS and Plot3D exports.
- CLI: [`shore cgns`](/reference/cgns), [`shore plot3d`](/reference/plot3d),
  [`shore export`](/reference/export).
- API: [`write_cgns`](/reference/api/shore.io.cgns#write_cgns),
  [`write_plot3d`](/reference/api/shore.io.plot3d#write_plot3d),
  [`Mesh.write_vtk`](/reference/api/shore.mesh.mesh#Mesh-write_vtk).
