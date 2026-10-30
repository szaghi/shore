---
title: shore.io.vtk
outline: [2, 3]
---

# `shore.io.vtk`

VTK export for ParaView.

Writes structured grids to VTK XML StructuredGrid format (.vts), which
ParaView reads natively.

O-grid periodicity
------------------
VTK StructuredGrid has no concept of periodic connectivity. The standard
approach for O-grids is to close the j-direction by repeating the first
j-column at the end, producing a grid of shape (ni, nj+1, nk). This
creates a visible seam in ParaView that is geometrically correct.

Jacobian cell data
------------------
The full per-cell Jacobian field is attached as `CellData` under the
array name `"Jacobian"`, enabling direct quality visualisation in
ParaView without any post-processing filter.

## `geo_to_pyvista` <Badge type="info" text="function" /> {#geo_to_pyvista}

```python
geo_to_pyvista(grid: np.ndarray, with_jacobian: bool = True) -> pv.StructuredGrid
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/vtk.py#L77">source</a>

Convert a structured volume grid to a PyVista StructuredGrid.

The j-direction is closed (first column repeated at end) so that
the O-grid forms a topologically complete shell.

**Parameters**

- **`grid`** — Volume grid, shape (ni, nj, nk, 3).
- **`with_jacobian`** — If True (default), attach the per-cell Jacobian field as `CellData["Jacobian"]`.

**Returns**

- `pyvista.StructuredGrid`

## `surface_to_pyvista` <Badge type="info" text="function" /> {#surface_to_pyvista}

```python
surface_to_pyvista(k0: np.ndarray) -> pv.StructuredGrid
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/vtk.py#L110">source</a>

Convert a k=0 surface layer to a PyVista StructuredGrid (nk=1 slab).

**Parameters**

- **`k0`** — Surface layer, shape (ni, nj, 3).

**Returns**

- `pyvista.StructuredGrid`

## `write_vts` <Badge type="info" text="function" /> {#write_vts}

```python
write_vts(path: str | Path, grid: np.ndarray, with_jacobian: bool = True) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/vtk.py#L137">source</a>

Export a volume grid to VTK XML StructuredGrid format (.vts).

**Parameters**

- **`path`** — Output file path. The `.vts` extension is appended if absent.
- **`grid`** — Volume grid, shape (ni, nj, nk, 3).
- **`with_jacobian`** — Attach per-cell Jacobian as CellData (default True).

## `write_surface_vts` <Badge type="info" text="function" /> {#write_surface_vts}

```python
write_surface_vts(path: str | Path, k0: np.ndarray) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/vtk.py#L156">source</a>

Export a k=0 surface layer to VTK XML StructuredGrid format (.vts).

**Parameters**

- **`path`** — Output file path. The `.vts` extension is appended if absent.
- **`k0`** — Surface layer, shape (ni, nj, 3).

## `open_to_pyvista` <Badge type="info" text="function" /> {#open_to_pyvista}

```python
open_to_pyvista(grid: np.ndarray, with_jacobian: bool = True) -> pv.StructuredGrid
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/vtk.py#L173">source</a>

Convert a non-periodic (open) structured block to a PyVista StructuredGrid.

Unlike `geo_to_pyvista`, no j-direction closing is performed.  Used for
polar cap blocks and the equatorial sub-blocks of the cubed-sphere topology.

**Parameters**

- **`grid`** — Volume grid, shape (ni, nj, nk, 3).
- **`with_jacobian`** — If True (default), attach the per-cell Jacobian field as `CellData["Jacobian"]` using the open (non-periodic) Jacobian.

**Returns**

- `pyvista.StructuredGrid`

## `write_caps_vtm` <Badge type="info" text="function" /> {#write_caps_vtm}

```python
write_caps_vtm(
    path: str | Path,
    blocks: list[np.ndarray],
    with_jacobian: bool = True,
    labels: list[str] | None = None,
    boxes: list | None = None,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/vtk.py#L203">source</a>

Export a list of structured grids as a PyVista MultiBlock (.vtm).

Originally written for the 6-block cubed-sphere topology — when *labels*
is None and *blocks* has exactly 6 entries, the cubed-sphere block labels
(`sub0..sub3`, `cap_north`, `cap_south`) are used.  Otherwise
*labels* must be provided with one entry per block.

**Parameters**

- **`path`** — Output path. The `.vtm` extension is used if no extension given.
- **`blocks`** — List of volume grids, each shape `(ni, nj, nk, 3)`.
- **`with_jacobian`** — Attach per-cell Jacobian as CellData on each block.
- **`labels`** — Optional list of block labels (used as the keys of the resulting `MultiBlock`).  Must match `len(blocks)` when given.
- **`boxes`** — Optional list of :class:`shore.io.boxes.Box` instances (the `cc.par` marker boxes for this topology).  When provided, each box is appended as an additional MultiBlock child with label `"box_<n>"` (1-based) so ParaView shows the meshes and the marker boxes together.  See :func:`write_boxes_vtm` for the standalone form.

## `write_boxes_vtm` <Badge type="info" text="function" /> {#write_boxes_vtm}

```python
write_boxes_vtm(boxes: list, path: str | Path) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/vtk.py#L283">source</a>

Export `cc.par` marker boxes as a PyVista MultiBlock (.vtm).

Each :class:`shore.io.boxes.Box` becomes one `vtkUnstructuredGrid`
with a single `VTK_HEXAHEDRON` cell, exposed under the
MultiBlock key `"box_<n>"` (1-based).  Cell data fields:
`associated_block` and `box_index`.

**Parameters**

- **`boxes`** — Non-empty list of :class:`Box` instances.
- **`path`** — Output path; `.vtm` is used when no extension is given.

**Returns**

- `Path` — The path actually written.

**Raises**

- `ValueError` — On an empty `boxes` list.
