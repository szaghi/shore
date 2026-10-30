---
title: shore.mesh.inspect
outline: [2, 3]
---

# `shore.mesh.inspect`

File-based inspection utilities: info, check, export, view.

These functions operate on already-written `.geo` / `.npz` / `.vtm`
files without requiring a live `Mesh` object in memory.  They mirror the
in-memory `Mesh` methods for offline or post-hoc inspection.

## `info_geo` <Badge type="info" text="function" /> {#info_geo}

```python
info_geo(path: str | Path) -> dict
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/inspect.py#L17">source</a>

Return metadata for a single-block `.geo` file.

**Parameters**

- **`path`** — Path to a `.geo` file.

**Returns**

- `dict with keys ``shape``, ``n_points``, ``n_cells``,`
- ```xyz_min``, ``xyz_max``.`

## `info_surface` <Badge type="info" text="function" /> {#info_surface}

```python
info_surface(path: str | Path) -> dict
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/inspect.py#L44">source</a>

Return metadata for a `.npz` surface file.

**Parameters**

- **`path`** — Path to a `.npz` file with a `'k0'` key.

**Returns**

- `dict with keys ``shape``, ``n_points``, ``xyz_min``, ``xyz_max``.`

## `check_geo` <Badge type="info" text="function" /> {#check_geo}

```python
check_geo(
    path: str | Path,
    warn_threshold: float = 1e-06,
    periodic_j: bool = False,
) -> dict
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/inspect.py#L74">source</a>

Check cell Jacobians layer-by-layer for a single-block `.geo` file.

**Parameters**

- **`path`** — Path to a `.geo` file.
- **`warn_threshold`** — Warn if jmin drops below this (but remains positive).
- **`periodic_j`** — Use the periodic-j Jacobian (`True`) or the open (non-periodic) Jacobian (`False`, default).  Cubed-sphere topology blocks are always open; pass `True` for single-block O-grid topology blocks.

**Returns**

- `dict with keys ``ok`` (bool), ``jmins`` (list of per-layer minimums),`
- ```any_negative``, ``any_warn``.`

## `export_geo` <Badge type="info" text="function" /> {#export_geo}

```python
export_geo(
    path: str | Path,
    output: str | Path | None = None,
    with_jacobian: bool = True,
) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/inspect.py#L120">source</a>

Export a `.geo` file to VTK StructuredGrid `.vts` for ParaView.

**Parameters**

- **`path`** — Input `.geo` file.
- **`output`** — Output `.vts` path.  Defaults to same stem with `.vts` extension.
- **`with_jacobian`** — Attach per-cell Jacobian as CellData.

**Returns**

- `Path` — The written `.vts` file path.

## `export_surface` <Badge type="info" text="function" /> {#export_surface}

```python
export_surface(path: str | Path, output: str | Path | None = None) -> Path
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/inspect.py#L151">source</a>

Export a `.npz` surface file to VTK StructuredGrid `.vts`.

**Parameters**

- **`path`** — Input `.npz` file with `'k0'` key.
- **`output`** — Output `.vts` path.  Defaults to same stem with `.vts` extension.

**Returns**

- `Path` — The written `.vts` file path.

## `view_geo_file` <Badge type="info" text="function" /> {#view_geo_file}

```python
view_geo_file(
    path: str | Path,
    *,
    colorby: str = 'Jacobian',
    show_edges: bool = True,
    opacity: float = 0.6,
    jupyter: bool = False,
    off_screen: bool = False,
) -> object
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/inspect.py#L184">source</a>

Open an interactive PyVista window for a `.geo` file.

**Parameters**

- **`path`** — Input `.geo` file.
- **`colorby`** — `"Jacobian"` (default) or `"none"`.
- **`show_edges`** — Draw cell edges.
- **`opacity`** — Grid opacity [0-1].
- **`jupyter`** — Use HTML backend for inline Jupyter display.
- **`off_screen`** — Render off-screen (for testing).

**Returns**

- `pyvista.Plotter`

## `view_surface_file` <Badge type="info" text="function" /> {#view_surface_file}

```python
view_surface_file(
    path: str | Path,
    *,
    colorby: str = 'radius',
    show_edges: bool = True,
    jupyter: bool = False,
    off_screen: bool = False,
) -> object
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/inspect.py#L228">source</a>

Open an interactive PyVista window for a `.npz` surface file.

**Parameters**

- **`path`** — Input `.npz` file with `'k0'` key.
- **`colorby`** — `"radius"` (default) or `"none"`.
- **`show_edges`** — Draw cell edges.
- **`jupyter`** — Use HTML backend for inline Jupyter display.
- **`off_screen`** — Render off-screen (for testing).

**Returns**

- `pyvista.Plotter`
