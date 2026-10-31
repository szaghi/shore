---
title: shore.mesh.mesh
outline: [2, 3]
---

# `shore.mesh.mesh`

Mesh: multi-block orchestrator — projection, marching, I/O.

## `check_blend_normals_k` <Badge type="info" text="function" /> {#check_blend_normals_k}

```python
check_blend_normals_k(blend_normals_k: int | None) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L59">source</a>

Refuse the removed `blend_normals_k` (ADR 0034); warn on an explicit 0.

## `Mesh` <Badge type="tip" text="class" /> {#Mesh}

```python
Mesh(body: Body | None, blocks: list[HexBlock], topology: Topology)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L260">source</a>

Multi-block structured mesh.

**Attributes**

- **`body`** — The triangulated body surface.
- **`blocks`** — Ordered list of blocks as returned by the topology builder.
- **`topology`** — The topology object that built the blocks.

### `Mesh.geometries` <Badge type="info" text="property" /> {#Mesh-geometries}

```python
geometries: dict[str, Geometry]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L291">source</a>

Geometry map for associations (`"body"` when the mesh has one).

### `Mesh.snap` <Badge type="info" text="method" /> {#Mesh-snap}

```python
snap(exact: bool = False) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L295">source</a>

Project every associated node onto its geometry (ADR 0003).

Faces, then edges, then vertices (vertex -&gt; edge -&gt; face priority).
`exact=True` uses the exact CAD geometry (tier 2) instead of its
tessellation; for STL bodies the two are the same.

### `Mesh.from_array` <Badge type="info" text="classmethod" /> {#Mesh-from_array}

```python
from_array(
    surface: np.ndarray,
    topology: Topology,
    nk: int,
    spacing_k: Spacing1D | None = None,
    **kwargs,
) -> Mesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L307">source</a>

Build a Mesh from a pre-built surface array (no STL).

Used by topologies whose input is a numeric array (e.g. a
cross-section polyline) rather than a triangulated body.  The
resulting :class:`Mesh` has `body = None`; methods that
depend on `self.body` (e.g. :meth:`snap`) will fail loudly if
called.

**Parameters**

- **`surface`** — Topology-specific surface array.  For :class:`CGridTipCap`: a `(ni, 3)` cross-section polyline.
- **`topology`** — The topology builder.
- **`nk`** — Topology-specific node count.  For :class:`CGridTipCap`, the cap's spanwise (k-axis) node count.
- **`spacing_k`** — Optional wall-normal / spanwise spacing law forwarded to `topology.build`.
- **`**kwargs`** — Forwarded verbatim to `topology.build`.

**Returns**

- `Mesh` — With `body=None` and the topology-built blocks.

### `Mesh.from_stl` <Badge type="info" text="classmethod" /> {#Mesh-from_stl}

```python
from_stl(path: str | Path, topology: Topology, *args: object, **kwargs: object) -> Mesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L353">source</a>

Alias of `from_file` (which also reads STEP/IGES/BREP).

### `Mesh.from_file` <Badge type="info" text="classmethod" /> {#Mesh-from_file}

```python
from_file(
    path: str | Path,
    topology: Topology,
    ni: int,
    nj: int,
    nk: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
    spacing_k: Spacing1D | None = None,
) -> Mesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L360">source</a>

Load a body, project to surface, wire blocks.  Does not march.

**Parameters**

- **`path`** — STL, or a CAD file (STEP/IGES/BREP, needs `shore-mesh[cad]`). For CAD the k=0 ray casts are snapped to the exact surface.
- **`topology`** — Wired topology builder (e.g. `CubedSphereTopology()`).
- **`ni, nj`** — Surface node counts (i and j directions, the array shape).
- **`nk`** — Wall-normal node count, including the wall row at k=0.  With `nk = 30` the marcher extrudes 29 cells from the wall. (SHORE convention: `ni`, `nj`, `nk` are always node counts; the only file that stores cell counts is the binary `.grd` header, where the conversion is automatic.)
- **`theta_cap / theta_cap_deg`** — Polar exclusion angle.
- **`spacing_k`** — Wall-normal spacing law.  Defaults to `Spacing1D(law='geometric', ds=1e-3)`.

**Returns**

- `Mesh` — Ready to call `.march()`.

### `Mesh.march` <Badge type="info" text="method" /> {#Mesh-march}

```python
march(
    ds: float = 0.001,
    growth: float = 1.15,
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    steps: np.ndarray | None = None,
    layer_callback: Callable[[int, float], None] | None = None,
    blend_normals_k: int | None = None,
    smooth_mode: str = 'limited',
    march_scheme: str = 'algebraic',
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L474">source</a>

Extrude all blocks simultaneously by hyperbolic marching.

Fills `block.nodes[:, :, 1:, :]` in place for every block.

**Parameters**

- **`ds`** — First-layer thickness (geometric law).
- **`growth`** — Geometric growth ratio (geometric law).
- **`smooth_strength`** — Laplacian blend factor per layer.
- **`smooth_iters`** — Number of Laplacian sweeps per layer.
- **`steps`** — Pre-computed step sizes `(nk-1,)`.  Overrides `ds`/`growth`.
- **`layer_callback`** — Optional `(k, jmin_global)` called after each accepted layer.
- **`blend_normals_k`** — **Removed** (ADR 0034): only `None` (the default) or 0 are accepted; 0 warns. It blended body normals into the first layers' marching directions, with no benefit (ADR 0019).
- **`smooth_mode`** — `'limited'` (default): the Laplacian smoothing motion is limited per node where it would damage the layer (it may pull a node back by at most a quarter of the step, and move it sideways by at most a quarter of the step); elsewhere it is unchanged, so grids that marched cleanly before are identical. `'legacy'`: never limit (kept to reproduce old grids); on curved walls the layer shrinks inward by about `smooth_strength * h**2 * curvature` and folds layers thinner than that.
- **`march_scheme`** — `'algebraic'` (default): extrude along the local normal, then Laplacian smoothing (`smooth_*`). `'hyperbolic'`: the Steger-Chaussee scheme (`volume.hyperbolic_scheme`): k-lines orthogonal to each layer, prescribed cell volumes, implicit coupling and dissipation that grow with step / cell size; the `smooth_*` arguments are not used.

**Raises**

- `MarchingError` — If a negative-Jacobian cell is detected.
- `NotImplementedError` — If the topology has no registered marcher.
- `ParameterError` — `blend_normals_k` other than 0 (removed, ADR 0034).

**Warns**

- `DeprecationWarning` — `blend_normals_k=0` given explicitly.

### `Mesh.write_geo` <Badge type="info" text="method" /> {#Mesh-write_geo}

```python
write_geo(stem: str | Path) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L1069">source</a>

Write each block as a separate `<stem>_<label>.geo` file.

### `Mesh.write_vtk` <Badge type="info" text="method" /> {#Mesh-write_vtk}

```python
write_vtk(stem: str | Path, with_jacobian: bool = True) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L1077">source</a>

Write all blocks as a single PyVista MultiBlock `.vtm` file.

**Parameters**

- **`stem`** — Output path (without extension).  Writes `<stem>.vtm`.
- **`with_jacobian`** — Attach per-cell Jacobian as CellData (default True).

### `Mesh.info` <Badge type="info" text="method" /> {#Mesh-info}

```python
info() -> dict
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L1101">source</a>

Return a dict of mesh metadata.

**Returns**

- `dict with keys:` — `n_blocks`, `labels`, and per-block dicts under each label with `shape`, `n_points`, `n_cells`, `xyz_min`, `xyz_max`.

### `Mesh.check` <Badge type="info" text="method" /> {#Mesh-check}

```python
check(warn_threshold: float = 1e-06) -> dict
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L1125">source</a>

Check cell Jacobians for all blocks layer by layer.

**Parameters**

- **`warn_threshold`** — Warn if jmin drops below this value (but stay positive).

**Returns**

- `dict with keys ``ok`` (bool) and ``blocks`` (per-block layer data).`

**Raises**

- `MarchingError` — If any block has a negative Jacobian layer.

### `Mesh.export` <Badge type="info" text="method" /> {#Mesh-export}

```python
export(stem: str | Path, with_jacobian: bool = True) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L1164">source</a>

Export the mesh to VTK MultiBlock format (.vtm) for ParaView.

Equivalent to :meth:`write_vtk`.

**Parameters**

- **`stem`** — Output path stem (extension set to `.vtm` automatically).
- **`with_jacobian`** — Attach per-cell Jacobian as CellData.

### `Mesh.view` <Badge type="info" text="method" /> {#Mesh-view}

```python
view(
    *,
    jupyter: bool = False,
    off_screen: bool = False,
    show_edges: bool = True,
    opacity: float = 0.6,
    colorby: str = 'Jacobian',
) -> object
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/mesh.py#L1182">source</a>

Open an interactive PyVista window showing all blocks.

**Parameters**

- **`jupyter`** — Use HTML backend for inline Jupyter display.
- **`off_screen`** — Render off-screen (for testing).
- **`show_edges`** — Draw cell edges.
- **`opacity`** — Block opacity [0-1].
- **`colorby`** — Cell data field to colour by (`"Jacobian"` or `"none"`).

**Returns**

- `pyvista.Plotter`
