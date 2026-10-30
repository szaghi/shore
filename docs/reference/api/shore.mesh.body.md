---
title: shore.mesh.body
outline: [2, 3]
---

# `shore.mesh.body`

Body: a geometry (STL or CAD) with structured sphere projection.

## `Body` <Badge type="tip" text="dataclass" /> {#Body}

```python
Body(
    mesh: trimesh.Trimesh | None = None,
    anchor: np.ndarray = field(default=None, compare=False, repr=False),
    geometry: Geometry = field(default=None, compare=False, repr=False),
    anchor_error: MeshInputError | None = field(default=None, compare=False, repr=False),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L26">source</a>

Body surface with ray-cast projection capability.

**Parameters**

- **`mesh`** — Triangulated body surface. For a CAD body (`geometry` given) it may be omitted and defaults to the geometry's merged tessellation, which serves the mesh-based consumers (anchor resolution, C-grid slicing).
- **`anchor`** — Ray-cast origin, shape (3,).  When *None*, defaults to the mesh centroid (set by `__post_init__`).
- **`geometry`** — The body as a `Geometry`. When *None*, an `StlGeometry` over *mesh*. Projections (k=0 ray casts, marching normals) go through it.
- **`anchor_error`** — Why no anchor could be resolved (`from_stl` / `from_file`): the body loads, and the error is raised where the anchor is first used, by `projector` (ADR 0018). Collars never ray-cast, so an open half-model body without an interior point still serves them.

### `Body.from_stl` <Badge type="info" text="classmethod" /> {#Body-from_stl}

```python
from_stl(
    path: str | Path,
    anchor_strategy: str = 'auto',
    skip_hygiene: bool = False,
) -> Body
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L69">source</a>

Load an STL and resolve the anchor point.

**Parameters**

- **`path`** — Path to the STL file.
- **`anchor_strategy`** — Anchor resolution strategy passed to `shore.surface.anchor.resolve_anchor`. `'auto'` (default) tries incenter then bbox with containment validation; `'centroid'` replicates the legacy behaviour.
- **`skip_hygiene`** — Bypass mesh hygiene (merge_vertices / fix_normals / remove degenerate faces).  Preserves legacy behaviour for regression tests.

### `Body.from_file` <Badge type="info" text="classmethod" /> {#Body-from_file}

```python
from_file(
    path: str | Path,
    anchor_strategy: str = 'auto',
    skip_hygiene: bool = False,
) -> Body
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L97">source</a>

Load an STL or a CAD file (STEP/IGES/BREP) and resolve the anchor.

STL goes through `from_stl`. CAD needs `shore-mesh[cad]`; the
anchor is resolved on the merged tessellation, and *skip_hygiene*
does not apply.

### `Body.is_cad` <Badge type="info" text="property" /> {#Body-is_cad}

```python
is_cad: bool
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L124">source</a>

True when the body has exact (non-STL) geometry.

### `Body.projector` <Badge type="info" text="property" /> {#Body-projector}

```python
projector: Projector | AnchoredCaster
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L129">source</a>

Anchored ray caster for the k=0 builders, one per Body.

STL: `geometry.projector(anchor)` (the legacy `Projector`, sharing
the geometry's BVH; results bit-identical to before). CAD: an
`AnchoredCaster` that casts on the tessellation and snaps hits to
the exact surface.

Raises the anchor resolution error when the body has no anchor.

### `Body.normals` <Badge type="info" text="method" /> {#Body-normals}

```python
normals(points: np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L149">source</a>

Unit outward surface normals at the closest body points (..., 3).

A CAD body gives its exact geometry normals. An STL body gives the
blend of trimesh's vertex normals (`_stl_normals`), smooth across
sharp edges, not the crease-aware normals of `StlGeometry` (ADR
0014, 0015). ADR 0005's acceptance spec measures k-line angles against
them; the marcher's normal blending that also read them is removed
(ADR 0034).

### `Body.project` <Badge type="info" text="method" /> {#Body-project}

```python
project(
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L168">source</a>

Ray-cast a structured lat-lon grid onto the body surface.

**Parameters**

- **`ni, nj`** — Latitude and longitude resolution.
- **`theta_cap`** — Polar exclusion angle in radians.
- **`theta_cap_deg`** — Polar exclusion angle in degrees (mutually exclusive with theta_cap).

**Returns**

- `np.ndarray, shape (ni, nj, 3)` — Coordinates of the k=0 structured surface layer.

### `Body.project_caps` <Badge type="info" text="method" /> {#Body-project_caps}

```python
project_caps(
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
    pole_eps_deg: float = 0.5,
) -> tuple[np.ndarray, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L212">source</a>

Unified ray-cast covering the full sphere for the cubed-sphere topology.

Returns the full projected surface plus the split index `nc_cap`.
Rows `0 .. nc_cap-1` are the south cap, rows `nc_cap .. ni_total-nc_cap`
are the equatorial band, rows `ni_total-nc_cap ..` are the north cap.

**Parameters**

- **`ni`** — Number of equatorial rows.
- **`nj`** — Longitude resolution (multiple of 4).
- **`theta_cap`** — Polar split angle in radians.  Default `0.05 * pi` (~9 deg).
- **`theta_cap_deg`** — Convenience alias in degrees.
- **`pole_eps_deg`** — Angular exclusion at each pole in degrees (default 0.5 deg).

**Returns**

- **`surface`** (np.ndarray, shape (ni_total, nj, 3))
- **`nc_cap`** (int)

### `Body.directions` <Badge type="info" text="method" /> {#Body-directions}

```python
directions(
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L264">source</a>

Return (ni, nj, 3) unit direction vectors without casting rays.

### `Body.save_surface` <Badge type="info" text="staticmethod" /> {#Body-save_surface}

```python
save_surface(surface: np.ndarray, path: str | Path) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L294">source</a>

Save a k=0 surface layer to a compressed .npz file.

### `Body.load_surface` <Badge type="info" text="staticmethod" /> {#Body-load_surface}

```python
load_surface(path: str | Path) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/mesh/body.py#L302">source</a>

Load a k=0 surface layer from a .npz file.
