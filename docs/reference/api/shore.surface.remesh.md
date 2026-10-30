---
title: shore.surface.remesh
outline: [2, 3]
---

# `shore.surface.remesh`

Unstructured STL surface → structured (ni, nj) quad grid via sphere projection.

Strategy
--------
1. Place a unit sphere centered at the body anchor point.
2. Build a structured (ni, nj) lat-lon grid on that sphere.
3. For each sphere point, cast a ray from the anchor outward and find
   the first intersection with the STL body.
4. The intersection points form the k=0 layer of the O-grid.

For the 6-block cubed-sphere topology a unified projection covers the full
sphere in one ray-cast (`project_caps_to_surface`).  The polar cap blocks
and equatorial band share the same `di` by construction — no Laplacian
seed needed.

Limitation: the approach requires that every ray from the anchor hits the
body exactly once (i.e. the body is star-shaped w.r.t. its anchor).  This
covers spheres, ellipsoids, fuselages, and most convex bodies.  Concave or
multiply-connected surfaces need block decomposition — deferred to future work.

## `sphere_latlon_directions` <Badge type="info" text="function" /> {#sphere_latlon_directions}

```python
sphere_latlon_directions(
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    theta_margin: float | None = None,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/remesh.py#L31">source</a>

Build (ni, nj) unit direction vectors on a lat-lon sphere.

The poles are excluded so the grid is open (no singular points).
The grid closes in the j-direction (periodic), not in i.

**Parameters**

- **`ni`** — Number of points in the latitude direction (pole-to-pole, exclusive).
- **`nj`** — Number of points in the longitude direction (full circle, periodic).
- **`theta_cap`** — Polar exclusion angle in radians: the equatorial band covers `[theta_cap, pi - theta_cap]`.  The cap boundary parallel has radius `sin(theta_cap)` on the unit sphere.  Must be in `(0, pi/2)`.  Default `0.05*pi` (~9°).
- **`theta_cap_deg`** — Convenience alias for `theta_cap` in degrees.  Mutually exclusive with `theta_cap` and `theta_margin`.
- **`theta_margin`** — Deprecated alias for `theta_cap` (radians).

**Returns**

- `numpy.ndarray, shape (ni, nj, 3)` — Unit direction vectors.

**Raises**

- `ParameterError` — If more than one of `theta_cap`, `theta_cap_deg`, `theta_margin` are supplied, or if the resolved angle is outside `(0, pi/2)`.

## `project_sphere_to_surface` <Badge type="info" text="function" /> {#project_sphere_to_surface}

```python
project_sphere_to_surface(
    mesh: trimesh.Trimesh,
    ni: int,
    nj: int,
    theta_cap: float | None = None,
    theta_cap_deg: float | None = None,
    theta_margin: float | None = None,
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
    j_cluster_lon: float | None = None,
    projector: object | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/remesh.py#L106">source</a>

Project a structured sphere grid onto a triangulated surface.

**Parameters**

- **`mesh`** — The body surface (raw STL).
- **`ni`** — Latitude resolution of the structured grid.
- **`nj`** — Longitude resolution of the structured grid.
- **`theta_cap`** — Polar exclusion angle in radians: the equatorial band covers `[theta_cap, pi - theta_cap]`.  Default `0.05*pi` (~9°).
- **`theta_cap_deg`** — Convenience alias for `theta_cap` in degrees.
- **`theta_margin`** — Deprecated alias for `theta_cap` (radians).
- **`projector`** — Optional pre-built `Projector` instance (shares its cached BVH). When *None*, a temporary projector using `mesh.centroid` is created.

**Returns**

- `numpy.ndarray, shape (ni, nj, 3)` — Coordinates of the k=0 structured surface layer.

**Raises**

- `ParameterError` — If mutually exclusive theta_cap arguments are combined.
- `MeshInputError` — If any ray misses the body.

## `project_caps_to_surface` <Badge type="info" text="function" /> {#project_caps_to_surface}

```python
project_caps_to_surface(
    mesh: trimesh.Trimesh,
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
    projector: object | None = None,
) -> tuple[np.ndarray, int]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/remesh.py#L173">source</a>

Project a unified lat-lon grid (cap blocks + equatorial) onto a surface.

Unlike `project_sphere_to_surface`, this function covers the full sphere
in a single ray-cast.  The equatorial band and both polar caps share the
same angular step `di` by construction — the seam between equatorial and
cap blocks is spacing-continuous.

**Parameters**

- **`mesh`** — Triangulated body surface.
- **`ni`** — Number of rows in the equatorial band.
- **`nj`** — Number of longitude points (multiple of 4 for the cubed-sphere topology).
- **`theta_cap`** — Requested polar exclusion angle in radians.  Default `0.05 * pi`.
- **`theta_cap_deg`** — Convenience alias in degrees.  Mutually exclusive with `theta_cap`.
- **`i_spacing`** — Latitude spacing law: `'uniform'` (default), `'tanh'`, `'tanh2'`.
- **`i_beta`** — Clustering strength for tanh/tanh2 i-spacing.
- **`j_spacing`** — Longitude spacing law: `'uniform'` or `'tanh2'`.
- **`j_beta`** — Clustering strength for tanh2 j-spacing.
- **`j_cluster_lon`** — Target longitude for tanh2 j-clustering (radians).
- **`pole_eps_deg`** — Angular exclusion at each pole in degrees (default 0.5°).
- **`projector`** — Optional pre-built `Projector` (shares cached BVH).

**Returns**

- **`surface`** (numpy.ndarray, shape (ni_total, nj, 3))
- **`nc_cap`** (int) — Number of rows in each polar cap block.

**Raises**

- `ParameterError` — If theta_cap and theta_cap_deg are both supplied.
- `MeshInputError` — If any ray misses the body.
