---
title: shore.surface.tfi
outline: [2, 3]
---

# `shore.surface.tfi`

Surface grids on one parametric face: TFI and elliptic smoothing in (u, v) (ADR 0007).

`surface_grid` fills a face from its four boundary rows:

1. the boundary nodes are projected once onto the SURFACE entity, which
   gives their parameters (u, v);
2. transfinite interpolation (TFI) of those parameters;
3. with `smoothing="elliptic"`, Winslow smoothing in (u, v):

       alpha (r_xixi + P r_xi) - 2 beta r_xieta + gamma (r_etaeta + Q r_eta) + c = 0

   for `r = (u, v)`. `alpha = |x_eta|^2`, `beta = x_xi . x_eta` and
   `gamma = |x_xi|^2` are physical (`x_xi = S_u u_xi + S_v v_xi`), and
   ``c = G^-1 [S_u, S_v]^T (alpha S(r_xi, r_xi) - 2 beta S(r_xi, r_eta) +
   gamma S(r_eta, r_eta))`` holds the parametrisation's second derivatives
   (`S(p, w) = S_uu p_u w_u + S_uv (p_u w_v + p_v w_u) + S_vv p_v w_v`,
   `G` the first fundamental form): the physical Winslow equations
   projected on the tangent plane, so the result does not depend on how
   the face is parametrised. `P` and `Q` are Thomas-Middlecoff control
   functions, `-(x_xi . x_xixi) / |x_xi|^2` along the two j-boundaries
   (`P`) and the two i-boundaries (`Q`) of the physical boundary,
   blended linearly across the face: on tensor-product boundary data they
   make the TFI grid an exact solution, so the boundary clustering is kept;
4. one evaluation on the surface at the end; the boundary rows are the
   input.

The surface metric comes from a parameter lattice evaluated once (its
finite-difference derivatives are interpolated bilinearly), so the number
of geometry queries is bounded per node, never per sweep (ADR 0002). The
nonlinear system is solved by Picard iteration, one SuperLU factorisation
per outer step for both u and v (deterministic).

## `EllipticParams` <Badge type="tip" text="dataclass" /> {#EllipticParams}

```python
EllipticParams(
    lattice_refine: int = 2,
    max_iterations: int = 400,
    tolerance: float = 1e-12,
    control_limit: float = 1.5,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/tfi.py#L49">source</a>

Defaults of the elliptic smoother (not per-case tuning knobs).

## `SurfaceGrid` <Badge type="tip" text="dataclass" /> {#SurfaceGrid}

```python
SurfaceGrid(points: np.ndarray, params: np.ndarray)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/tfi.py#L62">source</a>

A structured grid on one face.

**Parameters**

- **`points`** (np.ndarray, shape (ni, nj, 3)) — Node positions, on the surface; the boundary rows are the input.
- **`params`** (np.ndarray, shape (ni, nj, 2)) — The nodes' (u, v) on the surface.

## `surface_grid` <Badge type="info" text="function" /> {#surface_grid}

```python
surface_grid(
    geometry: Geometry,
    surface: EntityRef,
    boundary: np.ndarray,
    smoothing: str = 'elliptic',
    params: EllipticParams = DEFAULT,
) -> SurfaceGrid
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/tfi.py#L77">source</a>

Fill the face of *geometry*'s SURFACE entity *surface* from *boundary*.

**Parameters**

- **`geometry`** — A backend with parametric SURFACE entities (`OcctGeometry`).
- **`surface`** — The SURFACE entity the face lies on.
- **`boundary`** (np.ndarray, shape (ni, nj, 3)) — Only its four boundary rows are read; they must lie on the surface.
- **`smoothing`** — `"none"` (plain TFI of the parameters) or `"elliptic"`.
- **`params`** — Elliptic smoother settings.

**Returns**

- `SurfaceGrid`

**Raises**

- `ParameterError` — Unknown *smoothing*; *boundary* not `(ni, nj, 3)` with ni, nj &gt;= 3; boundary nodes off the surface (beyond 1e-6 of the model diagonal); an entity without (u, v) parameters.
