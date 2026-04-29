# Spacing laws

SHORE provides three 1D spacing distributions for controlling node density along each structured direction. These are implemented in `shore.volume.spacing` using pure numpy — no scipy dependency.

## The three laws

### Uniform

Equal spacing across the interval. This is the default for all directions.

$$x_i = x_0 + i \cdot \Delta x, \quad \Delta x = \frac{x_N - x_0}{N - 1}$$

Use when: the geometry is smooth and you have no reason to cluster cells in any particular region.

### Tanh (one-sided Roberts)

Clusters points toward the **start** of the interval (the wall in the k-direction, the cap boundary in the i-direction). The distribution follows the reversed Roberts formula:

$$t(\tau) = 1 - \frac{\tanh\bigl(\beta(1 - \tau)\bigr)}{\tanh(\beta)}$$

where $\tau \in [0, 1]$ is the normalised parameter and $\beta > 0$ controls clustering strength.

- $\beta = 1$: mild clustering
- $\beta = 3$: moderate (default)
- $\beta = 5$: aggressive

Use when: you need fine resolution at one boundary (wall, or cap edge) and are willing to accept coarser cells at the opposite end.

### Tanh2 (two-sided Roberts)

Clusters points toward **both ends** of the interval, or toward a specified interior point:

$$t(\tau) = \frac{1}{2}\left(1 + \frac{\tanh\bigl(\beta(\tau - \tau_c)\bigr)}{\tanh(\beta \tau_c / 2)}\right)$$

where $\tau_c$ is the normalised cluster target:
- $\tau_c = 0.5$: symmetric clustering at both ends (default when no target specified)
- $\tau_c \neq 0.5$: asymmetric clustering shifted toward that point

Use when: you need fine resolution at both boundaries, or at a specific interior location (e.g., the leading edge of an airfoil).

### Pinned

Returns a 1D distribution whose **first and last intervals exactly equal caller-supplied step sizes**, with the interior shaped by one of the three preceding laws (`uniform` / `tanh` / `tanh2`):

$$x_1 - x_0 = \Delta s_{\text{first}}, \qquad x_{N-1} - x_{N-2} = \Delta s_{\text{last}}$$

The remaining $N-2$ interior points are distributed by the chosen interior law on $[x_0 + \Delta s_{\text{first}},\, x_{N-1} - \Delta s_{\text{last}}]$.

Use when: a node row has to enter and leave the interval at specific physical step sizes — typically because a multi-block seam fixes the step on the partner block. This is exactly what the cubed-sphere topology does internally at its cap–equator C1 seam (each meridian's first and last interval are pinned to the cap-edge first cell). The `pinned` law exposes the same construction as user-facing public API, available wherever `Spacing1D` is accepted.

```python
from shore.mesh.spacing import Spacing1D

sp = Spacing1D(
    law="pinned",
    first_step=0.05,   # required
    last_step=0.10,    # required
    interior="tanh",   # default "uniform"; controls the middle N-2 cells
    beta=4.0,          # tanh / tanh2 strength
)
sp.build(8, start=0.0, stop=1.0)   # → length-8 array; first and last
                                   # intervals exactly 0.05 and 0.10
```

Validation errors (`ParameterError`):
- `first_step` or `last_step` not strictly positive
- `first_step + last_step > stop - start` (the two pins overlap)
- `interior` is not one of `uniform`/`tanh`/`tanh2`
- `n < 2`

## Where spacing applies

SHORE has three independent spacing axes, each with its own set of allowed laws:

| Direction | CLI flag | Allowed laws | Default |
|---|---|---|---|
| Latitude (i) | `--surface-i-spacing` | `uniform`, `tanh`, `tanh2` | `uniform` |
| Longitude (j) | `--surface-j-spacing` | `uniform`, `tanh2` | `uniform` |
| Wall-normal (k) | `--volume-k-spacing` | `geometric`, `tanh`, `tanh2` | `geometric` |

### Latitude (i-direction)

Controls the meridional distribution of equator nodes between the south
cap seam and the north cap seam. Each j-column of the equator is a
great-circle arc parametrised in colatitude on
$[\theta_{\text{cap}},\, \pi - \theta_{\text{cap}}]$.

**"Start" and "end" in physical space:**

- `i = 0` = the **south cap seam** (colatitude $\theta_{\text{cap}}$, on the
  cap_south boundary).
- `i = ni - 1` = the **north cap seam** (colatitude
  $\pi - \theta_{\text{cap}}$, on the cap_north boundary).
- `i = ni // 2` = the **equator** ($\theta = \pi/2$).

**What clusters where:**

| Law | Dense region | Coarse region | Typical use |
|---|---|---|---|
| `tanh` | Near both cap seams | Equator (mid-meridian) | Boundary-layer-style refinement at both cap seams |
| `tanh2` | Near both cap seams | Equator | Symmetric variant of `tanh` |
| `uniform` | Everywhere | — | Default |

::: info Per-meridian seam pinning
The equator's first and last meridional intervals are **pinned** to the
cap-edge first cell at the corresponding seam vertex (per-j), giving C1
spacing continuity at every seam node when the cap and equator radial
step sizes are within ~2x of each other. The chosen `i_spacing` law
shapes the interior `ni − 2` cells between the two pinned endpoints. See
the [cubed-sphere topology page](/topologies/cubed-sphere) for the compatibility
window and fallback behaviour.
:::

### Longitude (j-direction)

Controls the **cap-edge tangent-plane distribution** along the
$[-\tan(\theta_{\text{cap}}),\, +\tan(\theta_{\text{cap}})]$ interval used
to build the gnomonic flat-square cap seed face. The cap edges are then
written bit-exactly into the equator's pole rows (DIRICHLET binding), so
the equator inherits the cap's j-distribution at the seam. Only `uniform`
and `tanh2` are available — `tanh` (one-sided) is not meaningful for a
4-corner symmetric distribution.

| Law | Effect | Use |
|---|---|---|
| `uniform` | Equispaced cap-edge nodes | Default; produces square-cell caps. |
| `tanh2` | Cluster at cap corners (sub-block junctions) | Refine the seam near the 4-block junctions. |

### Wall-normal (k-direction)

Controls the physical step sizes during hyperbolic marching. The interval spans $[0, L]$ where $L$ is the total wall-normal thickness (`--volume-k-thickness`).

**"Start" and "end" in physical space:**

- `k = 0` (start) = the **body wall** (the projected STL surface)
- `k = nk-1` (end) = the **far-field boundary** (outer limit of the volume)

**What clusters where:**

| Law | Dense region | Coarse region | Requires |
|---|---|---|---|
| `geometric` | Wall (k=0) | Far-field | `--ds`, `--growth` |
| `tanh` | Wall (k=0) | Far-field | `--volume-k-thickness` |
| `tanh2` | Wall (k=0) *and* far-field (k=nk-1) | Mid-layer region | `--volume-k-thickness` |

::: warning Geometric vs. tanh at the wall
Both `geometric` and `tanh` cluster near the wall, but they differ in how the outer boundary is determined:
- **Geometric**: no fixed thickness. The outer boundary position depends on `ds`, `growth`, and `nk`. More layers = farther out.
- **Tanh**: fixed total thickness `L`. The outer boundary is always at distance `L` from the wall. The first-layer thickness is determined by `nk`, `L`, and `beta` — `--ds` and `--growth` are ignored.
:::

Use `tanh` or `tanh2` when you need to control the far-field distance (e.g., matching an outer domain boundary). Use `geometric` when you care primarily about the first-layer thickness (e.g., $y^+$ resolution for turbulence modelling).

## Per-edge spacing

The four laws above are normally chosen *per axis*: one law for i, one for j, one for k, applied uniformly to the whole block. SHORE also lets each of a block's twelve edges carry its own `Spacing1D`. When the four parallel edges of an axis disagree, construction switches to **transfinite interpolation (TFI)** — the axis coordinate becomes a bilinear blend of the four edge distributions over the transverse face.

### Mechanism

Each `HexBlock.edge(name).spacing` field stores a `Spacing1D`. The constructor stamps the same instance onto all four parallel edges of an axis (the per-axis default), so an edge whose `spacing` differs from `block._default_spacing_*` has been overridden by the caller. Consumers — primitive constructors and `Mesh.march()` — read the per-edge values when they detect an override, and fall back to the per-axis default otherwise.

The `Spacing1D` instance can also carry `first_step` / `last_step` (used directly by the `pinned` law, ignored by other laws but readable as metadata).

### TFI blend

For axis i with four parallel edges keyed by their (j, k) endpoint corners, the parametric coordinate at $(i, j, k)$ is

$$t(i, j, k) = (1-\eta)(1-\zeta)\, t_{00}(i) + \eta(1-\zeta)\, t_{10}(i) + (1-\eta)\zeta\, t_{01}(i) + \eta\zeta\, t_{11}(i)$$

with $\eta = j / (n_j - 1)$, $\zeta = k / (n_k - 1)$. Each $t_{\alpha\beta}(i)$ is the corresponding edge's distribution, normalised to $[0, 1]$. The world coordinate is $\text{start} + (\text{stop} - \text{start})\, t(i, j, k)$. When all four edges carry the same `Spacing1D`, $t_{00} = t_{10} = t_{01} = t_{11}$ and the TFI collapses back to the legacy 1D distribution broadcast across the transverse face.

### Where it's wired

Per-edge support varies by topology and primitive — every case that *isn't* fully wired is so because of a structural constraint, not a missing implementation:

| Geometry | i / j | k |
|---|---|---|
| `BoxMesh` | full TFI via `edge_pins` | full TFI via `edge_pins` |
| `AnnularCylinderMesh` | full TFI on `r`, on `θ` for partial sectors. Periodic full annulus rejects per-edge `j`. | full TFI |
| `FlatCapsCylinder` | rejected — radial (i) and tangential (c) are pinned by the geometric construction | rejected — k-coupled SHARED seams across the 5 blocks must remain coincident layer by layer |
| `OGridTopology` | rejected — j-periodic, the four parallel i-edges are at adjacent azimuths under the wrap, not at distinct corners | per-edge k via `Mesh.march()` resolver |
| `CubedSphereTopology` | deferred — `build_cap_k0` / `build_equator` are not edge-aware (the C1 cap-equator pinning *is* per-edge, but it's geometric, not user-driven) | per-edge k via `Mesh.march()` resolver, with the constraint that all 6 blocks share one schedule |

### Per-edge k via `Mesh.march()`

For topologies that go through hyperbolic marching (cubed-sphere, O-grid), per-edge k-spacing is consumed by the marcher rather than at construction time. The user replaces the `spacing` field of every k-edge of every block with the same `Spacing1D` instance, and `Mesh.march()` reads it instead of its `(ds, growth)` arguments:

```python
from shore.mesh.spacing import Spacing1D

override = Spacing1D(law="geometric", ds=5e-4, growth=1.2)
for blk in mesh.blocks:
    for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
        blk.edge(name).spacing = override

mesh.march(ds=99, growth=99)   # ds, growth ignored — override wins
```

The resolver (`shore.mesh.mesh._resolve_k_steps`) detects the override by *identity* (each block's k-edges no longer point at `block._default_spacing_k`). Cross-block disagreement triggers a `UserWarning` and block 0's schedule wins — the constraint is geometric: all blocks march with one schedule because k-coupled SHARED / DIRICHLET seams must remain coincident layer by layer.

### Per-edge i/j on primitives

`BoxMesh` and `AnnularCylinderMesh` accept an `edge_pins` kwarg at construction:

```python
from shore.primitive import BoxMesh
from shore.mesh.spacing import Spacing1D

box = BoxMesh.build(
    vmin=(0, 0, 0), vmax=(10, 5, 3),
    ni=40, nj=30, nk=20,
    edge_pins={
        "j_lo_k_lo": Spacing1D(law="tanh", beta=4.0),  # bottom-front i-edge
        # other i-edges fall back to the global i_law (default uniform)
    },
)
```

Edge names are the canonical 12 from `shore.mesh.edge`. The user pins as many or as few edges as needed; unpinned edges inherit the per-axis `*_law` default. After construction, `block.edge(name).spacing` returns the actual `Spacing1D` used (either the user override or the axis default).

See the worked examples:
- `examples/03-per-edge-spacing/per_edge_ij_box.py` — per-edge i on a Box, with side-by-side comparison to the uniform tensor-product baseline
- `examples/03-per-edge-spacing/per_edge_k_spacing_cubed_sphere.py` — per-edge k on the cubed sphere
- `examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py` — per-edge k on the O-grid

## Choosing $\beta$

The clustering strength $\beta$ trades resolution in the clustered region against degradation in the coarse region:

| $\beta$ | Effect |
|---|---|
| 1–2 | Gentle — ratio of finest to coarsest cell ~2:1 |
| 3 (default) | Moderate — good balance for most geometries |
| 4–5 | Aggressive — tight clustering, but cells in the coarse region become very large |
| > 5 | Usually too extreme — cells degenerate in the coarse region |

A practical check: after generating the grid, run `shore check` and inspect the per-layer Jacobians. If `jmin` drops sharply in the transition zone between fine and coarse cells, reduce $\beta$ or switch to a milder law.

## Smoothing

Laplacian smoothing is applied after each hyperbolic marching step to suppress high-frequency wiggles introduced by the normal-based extrusion. It is controlled independently of spacing:

- `--smooth` (blend strength, 0–1): how far each node moves toward the average of its four neighbours. 0 disables smoothing entirely.
- `--smooth-iters` (integer): number of Laplacian sweeps per layer. More sweeps give a smoother result but increase run time.

Boundary nodes (non-periodic directions and the k=0 wall) are always clamped — they never move during smoothing. This preserves conformity between adjacent blocks in multi-block topologies.

## Practical examples

### High-res boundary layer with controlled far-field

A common CFD requirement: $y^+ \approx 1$ at the wall with a fixed outer boundary at 5 body lengths.

```bash
shore mesh fuselage.stl \
  --ni 80 --nj 120 --nk 60 \
  --volume-k-spacing tanh --volume-k-beta 4.0 --volume-k-thickness 5.0 \
  -o fuselage_bl.geo
```

### Cap-corner refinement

Cluster cells near the cap corners (sub-block junctions on the seam) using
tanh2 j-spacing. The equator's seam pole rows inherit the same
distribution by DIRICHLET binding, so the j-clustering is consistent
between cap and equator.

```bash
shore mesh body.stl \
  --ni 40 --nj 60 --nk 30 \
  --surface-j-spacing tanh2 --surface-j-beta 4.0 \
  --theta-cap-deg 30 \
  -o body.geo
```

### Cubed-sphere topology with meridional refinement

Concentrate equator cells near both cap seams using tanh i-spacing,
combined with tanh wall-normal clustering.

```bash
shore mesh body.stl \
  --ni 60 --nj 80 --nk 30 \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  --volume-k-spacing tanh --volume-k-beta 3.0 --volume-k-thickness 2.0 \
  --theta-cap-deg 30 \
  -o body.geo
```

## Python API

```python
from shore.volume.spacing import build_i_coords, build_k_steps, tanh_two_sided

# Latitude (equator meridian): tanh clustering near both cap seams
i_coords = build_i_coords(60, theta_start=0.524, theta_stop=2.618,
                           spacing="tanh", beta=4.0)

# Cap-edge tangent-plane (j-direction): tanh2 clustering at corners
half_side = 0.577  # tan(30 deg)
j_coords = tanh_two_sided(16, -half_side, half_side, beta=3.0)

# Wall-normal: tanh within fixed thickness
k_steps = build_k_steps(40, ds=0.01, spacing="tanh", beta=3.0, total_thickness=3.0)

# Wall-normal: tanh2 — fine at wall AND far-field
k_steps = build_k_steps(40, ds=0.01, spacing="tanh2", beta=4.0, total_thickness=3.0)
```

## See also

- [`shore mesh`](/reference/mesh) — full pipeline with all spacing flags
- [k=0 mesh construction](k0-construction) — how `j_spacing` and
  `i_spacing` interact at the cap-equator seam
- [Hyperbolic marching](hyperbolic-marching) — where `volume_k_spacing`
  controls the wall-normal step sizes
- [Cubed-sphere topology](/topologies/cubed-sphere) — full pipeline view
