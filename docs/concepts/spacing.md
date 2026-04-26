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
the [caps topology page](/concepts/sphere-cap) for the compatibility
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

### Caps topology with meridional refinement

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
- [Caps topology](/concepts/sphere-cap) — how cap and equator spacing
  interact at the seam
- [Hyperbolic marching](/concepts/hyperbolic-marching) — the extrusion algorithm
