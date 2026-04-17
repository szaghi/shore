# Spacing laws

SHORE provides three 1D spacing distributions for controlling node density along each structured direction. These are implemented in [`shore.volume.spacing`](/reference/api-spacing) using pure numpy — no scipy dependency.

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

Controls the angular distribution of ray directions during surface projection. The interval spans $[\theta_{\text{cap}},\, \pi - \theta_{\text{cap}}]$ in colatitude.

**"Start" and "end" in physical space:**

- `i = 0` (start of interval) = the **north cap boundary** (colatitude $\theta_{\text{cap}}$, near the north pole)
- `i = ni-1` (end of interval) = the **south cap boundary** (colatitude $\pi - \theta_{\text{cap}}$, near the south pole)
- `i = ni//2` (midpoint) = the **equator** ($\theta = \pi/2$)

**What clusters where:**

| Law | Dense region | Coarse region | Typical use |
|---|---|---|---|
| `tanh` | Near *both* cap boundaries (i=0 and i=ni-1) | Equatorial band | Caps topology: concentrate cells where equatorial sub-blocks meet the polar caps |
| `tanh2` (default, no target) | Near *both* cap boundaries | Equatorial band | Same as `tanh` but symmetric — no bias toward north or south |
| `uniform` | Everywhere | — | Default for smooth bodies with no polar singularity concerns |

::: tip Why does `tanh` cluster at *both* ends?
The i-direction interval spans a colatitude range symmetric about $\pi/2$. Since `tanh` clusters at the *start* ($\theta_{\text{cap}}$), the mirror-symmetric end ($\pi - \theta_{\text{cap}}$) receives the same clustering by the geometry of the sphere. Use `tanh` when you want finer cells at the cap boundaries and don't care about the asymmetry direction.
:::

### Longitude (j-direction)

Controls the azimuthal distribution of ray directions. The interval spans $[0, 2\pi)$ (periodic). One-sided tanh is not meaningful for a periodic axis, so only `uniform` and `tanh2` are available.

**What `--j-cluster-lon` controls:**

The `--j-cluster-lon` flag specifies the longitude $\phi_0$ (in radians) where the dense region is centred. Cells are finest in a band around $\phi = \phi_0$ and coarsest at the antipodal longitude $\phi = \phi_0 + \pi$.

| `--j-cluster-lon` value | Physical location | Use case |
|---|---|---|
| `0.0` | $\phi = 0$ (positive x-axis) | Leading edge of an airfoil aligned along x |
| `pi / 2` ($\approx 1.571$) | $\phi = \pi/2$ (positive y-axis) | Feature at 90 deg from the reference axis |
| `pi` ($\approx 3.141$) | $\phi = \pi$ (negative x-axis) | Trailing edge |
| *(omit)* | $\phi = \pi$ (default) | Symmetric clustering at $\phi = 0$ and $\phi = \pi$ |

```bash
# Cluster at the leading edge (phi=0)
--surface-j-spacing tanh2 --j-cluster-lon 0.0

# Cluster at a 45-deg feature (phi=pi/4)
--surface-j-spacing tanh2 --j-cluster-lon 0.7854

# Default: symmetric at phi=pi (dense at 0 and pi, coarse at +-pi/2)
--surface-j-spacing tanh2
```

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

### Leading-edge refinement on an airfoil

Cluster longitude cells around $\phi = 0$ where the leading edge is, combined with tanh wall-normal spacing.

```bash
shore mesh airfoil.stl \
  --ni 40 --nj 80 --nk 40 \
  --surface-j-spacing tanh2 --surface-j-beta 5.0 --j-cluster-lon 0.0 \
  --volume-k-spacing tanh --volume-k-beta 3.0 --volume-k-thickness 3.0 \
  -o airfoil_le.geo
```

### Caps topology with boundary-layer refinement

Concentrate latitude cells at the cap boundaries where the equatorial and polar blocks meet.

```bash
shore mesh body.stl \
  --topology caps --ni 60 --nj 80 --nk 30 \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  --cap-pull 0.4 \
  -o body_caps.geo
```

### Experimenting with spacing via the two-stage pipeline

Project once, then try different volume distributions without re-projecting:

```bash
# Stage 1: project once
shore project wing.stl --ni 60 --nj 90 \
  --surface-j-spacing tanh2 --j-cluster-lon 0.0 \
  -o wing_k0.npz

# Stage 2a: geometric spacing (default)
shore march wing_k0.npz --nk 40 --ds 1e-3 --growth 1.10 -o wing_geom.geo

# Stage 2b: tanh spacing within thickness=2
shore march wing_k0.npz --nk 40 --volume-k-spacing tanh --volume-k-beta 3.0 --volume-k-thickness 2.0 -o wing_tanh.geo

# Stage 2c: tanh2 for bounded domain
shore march wing_k0.npz --nk 40 --volume-k-spacing tanh2 --volume-k-beta 4.0 --volume-k-thickness 2.0 -o wing_tanh2.geo

# Compare quality
shore check wing_geom.geo
shore check wing_tanh.geo
shore check wing_tanh2.geo
```

## Python API

```python
from shore.volume.spacing import build_i_coords, build_j_coords, build_k_steps

# Latitude: tanh clustering near cap boundaries
i_coords = build_i_coords(60, theta_start=0.1, theta_stop=3.04, spacing="tanh", beta=4.0)

# Longitude: tanh2 clustering at the leading edge (phi = 0)
j_coords = build_j_coords(80, spacing="tanh2", beta=5.0, cluster_lon=0.0)

# Wall-normal: tanh within fixed thickness
k_steps = build_k_steps(40, ds=0.01, spacing="tanh", beta=3.0, total_thickness=3.0)

# Wall-normal: tanh2 — fine at wall AND far-field
k_steps = build_k_steps(40, ds=0.01, spacing="tanh2", beta=4.0, total_thickness=3.0)
```

## See also

- [`shore mesh`](/reference/mesh) — full pipeline with all spacing flags
- [`shore project`](/reference/project) — surface spacing (stage 1)
- [`shore march`](/reference/march) — volume spacing (stage 2)
- [Hyperbolic marching](/concepts/hyperbolic-marching) — the extrusion algorithm
