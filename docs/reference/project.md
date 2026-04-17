# shore project

Project a structured sphere grid onto an STL surface, producing the k=0 layer only.

Use this command to inspect the surface projection before committing to a full march, or when you want to cache stage 1 and experiment with different marching parameters without repeating the ray casting.

## Synopsis

```
shore project [OPTIONS] INPUT
```

## Arguments

| Argument | Description |
|---|---|
| `INPUT` | Path to the input `.stl` file |

## Options

### Grid resolution

| Option | Default | Description |
|---|---|---|
| `--ni INTEGER` | 40 | Latitude resolution |
| `--nj INTEGER` | 60 | Longitude resolution |

### Topology

| Option | Default | Description |
|---|---|---|
| `--theta-cap FLOAT` | — | Polar exclusion angle in radians |
| `--theta-cap-deg FLOAT` | — | Same as `--theta-cap` but in degrees |

### Surface spacing

Surface spacing redistributes the angular positions of the projection rays. Once baked into the `.npz`, it cannot be changed without re-projecting.

| Option | Default | Description |
|---|---|---|
| `--surface-i-spacing TEXT` | `uniform` | Latitude spacing: `'uniform'`, `'tanh'` (cluster at cap boundaries), `'tanh2'` (both ends) |
| `--surface-i-beta FLOAT` | 3.0 | Clustering strength for latitude tanh/tanh2 |
| `--surface-j-spacing TEXT` | `uniform` | Longitude spacing: `'uniform'`, `'tanh2'` (cluster at a target longitude) |
| `--surface-j-beta FLOAT` | 3.0 | Clustering strength for longitude tanh2 |
| `--j-cluster-lon FLOAT` | — | Target longitude $\phi_0$ in radians where cells are finest (tanh2 only). `0.0` = +x axis, `pi/2` = +y axis, `pi` = -x axis. Defaults to `pi` when omitted. |

### I/O and display

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | `<INPUT stem>_k0.npz` | Output `.npz` file |
| `-v` / `--verbose` | off | Show Rich progress |
| `--help` | — | Show this help and exit |

## Output format

The output is a compressed NumPy archive (`.npz`) with a single array `k0` of shape `(ni, nj, 3)`:

```python
import numpy as np
data = np.load("sphere_k0.npz")
k0 = data["k0"]          # shape (ni, nj, 3)
print(k0.shape)           # (40, 60, 3)
```

The array holds the `(x, y, z)` coordinates of each structured surface node. Row `i` spans the latitude direction, column `j` spans the longitude direction (periodic).

Surface spacing is baked into the projection at this stage — the angular distribution of ray directions is set here and cannot be changed later without re-running `shore project`.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input file not found or projection failed (rays missed body) |

## Examples

Project a sphere and inspect in Python:

```bash
shore project sphere.stl --ni 40 --nj 60 -o sphere_k0.npz -v
```

```python
import numpy as np
k0 = np.load("sphere_k0.npz")["k0"]
radii = np.linalg.norm(k0, axis=-1)
print(f"mean radius = {radii.mean():.4f}")   # ≈ 1.0
print(f"std  radius = {radii.std():.4f}")    # ≈ 0.0 for a perfect sphere
```

Verify the projected layer looks correct before marching:

```bash
shore project fuselage.stl --ni 60 --nj 90 -o fuselage_k0.npz
# Then inspect fuselage_k0.npz in pyvista or matplotlib before marching
shore march fuselage_k0.npz --nk 40 --ds 1e-4 -o fuselage.geo -v
```

Latitude clustering near the cap boundaries (finer cells at i=0 and i=ni-1, coarser at the equator):

```bash
shore project body.stl \
  --ni 60 --nj 80 \
  --surface-i-spacing tanh --surface-i-beta 4.0 \
  -o body_k0.npz
```

Longitude clustering toward a geometric feature. The `--j-cluster-lon` value is the azimuthal angle in radians of the feature you want to refine. Cells will be finest in a band around that longitude and coarsest at the antipodal point ($\phi_0 + \pi$).

```bash
# Leading edge at phi=0 (positive x-axis)
shore project airfoil.stl \
  --ni 40 --nj 80 \
  --surface-j-spacing tanh2 --surface-j-beta 5.0 --j-cluster-lon 0.0 \
  -o airfoil_k0.npz

# Feature at 45 degrees (phi = pi/4)
shore project nacelle.stl \
  --ni 40 --nj 80 \
  --surface-j-spacing tanh2 --surface-j-beta 4.0 --j-cluster-lon 0.7854 \
  -o nacelle_k0.npz
```

## Algorithm

See [Pipeline overview — Stage 1](/concepts/pipeline) for a full description of the sphere projection algorithm.

The projection uses trimesh's `RayMeshIntersector` with `multiple_hits=False` — only the first intersection along each ray (closest to the centroid) is recorded. Ray origins are all placed at the body centroid; directions are the unit vectors of a lat-lon sphere grid.

When surface spacing is non-uniform, the lat-lon directions are redistributed using the selected spacing law before ray casting. See [Spacing laws](/concepts/spacing) for details.

## See also

- [`shore march`](march) — extrude the saved k=0 layer
- [`shore mesh`](mesh) — run the full pipeline
- [Spacing laws](/concepts/spacing) — clustering distributions explained
- [Pipeline overview](/concepts/pipeline)
