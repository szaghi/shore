# Pipeline overview

SHORE's pipeline transforms a raw triangulated surface into a structured volume block in two stages. Each stage is independent and can be run, inspected, and re-run separately.

## Stages

```
┌─────────────────────────────────────────────────────────────────┐
│  Input: body.stl  (unstructured triangles, any complexity)       │
└──────────────────────────────┬──────────────────────────────────┘
                               │
                 ┌─────────────▼─────────────┐
                 │   Stage 1 — Projection     │  shore project
                 │                           │
                 │  Build a (ni×nj) lat-lon  │
                 │  sphere around centroid.  │
                 │  Cast ni×nj rays inward.  │
                 │  Intersect with STL BVH.  │
                 └─────────────┬─────────────┘
                               │  k=0 layer: (ni, nj, 3)
                 ┌─────────────▼─────────────┐
                 │   Stage 2 — Marching       │  shore march
                 │                           │
                 │  For k = 0 … nk-2:        │
                 │    compute normals(k)      │
                 │    advance k+1            │
                 │    smooth k+1             │
                 │    check Jacobian(k, k+1) │
                 └─────────────┬─────────────┘
                               │  volume: (ni, nj, nk, 3)
                 ┌─────────────▼─────────────┐
                 │   Output: body.geo         │  write_geo
                 │   Xall_GPU format          │
                 └───────────────────────────┘
                               │
              ┌────────────────┴────────────────┐
              │                                 │
┌─────────────▼─────────────┐   ┌──────────────▼──────────────┐
│  ParaView export           │   │  Interactive viewer          │
│  shore export body.geo     │   │  shore view body.geo         │
│  → body.vts (.vts format)  │   │  → PyVista window / Jupyter  │
└───────────────────────────┘   └─────────────────────────────┘
```

## Data shapes

| Quantity | Shape | Description |
|---|---|---|
| k=0 surface layer | `(ni, nj, 3)` | Structured surface on the body (x, y, z per node) |
| Volume grid | `(ni, nj, nk, 3)` | Full O-grid block |
| Geo file on disk | — | ASCII: header `ni nj nk`, then all x, all y, all z in Fortran order |

The j-index is periodic (O-grid topology): `j=nj-1` and `j=0` are adjacent neighbours.

## Intermediate format

Stage 1 saves its output as a compressed NumPy file (`.npz`):

```python
np.savez_compressed("body_k0.npz", k0=surface_layer)
# Load:
k0 = np.load("body_k0.npz")["k0"]   # shape (ni, nj, 3)
```

This allows inspecting the surface projection before committing to a full march, and re-running stage 2 with different parameters without repeating stage 1.

## Visualization

Two independent visualization paths are available after the volume is generated.

### ParaView (VTK export)

```bash
shore export body.geo        # writes body.vts
paraview body.vts
```

`shore export` converts the O-grid to a VTK StructuredGrid and attaches a `Jacobian` cell-data array. The j-direction is closed (`nj → nj+1`) because VTK has no concept of periodic connectivity.

### PyVista (interactive / Jupyter)

```bash
shore view body.geo          # desktop window
```

From Python or Jupyter:

```python
from shore.vis import view_geo
from shore.io.geo import read_geo

grid = read_geo("body.geo")
view_geo(grid, jupyter=True).show()
```

The wall layer is highlighted separately and the outer layers are coloured by Jacobian magnitude.

## Limitations

### Star-shaped constraint

Stage 1 requires that every ray from the body centroid hits the STL surface exactly once. This is called the **star-shaped** (or star-convex) property. It is satisfied by:

- Spheres, ellipsoids, cylinders
- Fuselages, nacelles, bluff bodies
- Wings and airfoils (approximately, if the centroid is inside the section)

It fails for:

- Bodies with re-entrant cavities (e.g. engine inlets seen from outside)
- Multi-component geometries (wing + fuselage as a single mesh)
- Surfaces with holes (non-watertight STL)

When the constraint fails, `project_sphere_to_surface` raises `ValueError` with a count of the rays that missed or double-hit. The fix is either to repair the STL (fill holes, split components) or to shift the projection centroid manually.

### No smooth far field

The outer boundary of the volume (k=nk-1) currently inherits the shape of the sphere projection template — it is not morphed to a box or matched to the Xall_GPU background grid. Box morphing is the next planned feature.

### Single block

SHORE currently generates one `(ni, nj, nk)` block per invocation. Multi-block decomposition for complex surfaces is on the roadmap.
