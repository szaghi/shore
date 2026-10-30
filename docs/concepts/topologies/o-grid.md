# O-grid topology

The O-grid is a single block wrapped around the body: closed (periodic)
in one surface direction, so the grid lines around a cross-section look
like the letter O. In SHORE it is a lat-lon grid ray-cast onto the body
and marched outward, built by `shore mesh --topology ogrid` (or
`[grid] topology = "ogrid"` in a config). It writes one `<stem>.geo`.

<Figure name="ogrid-single" caption="A single-block O-grid: the wall grid and one j-sheet of the volume." />

## Block

| Index | Direction | Range | Periodic |
|---|---|---|---|
| i | colatitude | $[\theta_{cap},\ \pi - \theta_{cap}]$ | no |
| j | longitude, `nj` nodes around the body | $[0, 2\pi)$ | yes: `j = nj − 1` is adjacent to `j = 0` |
| k | wall normal | wall to outer layer | no |

The directions come from `sphere_latlon_directions` and are cast from the
anchor in one batch ([Surface projection](/concepts/projection)). The
block's `k_lo` is `WALL`; normals, smoothing and the Jacobian check wrap
around j with modular indexing; the i boundaries use one-sided tangents.

::: warning The poles are not covered
The first and last i-rows are rings at colatitude `theta_cap` (default
9°) from each pole, not points. The two polar regions inside those rings
have no cells. The O-grid is therefore an open band around the body: in
an overset assembly the poles must be covered by another grid. Use the
[cubed sphere](cubed-sphere) for a closed wall grid.
:::

## Compared with the cubed sphere

| | O-grid | Cubed sphere |
|---|---|---|
| Blocks | 1 | 6 |
| Poles | uncovered (rings at `theta_cap`) | covered by two caps |
| j direction | periodic | four sub-blocks with `SHARED` seams |
| Output | `<stem>.geo` | six `<stem>_<label>.geo` |
| Use | a single periodic block for tools that expect one | wall grids of closed bodies |

## Parameters

All `shore mesh` surface and marching flags apply: `--ni`, `--nj` (the
block's j node count directly), `--nk`, `--theta-cap` /
`--theta-cap-deg`, the surface i/j and wall-normal spacing laws,
`--march-scheme`, `--smooth*`. `--surface-j-spacing tanh2` clusters
toward longitude $\pi$. Per-edge spacing: only the four k-edges can be
overridden (j is periodic, so i-edges are not at distinct corners);
`Mesh.march` reads the override
([Spacing laws](/concepts/spacing-laws#per-edge-spacing)).

## Python

```python
from shore.mesh import Mesh, OGridTopology, Spacing1D

m = Mesh.from_stl(
    "sphere.stl",
    topology=OGridTopology(),
    ni=40, nj=60, nk=30,
    theta_cap_deg=9.0,
    spacing_k=Spacing1D(law="geometric", ds=0.02, growth=1.1),
)
m.march(ds=0.02, growth=1.1)
m.write_geo("sphere")  # sphere_ogrid.geo (Mesh.write_geo appends the block label)
```

## See also

- [Cubed sphere](cubed-sphere): the closed six-block alternative
- [Surface projection](/concepts/projection), [Marching](/concepts/hyperbolic-marching)
- [`shore mesh`](/reference/mesh), [`OGridTopology`](/reference/api/shore.mesh.topology#OGridTopology)
