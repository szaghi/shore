# shore.surface.cap & shore.surface.equator

The two k=0 builders consumed by [`CubedSphereTopology`](./api-topology):

- `shore.surface.cap.build_cap_k0` — gnomonic flat-square cap k=0
  seed (one cap face per pole).
- `shore.surface.equator.build_equator` — great-circle meridian
  equator k=0 from the two cap seam rings, with per-j C1 step
  matching at the cap-equator seam.

Both consume a [`Projector`](./api-projector) and emit a coordinate
array on the body surface. They are decoupled from the `HexBlock`
wiring — `CubedSphereTopology.build` calls them, then stitches the
output into the 6-block ring + cap allocations.

For the algorithm, see
[k=0 mesh construction](/algorithm/k0-construction).

## `build_cap_k0`

```python
def build_cap_k0(
    theta_cap: float,
    nc: int,
    pole: str,
    projector: object,
    j_spacing: str = "uniform",
    j_beta: float = 3.0,
    phi_offset: float = 0.0,
) -> np.ndarray
```

Build the `(nc, nc, 3)` cap k=0 seed face by gnomonic flat-square
projection: a flat diamond is placed in the tangent plane at the
cap pole with corners on the cardinal axes at distance
`tan(theta_cap)` (4 corners at azimuths `{0°, 90°, 180°, 270°}`,
aligned with sub-block junctions). Interior nodes by bilinear
interpolation, then radially projected through the cached BVH onto
the STL.

| Parameter | Description |
|---|---|
| `theta_cap` | Polar exclusion angle in radians. Cap corners sit at colatitude `theta_cap`. Cap-edge midpoints sit at colatitude `atan(tan(theta_cap)/sqrt(2))` — slightly poleward of the corners (the natural shape of a great-circle arc between two parallel-circle corners). |
| `nc` | Cap resolution. Must equal `nj // 4 + 1` so the cap edges have the same node count as the sub-block pole rows. |
| `pole` | `"north"` (cap_north, `z<0` of body) or `"south"` (cap_south, `z>0`). The two caps differ in their `(z, y)` sign convention so the index winding matches `_overwrite_cap_boundaries`. |
| `projector` | A [`Projector`](./api-projector) with cached BVH and the body's anchor. |
| `j_spacing` | Tangent-plane node distribution along cap edges: `"uniform"` (default) or `"tanh2"`. Inherited by the equator seam. |
| `j_beta` | Clustering strength for `"tanh2"`. |
| `phi_offset` | Azimuthal rotation (radians) of the gnomonic diamond about the polar axis. Default `0.0` places cap corners at world azimuths `{0, pi/2, pi, 3*pi/2}` (the `CubedSphereTopology` convention). `OHTopology` passes `pi/4` so the four sub-block *centres* face the world cardinal axes — what makes the 1:1 sub-block-to-box-face mapping possible for the wrap-channel hybrid. |

**Returns** `numpy.ndarray (nc, nc, 3)` — cap k=0 layer with nodes
on the body surface.

The 4-fold symmetry of the gnomonic square gives **square cap
cells by construction** (no concentric pole pattern, AR ≈ 1.0); see
[`Cubed-sphere topology`](/topologies/cubed-sphere) and
[k=0 mesh construction — gnomonic flat-square cap](/algorithm/k0-construction).

## `build_equator`

```python
def build_equator(
    cap_south: np.ndarray,
    cap_north: np.ndarray,
    ni: int,
    theta_cap: float,
    projector: object,
    i_spacing: str = "uniform",
    i_beta: float = 3.0,
) -> np.ndarray
```

Build the `(ni, nj, 3)` equator k=0 layer from the two cap seam
rings.

- Row `i=0` is the south seam (cap_south boundary, traversed in
  sub-block order); row `i=ni-1` is the north seam. Both are
  bit-exact copies of the cap edges (no projector round-trip).
- Interior rows: per-j slerp meridian from south seam to north seam,
  ray-cast through `projector` onto the STL.
- Each meridian's first/last interval is **per-j pinned** to the
  cap-edge first cell at the corresponding seam vertex (C1 spacing
  match at every seam vertex — including the 3-edge cap-corner
  junctions). The pinning is enabled when the cap and equator radial
  step magnitudes are comparable; if the equator's natural first
  step disagrees with the cap by more than ~2×, the equator falls
  back to its unpinned `i_spacing` law and emits a `UserWarning`.

| Parameter | Description |
|---|---|
| `cap_south`, `cap_north` | Cap k=0 layers, `(nc, nc, 3)` each, as returned by `build_cap_k0`. |
| `ni` | Number of equator latitude rows including both seam rows. |
| `theta_cap` | Polar exclusion angle in radians (matches the value used to build the caps). |
| `projector` | The same [`Projector`](./api-projector) used for the cap build. |
| `i_spacing` | Latitude spacing law: `"uniform"` (default) / `"tanh"` / `"tanh2"`. |
| `i_beta` | Clustering strength for tanh / tanh2. |

**Returns** `numpy.ndarray (ni, nj, 3)` with `nj = 4 * (nc - 1)`.

## Example

```python
from shore.surface.io import load_stl
from shore.surface.anchor import resolve_anchor
from shore.surface.projector import Projector
from shore.surface.cap import build_cap_k0
from shore.surface.equator import build_equator

mesh = load_stl("body.stl")
proj = Projector(mesh, resolve_anchor(mesh))

theta_cap, nc, ni = 0.5236, 16, 40   # ~30 deg, 16 nodes per cap edge

cap_s = build_cap_k0(theta_cap, nc, pole="south", projector=proj)
cap_n = build_cap_k0(theta_cap, nc, pole="north", projector=proj)
equator = build_equator(cap_s, cap_n, ni=ni, theta_cap=theta_cap, projector=proj)

print(cap_s.shape, cap_n.shape, equator.shape)
# (16, 16, 3) (16, 16, 3) (40, 60, 3)
```

In practice this orchestration is hidden inside
`CubedSphereTopology.build` — call it via
[`Mesh.from_stl(... topology=CubedSphereTopology())`](./api-mesh).
The two functions are documented separately for users who want to
read the kernel or reuse the cap k=0 / equator k=0 algorithm in a
new topology.

## See also

- [`shore.mesh.topology.CubedSphereTopology`](./api-topology) — the
  consumer that wires both outputs into the 6-block layout.
- [`shore.surface.projector`](./api-projector) — `Projector`.
- [k=0 mesh construction](/algorithm/k0-construction) — gnomonic
  cap, slerp equator, C1 seam pinning.
- [Cubed-sphere topology](/topologies/cubed-sphere) — block layout.
