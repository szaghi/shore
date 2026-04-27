# shore.topology.sphere_cap

Marching primitives shared by the [cubed-sphere topology](/topologies/cubed-sphere):
single-layer hyperbolic step + Laplacian smooth, cap boundary pinning,
and Jacobian guard. The high-level orchestration lives in
`shore.mesh.Mesh._march_caps`; the cap k=0 seed and the equator k=0
construction live in `shore.surface.cap` and `shore.surface.equator`
respectively.

## Public API

There is no public function in this module — application code uses the
OOP entry point [`shore.mesh.Mesh.from_stl`](mesh) and `Mesh.march()`.
The functions below are documented for users who want to subclass the
marcher or call the primitives directly.

## Internal helpers

### `_adjust_nj`

```python
def _adjust_nj(nj: int) -> int
```

Round `nj` down to the nearest multiple of 4 so the equatorial band can
be split into 4 equal sectors. Emits `warnings.warn` if adjusted; returns
`nj` unchanged if already divisible by 4.

---

### `_overwrite_cap_boundaries`

```python
def _overwrite_cap_boundaries(
    cap_layer: np.ndarray,
    eq_layers: list[np.ndarray],
    pole: str = "north",
) -> None
```

Overwrite the four boundary edges of a cap layer with the equatorial
sub-block pole rows. The mapping differs between the north and south
caps due to opposite face windings:

**North cap (z<0 of the body):**

| Cap edge | Equatorial source |
|---|---|
| `cap[0, :]` | `sub0[ni-1, :]` |
| `cap[:, -1]` | `sub1[ni-1, :]` |
| `cap[-1, :]` | `sub2[ni-1, ::-1]` |
| `cap[:, 0]` | `sub3[ni-1, ::-1]` |

**South cap (z>0 of the body):**

| Cap edge | Equatorial source |
|---|---|
| `cap[:, 0]` | `sub0[0, :]` |
| `cap[-1, :]` | `sub1[0, :]` |
| `cap[:, -1]` | `sub2[0, ::-1]` |
| `cap[0, :]` | `sub3[0, ::-1]` |

---

### `_march_one_layer`

```python
def _march_one_layer(
    layer: np.ndarray,
    step: float,
    smooth_strength: float,
    smooth_iters: int,
    periodic: tuple[bool, bool],
    normals_override: np.ndarray | None = None,
    ghost_rows: dict[str, np.ndarray] | None = None,
) -> np.ndarray
```

Advance a single block one wall-normal layer:

1. Compute outward normals via `_compute_normals(layer, periodic,
   ghost_rows)`. When `ghost_rows` is supplied (one entry per face from
   `{'i_lo', 'i_hi', 'j_lo', 'j_hi'}`), the corresponding boundary tangent
   uses a centred difference against the ghost row instead of a one-sided
   fallback. This is what makes block-seam spacing continuous during
   marching (see [seam-aware normals](/topologies/cubed-sphere#hyperbolic-marching-with-seam-aware-normals)).
2. Add `step * normals` to the layer.
3. Apply `smooth_iters` Laplacian sweeps with blend factor
   `smooth_strength`; non-periodic boundary rows/columns are clamped.

Used by `Mesh._march_caps` for each of the 4 equatorial sub-blocks.

---

### `_march_cap_layer`

```python
def _march_cap_layer(
    layer: np.ndarray,
    step: float,
    smooth_strength: float,
    smooth_iters: int,
    eq_rows: list[np.ndarray],
    pole: str = "north",
    normals_override: np.ndarray | None = None,
    ghost_rows: dict[str, np.ndarray] | None = None,
) -> np.ndarray
```

Advance one cap layer:

1. Compute outward normals (with optional `ghost_rows` for seam-aware
   centred differences against the partner sub-blocks).
2. Add `step * normals` to the layer.
3. Apply `smooth_iters` Laplacian sweeps; all 4 boundary rows clamped.
4. **Dirichlet pin** — overwrite the 4 boundary edges with `eq_rows`
   (the next-layer pole rows from the 4 equatorial sub-blocks) using
   `_overwrite_cap_boundaries`.

The Dirichlet pin in step 4 is what enforces seam C0 conformity at every
layer; the seam-aware normals in step 1 enforce continuous spacing
across the seam.

---

### `_check_jacobian`

```python
def _check_jacobian(
    layer: np.ndarray,
    new_layer: np.ndarray,
    periodic_j: bool,
    label: str,
    k: int,
    exclude_corners: bool = False,
    tol: float = 0.0,
) -> float
```

Compute the minimum cell Jacobian between two consecutive layers and
raise `MarchingError` if it is below `tol`. When `exclude_corners=True`,
the 4 corner cells of the layer are masked out — these cells border the
4-sub-block junction and have near-zero signed volume by construction.

`Mesh._march_caps` calls this after each layer with `tol=-1e-6` to
absorb floating-point noise around the corner junctions.

::: tip Internal functions
Functions prefixed with `_` are internal implementation details. Use
[`shore.mesh.Mesh.from_stl`](mesh) for all application code.
:::

## Related modules

- `shore.surface.cap` — gnomonic flat-square cap k=0 seed
  (`build_cap_k0`).
- `shore.surface.equator` — great-circle meridian equator k=0 with
  per-j seam pinning (`build_equator`, `_constrained_distribution`).
- `shore.mesh.topology.CubedSphereTopology` — builds and wires the 6 blocks.
- `shore.mesh.Mesh._march_caps` — orchestrates the per-layer march
  with seam-aware ghost rows for all 6 blocks.

## VTK export

The 6 blocks can be exported to a ParaView MultiBlock file via
`shore.io.vtk.write_caps_vtm`:

```python
from shore.io.vtk import write_caps_vtm

write_caps_vtm("sphere.vtm", [b.nodes for b in mesh.blocks], with_jacobian=True)
```

This writes 6 `.vts` files plus one `.vtm` MultiBlock that groups them
under names `sub0`, `sub1`, `sub2`, `sub3`, `cap_north`, `cap_south`.

## See also

- [Cubed-sphere topology](/topologies/cubed-sphere) — conceptual overview
- [`shore.volume.hyperbolic`](api-hyperbolic) — `_compute_normals`,
  `_laplacian_smooth`, `_stl_normals` (with `ghost_rows` support)
- [`shore.volume.jacobian`](api-jacobian) — `jacobian_field`,
  `_jacobian_field_open`, `min_jacobian_two_layers`
- [`shore.io.vtk`](api-vtk) — VTK export including `write_caps_vtm`
