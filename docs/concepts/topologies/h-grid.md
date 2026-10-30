# H-grid topology

An H-grid fills the inside of a duct: i runs streamwise from the inlet
face to the outlet face, j and k span the cross-section, and the
cross-section faces are the duct walls. SHORE's H-grids are analytic
primitives without an internal body. A body in a channel is handled the
overset way: an H-grid for the channel plus a body-fitted grid around the
body ([Chimera mesh generation](/concepts/chimera)).

<Figure name="h-grid" caption="H-grid duct primitive: a single block with the walls on the cross-section faces." />

## Variants

`shore.primitive` exposes two analytic H-grid primitives:

| Variant | CLI | Cross-section | Block layout | Use case |
|---|---|---|---|---|
| **Rectangular** | `shore primitive h-grid-rect` | rectangle | 1 hex block, no singular axis | prismatic rectangular ducts; channels |
| **Circular** | `shore primitive h-grid-tube` | circle | 1 hex block, j-periodic, polar singular axis at the centerline | straight cylindrical ducts; simple SRM chambers |

Ducts of varying or non-circular cross-section are not covered by a
primitive; build them as a [blocking](/concepts/blocking) on the duct
geometry.

## Conventions (both variants)

- Block axes: ``i`` is **streamwise**, ``j`` and ``k`` are the two
  cross-section directions.
- World alignment: ``streamwise_axis`` (``"x"``, ``"y"``, or ``"z"``)
  selects which world axis the duct length runs along.  Default ``x``.
- Default face BCs: streamwise faces (``i_lo``, ``i_hi``) are
  ``FREE`` (inflow / outflow); cross-section faces are ``WALL``.
- The ``j_lo`` / ``j_hi`` faces are SHARED by ``periodic_j = True``
  on the circular tube (closed circumferential direction); on the
  rectangular duct they are independent ``WALL`` faces.

## Rectangular duct

Single hex block of shape ``(ni, nj, nk, 3)``.  ``length`` / ``width``
/ ``height`` set the three Cartesian extents (``length`` along the
streamwise axis, ``width`` along block ``j``, ``height`` along block
``k``).  Tanh / tanh2 spacing on the cross-section axes clusters
nodes near the duct walls — the natural choice for boundary-layer
resolution.

### Face boundary conditions

| Face | BC | Meaning |
|---|---|---|
| `i_lo` | `FREE` | Streamwise inflow |
| `i_hi` | `FREE` | Streamwise outflow |
| `j_lo` / `j_hi` | `WALL` | Two of the four duct walls |
| `k_lo` / `k_hi` | `WALL` | The other two duct walls |

In `cc.par` the streamwise `FREE` faces take the block's
`free_face_bc`; set a positive natural code for them, per block or per face
(3 / 4 for inflow / outflow).

## Circular tube

Single hex block of shape ``(n_axial, n_circ, n_radial, 3)`` with a
**polar 2D cross-section**:

- ``j`` = circumferential angle (``0`` → ``2π``), **periodic** by
  construction (`block.periodic_j = True`).
- ``k`` = radial coordinate from the centerline (``k=0``) to the
  outer wall (``k=n_radial-1``).
- ``i`` = streamwise.

### Polar singular axis

The ``k = 0`` row at the duct centerline is a **polar singular axis**:
all ``n_circ`` angular nodes collapse to a single point at the
centerline of the duct.  Formally, the cells at ``k=0`` are degenerate
(zero volume on the singular axis). The face is tagged `WALL`, as the
code sets it; a solver without singular-axis handling needs the
`h-grid-rect` instead.

### Face boundary conditions

| Face | BC | Meaning |
|---|---|---|
| `i_lo` | `FREE` | Streamwise inflow |
| `i_hi` | `FREE` | Streamwise outflow |
| `j_lo` / `j_hi` | SHARED (periodic) | Circumferential periodicity |
| `k_lo` | `WALL` | Centerline singular axis (degenerate) |
| `k_hi` | `WALL` | Outer wall (real wall) |

## Parameters

| Command | Sizes | Nodes | Spacing |
|---|---|---|---|
| `h-grid-rect` | `--length` (1.0), `--width` (0.1), `--height` (0.1) | `--ni` 40, `--nj` 20, `--nk` 20 | `--i/j/k-spacing` `uniform` / `tanh` / `tanh2`, `--i/j/k-beta` 3.0 |
| `h-grid-tube` | `--length` (1.0), `--radius` (0.05) | `--n-axial` 40, `--n-circ` 32, `--n-radial` 8 | `--i-spacing`, `--k-spacing` (radial, centreline to wall), betas 3.0 |

Both take `--streamwise-axis` (`x` default). Two-sided `tanh2` on the
cross-section axes clusters toward the walls. In Python they are
`HGridRectMesh.build(...)` and `HGridTubeMesh.build(...)`, with the same
`PrimitiveMesh` contract as the other [primitives](primitives). Commands:
[`shore primitive`](/reference/primitive).

## See also

- [Topology catalogue](./) — pick the right topology for your case.
- [`shore primitive`](/reference/primitive) — full CLI reference for all primitives.
- [Chimera assembly tutorial](/tutorials/chimera-assembly) — a body-fitted grid inside a background.
- [Per-face natural BC codes](/concepts/load-balancing#per-face-natural-bc-codes) — different codes on the inlet and outlet faces.
