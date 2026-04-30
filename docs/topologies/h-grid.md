# H-grid topology

An H-grid is a structured grid topology for **internal flow in a
duct**: the mesh fills the duct interior, with the streamwise direction
running through the inlet and outlet faces and the cross-section
faces meshing the duct walls.  No body wraps the inside — the H-grid
is empty inside, by design.

The name comes from the cross-section view of an H-grid wrapping a
body in a channel (the wraps for the body and the channel walls,
viewed end-on, look like the letter H).  This page covers the
**body-less H-grid baseline**: a duct without an internal body.  When
a body in a channel is needed, the standard CFD pattern is to compose
an H-grid (filling the channel) with an O-grid or C-grid (wrapping
the body), via SHORE's existing
[Chimera assembly](/guide/cc-par-pipeline) pipeline.

## Variants

`shore.primitive` exposes two analytic H-grid primitives:

| Variant | CLI | Cross-section | Block layout | Use case |
|---|---|---|---|---|
| **Rectangular** | `shore primitive h-grid-rect` | rectangle | 1 hex block, no singular axis | prismatic rectangular ducts; channels |
| **Circular** | `shore primitive h-grid-tube` | circle | 1 hex block, j-periodic, polar singular axis at the centerline | straight cylindrical ducts; simple SRM chambers |

A third variant — **STL-driven H-grid** — for variable-cross-section
or non-circular ducts is on the roadmap (`shore mesh --topology
hgrid`).  The user supplies an inner-surface STL of the duct walls;
the slicer extracts cross-section polylines per streamwise station
and fits a structured 2D mesh to each.  Star-shaped cross-sections
from the centroid are the planned scope.

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

```bash
shore primitive h-grid-rect \
    --length 1.0 --width 0.1 --height 0.2 \
    --ni 80 --nj 20 --nk 24 \
    --streamwise-axis x \
    --j-spacing tanh2 --j-beta 4.0 \
    --k-spacing tanh2 --k-beta 4.0 \
    -o duct
# → duct.geo
```

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

For a chimera assembly the streamwise FREE faces become inflow / outflow
BCs in `cc.par`.

## Circular tube

```bash
shore primitive h-grid-tube \
    --length 0.5 --radius 0.05 \
    --n-axial 60 --n-circ 32 --n-radial 10 \
    --streamwise-axis x \
    --k-spacing tanh --k-beta 3.0 \
    -o tube
# → tube.geo
```

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
(zero-volume on the singular axis); finite-volume CFD codes that
expect singular-axis handling (Xall's sphere-grid path is one) deal
with this transparently.  Codes that don't can be fed an
``h-grid-rect`` instead.

### Face boundary conditions

| Face | BC | Meaning |
|---|---|---|
| `i_lo` | `FREE` | Streamwise inflow |
| `i_hi` | `FREE` | Streamwise outflow |
| `j_lo` / `j_hi` | SHARED (periodic) | Circumferential periodicity |
| `k_lo` | `WALL` | Centerline singular axis (degenerate) |
| `k_hi` | `WALL` | Outer wall (real wall) |

## Python API

```python
from shore.primitive import HGridRectMesh, HGridTubeMesh

# Rectangular duct.
duct = HGridRectMesh.build(
    length=1.0, width=0.1, height=0.2,
    ni=80, nj=20, nk=24,
    streamwise_axis="x",
    j_law="tanh2", j_beta=4.0,
    k_law="tanh2", k_beta=4.0,
)
duct.write_geo("duct.geo")
duct.write_vtk("duct.vtm")

# Circular tube.
tube = HGridTubeMesh.build(
    length=0.5, radius=0.05,
    n_axial=60, n_circ=32, n_radial=10,
    streamwise_axis="x",
    k_law="tanh", k_beta=3.0,
)
tube.write_geo("tube.geo")
tube.write_vtk("tube.vtm")
```

Both primitives follow the same `PrimitiveMesh` contract as
[`BoxMesh`](primitives) / `AnnularCylinderMesh` / `FlatCapsCylinder`:
``self.blocks`` is a one-element list of `HexBlock` with the BCs
already set.

## Examples

End-to-end:

- [`examples/07-h-grid/visualize_h_grid.py`](https://github.com/szaghi/shore/blob/main/examples/07-h-grid/visualize_h_grid.py)
  — Python script generating both variants with optional PyVista viewer.
- [`examples/07-h-grid/visualize_h_grid.sh`](https://github.com/szaghi/shore/blob/main/examples/07-h-grid/visualize_h_grid.sh)
  — bash equivalent, `shore primitive h-grid-rect` + `h-grid-tube` + `shore export`.

## See also

- [Topology catalogue](./) — pick the right topology for your case.
- [`shore primitive`](/reference/primitive) — full CLI reference for all primitives.
- [Chimera assembly walkthrough](/guide/cc-par-pipeline) — how to compose an H-grid background with body-fitted near-mesh blocks.
