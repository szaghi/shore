---
title: shore.primitive
outline: [2, 3]
---

# `shore.primitive`

Parametric structured-grid primitives for Chimera/overset assembly.

A *primitive* is a mesh built from simple geometric inputs (corners,
radii, axis) — no STL, no projection, no marching.  Each primitive
returns a :class:`shore.primitive.base.PrimitiveMesh` whose blocks
follow the same :class:`shore.mesh.block.HexBlock` contract used by
body-fitted meshes, so downstream tooling (`.geo` writer, VTK
export, Jacobian inspection) is shared.

Available primitives
--------------------

- :class:`shore.primitive.box.BoxMesh` — single-block axis-aligned
  Cartesian box.  Use case: outer background, simple buffer regions.
- :class:`shore.primitive.annulus.AnnularCylinderMesh` — single-block
  annular cylinder, optionally a partial sector.  Use case: buffer
  around a near-body sphere/airfoil mesh; cylindrical pipe sections.
- :class:`shore.primitive.flat_caps.FlatCapsCylinder` — 5-block
  flat-caps cylinder (central square + 4 trapezoidal sub-blocks).
  Use case: cylindrical background without Chimera fringes (high-order
  schemes, conservation, GPU memory locality).
- :class:`shore.primitive.h_grid_rect.HGridRectMesh` — single-block
  rectangular-duct H-grid: streamwise faces FREE (inflow / outflow),
  the four cross-section walls WALL.  Use case: prismatic ducts,
  internal flow without an internal body.
- :class:`shore.primitive.h_grid_tube.HGridTubeMesh` — single-block
  circular-bore H-grid (lat-lon polar cross-section): streamwise faces
  FREE, j-periodic, centerline singular axis + outer wall WALL.  Use
  case: straight cylindrical ducts, simple SRM chambers.

**Modules**

- [`shore.primitive.annulus`](./shore.primitive.annulus) — Annular cylinder primitive — a single block with optional theta-periodicity
- [`shore.primitive.base`](./shore.primitive.base) — Primitive mesh base class and shared utilities
- [`shore.primitive.box`](./shore.primitive.box) — Cartesian box primitive — a single axis-aligned hexahedral block
- [`shore.primitive.flat_caps`](./shore.primitive.flat_caps) — Flat-caps cylinder primitive — 5-block conforming alternative to a hollow Chimera background
- [`shore.primitive.h_grid_rect`](./shore.primitive.h_grid_rect) — Rectangular-duct H-grid primitive
- [`shore.primitive.h_grid_tube`](./shore.primitive.h_grid_tube) — Circular-bore H-grid primitive

**Re-exported**

- [`AnnularCylinderMesh`](./shore.primitive.annulus#AnnularCylinderMesh)
- [`AnnulusExtent`](./shore.primitive.annulus#AnnulusExtent)
- [`BoxExtent`](./shore.primitive.box#BoxExtent)
- [`BoxMesh`](./shore.primitive.box#BoxMesh)
- [`FlatCapsCylinder`](./shore.primitive.flat_caps#FlatCapsCylinder)
- [`FlatCapsExtent`](./shore.primitive.flat_caps#FlatCapsExtent)
- [`HGridRectMesh`](./shore.primitive.h_grid_rect#HGridRectMesh)
- [`HGridTubeMesh`](./shore.primitive.h_grid_tube#HGridTubeMesh)
- [`PrimitiveMesh`](./shore.primitive.base#PrimitiveMesh)
- [`axis_frame`](./shore.primitive.base#axis_frame)
- [`coords_1d`](./shore.primitive.base#coords_1d)
- [`resolve_axis`](./shore.primitive.base#resolve_axis)
