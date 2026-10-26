# 0015 — STL normals for the marcher's normal blending

- **Status:** Accepted
- **Date:** 2026-10-26

## Context

ADR 0014 made `StlGeometry`'s normals crease-aware. It stated that the
legacy `shore mesh` path (`volume.hyperbolic._stl_normals`) was untouched
and that existing meshes stayed bit-identical. **That statement was
wrong.** Since 8893973 (CAD input), the legacy marchers (cubed-sphere,
O-grid, C-grid, OH/CH) read STL normals through `Body.normals`, which
called `StlGeometry.project`. `_stl_normals` had no caller left. So the
ADR 0014 commit changed every legacy march with `blend_normals_k > 0` on
an STL with sharp edges. With the default `blend_normals_k = 0`, STL
normals are never read, and meshes were unchanged; that is all the
existing tests covered.

Measured on a 96-section cylinder STL (ds 0.01, growth 1.15, blend over 5
layers):

| Topology | Old normals (vertex blend) | Crease-aware normals |
|---|---|---|
| O-grid | marches; first layer up to 50.7° off the true normal (21.9° unblended) | marches; 18.4° |
| cubed-sphere | marches; 51.5° | **`MarchingError` at k = 1, for blend windows of 2, 5 and 10** |

The fold is at the cylinder's convex 90° rim. Neighbouring grid rows
0.016 apart lie on either side of it: the side row marches outward, the
cap row straight down. The blend applies the STL normals at full weight
at the wall, so the 90° jump between neighbouring marching directions,
with the layer smoothing, inverts 7 cells. The old normals averaged
across the rim. They were wrong as surface normals (up to 45-55°) but
robust as marching directions.

## Decision

User decision of 2026-10-26: **restore the old normals for the blend.**
`Body.normals` returns `_stl_normals` (trimesh's vertex normals, blended
over the closest triangle) for STL bodies, and the exact geometry normals
for CAD bodies, as before ADR 0014. It is bit-identical to the
pre-ADR-0014 `StlGeometry` normals, checked on a box and a cylinder.

The two kinds of normal answer different questions:
- crease-aware normals (`StlGeometry`, ADR 0014) answer "what is the
  surface normal here", the question of the protocol kernels (the tangent
  plane of `march_surface`, the direction of the intersection tracer);
- a marching direction must vary smoothly across a sharp edge, the
  question of the blend.

`_stl_normals` stays (it is not dead code, contrary to the first plan).

## Verification

Unit tests (`tests/test_geometry.py`, no spec: the blend's normals are the
pre-ADR-0014 ones, bit for bit):
- `Body.normals` on an STL box equals `_stl_normals`, bit for bit;
- the blended cubed-sphere march on the cylinder completes.

Both fail without the change. The ADR 0014 spec still reports 17 passed;
full suite: 1414 passed, 0 xfailed.

**Out of scope:**
- an edge-band treatment that gives correct *and* smooth marching
  directions across sharp edges (bisector at the edge, blended over a
  band). That would let the blend use correct normals, and needs its own
  spec.
