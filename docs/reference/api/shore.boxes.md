---
title: shore.boxes
outline: [2, 3]
---

# `shore.boxes`

Generate :class:`shore.io.boxes.Box` instances from an STL surface.

The Xall `cc.par` boxes section needs hexahedra that *cover the
body interior*: every cell of the background mesh whose centre falls
inside the solid body must lie inside at least one box, so the
chimera hole-cutter can mark it as donors-only.  See
:mod:`shore.io.boxes` for the role boxes play.

This module turns a triangulated body into a set of such hexahedra.
The default method (`"voxel-fill"`) voxelises the body's AABB at
a resolution chosen from `max_outward_gap` and emits each voxel
inside the body as an axis-aligned :class:`Box`.  Coverage is
guaranteed by construction; the only failure mode is excessive box
count for tiny gap budgets relative to the body size.

Methods
-------
- `"voxel-fill"` (default) — voxelise the body's AABB at a
  resolution such that every body-overlapping voxel has corners
  within `max_outward_gap` of the body surface (outward only).
  Each voxel that overlaps the body becomes one axis-aligned
  :class:`Box`; together they form a watertight cover of the
  body (interior + surface shell), suitable for chimera
  hole-cutting.
- `"aabb"`  — single world-axis-aligned bounding box of the body.
  Best effort: emits a `UserWarning` when the measured outward
  gap exceeds `max_outward_gap`.  Useful as a quick preview only.
- `"obb"`   — single oriented bounding box (area-weighted PCA).
  Best effort with the same warning behaviour as `aabb`.

Outward gap (best-effort methods only)
--------------------------------------
`aabb` and `obb` measure the actual distance from sampled
box-surface points (8 corners + a parametric `kxk` grid on each
face) to the STL surface (outward only) and compare to
`max_outward_gap`.  `voxel-fill` enforces the gap by its voxel-
size choice and does not need to measure.

Watertightness
--------------
`voxel-fill` calls `trimesh.contains()` which assumes a
watertight body STL.  Non-watertight inputs return unreliable
inside/outside classifications; :func:`shore.surface.io.load_stl`
already warns on non-watertight meshes, and :func:`boxes_from_stl`
re-checks and raises if voxel-fill is requested on one.

API
---
.. autofunction:: boxes_from_stl

## `boxes_from_stl` <Badge type="info" text="function" /> {#boxes_from_stl}

```python
boxes_from_stl(
    stl_path: str | Path,
    associated_block: str,
    *,
    max_outward_gap: float,
    method: Method = 'voxel-fill',
    face_samples_per_axis: int = 4,
) -> list[Box]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/boxes.py#L72">source</a>

Build a list of :class:`Box` covering the STL body interior.

**Parameters**

- **`stl_path`** — Path to the body's STL file.  Loaded via :func:`shore.surface.io.load_stl` (with hygiene).
- **`associated_block`** — SHORE block label every returned box references.  Resolved to a 1-based index by :func:`shore.io.cc_par.write_cc_par`.
- **`max_outward_gap`** — Maximum allowed distance, in STL world units, from any point on a box's surface to the nearest STL surface point — measured **outward only** (points on the box surface that lie inside the body don't count).  Strictly positive.  Should be smaller than the body's near-mesh wall-normal thickness so the boxes fit inside the wall.
- **`method`** — Algorithm: `"voxel-fill"` (default — guaranteed body-volume coverage by voxelising the AABB), `"aabb"` / `"obb"` (single-box, best-effort).
- **`face_samples_per_axis`** — `"aabb"` / `"obb"` post-check sampling density.  Each of the 6 box faces is sampled on a `kxk` grid (total = 8 corners + 6\*k\*k samples per box) for the outward-gap measurement.  Default 4 (104 samples/box).  Ignored by `voxel-fill`.

**Returns**

- `list[Box]` — One or more boxes, all sharing `associated_block`.

**Raises**

- `ParameterError` — If parameters are out of range, or `method` is unknown.
- `MeshInputError` — Forwarded from :func:`load_stl`, or raised when `method='voxel-fill'` on a non-watertight STL.
