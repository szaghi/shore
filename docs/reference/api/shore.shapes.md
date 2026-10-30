---
title: shore.shapes
outline: [2, 3]
---

# `shore.shapes`

3D body STL builders backed by trimesh.

This module produces `trimesh.Trimesh` objects (and writes `.stl`
files via the standard `mesh.export(path)` call).  Two families:

- **Analytic primitives** for testing and quick CFD studies:
  :func:`sphere`, :func:`cylinder`, :func:`torus`, :func:`box`,
  :func:`capsule`.  Thin wrappers over `trimesh.creation`.

- **Section extrusion**: :func:`extrude_section` takes a 2D
  Selig-ordered profile (e.g. from :mod:`shore.profiles`) plus
  span / taper / sweep / twist parameters and produces a wing STL.

For the Selig 2D format that :func:`extrude_section` consumes, see
:mod:`shore.profiles`.

## `box` <Badge type="info" text="function" /> {#box}

```python
box(extents: tuple[float, float, float] = (1.0, 1.0, 1.0)) -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/shapes/_primitives.py#L119">source</a>

Axis-aligned rectangular box centred at the origin.

**Parameters**

- **`extents`** — `(dx, dy, dz)` total side lengths.  All must be positive.

**Returns**

- `trimesh.Trimesh` — 12-triangle watertight box.

## `capsule` <Badge type="info" text="function" /> {#capsule}

```python
capsule(radius: float = 0.5, height: float = 1.0, sections: int = 16) -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/shapes/_primitives.py#L137">source</a>

Capsule (cylinder with hemispherical caps), axis along z.

**Parameters**

- **`radius`** — Hemispherical cap radius and cylinder radius.
- **`height`** — Length of the cylindrical middle section (caps add `2 * radius` to total length).
- **`sections`** — Segments around the circumference; the same value is used for both circumferential and meridional resolution on the caps.

**Returns**

- `trimesh.Trimesh` — Watertight capsule.

## `cylinder` <Badge type="info" text="function" /> {#cylinder}

```python
cylinder(
    radius: float = 1.0,
    height: float = 1.0,
    sections: int = 32,
) -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/shapes/_primitives.py#L40">source</a>

Right circular cylinder centred at the origin, axis along z.

**Parameters**

- **`radius`** — Cylinder radius.  Must be positive.
- **`height`** — Total length along the z-axis (the cylinder spans `[-height/2, +height/2]`).
- **`sections`** — Number of segments around the circumference.  Higher = smoother.

**Returns**

- `trimesh.Trimesh` — Watertight cylinder with circular caps.

## `extrude_section` <Badge type="info" text="function" /> {#extrude_section}

```python
extrude_section(
    section_2d: np.ndarray,
    *,
    span: float,
    n_stations: int,
    taper: float = 1.0,
    sweep_deg: float = 0.0,
    twist_deg: float = 0.0,
    span_axis: str = 'z',
    end_caps: bool = False,
) -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/shapes/_extrude.py#L42">source</a>

Extrude a 2D section along a span axis into a 3D wing.

**Parameters**

- **`section_2d`** — `(N, 2)` Selig-ordered closed-or-open section polyline. For closed sections (sharp TE) `section_2d[0]` and `section_2d[-1]` should coincide; for open / blunt TE they differ.  No reordering is applied.
- **`span`** — Total wing span (length along `span_axis`).  Must be positive.
- **`n_stations`** — Number of span stations including root and tip.  `>= 2`.
- **`taper`** — Tip-to-root chord ratio.  `1.0` (default) = constant chord; `0.5` = tip chord half the root chord.
- **`sweep_deg`** — Sweep angle of the leading edge, in degrees, about the thickness axis.  Positive = sweepback.  `0.0` (default).
- **`twist_deg`** — Geometric twist (washout / washin) at the tip, in degrees, about the LE-thickness axis.  Linearly interpolated from `0.0` at the root.  Positive = nose-up.
- **`span_axis`** — `"x"`, `"y"`, or `"z"`.  Default `"z"`.
- **`end_caps`** — If `True`, add a centroid-fan cap at each span end so the wing is watertight.  Default `False` — the SHORE C-grid topology relies on open ends to slice cleanly.

**Returns**

- `trimesh.Trimesh` — Wing STL.

**Raises**

- `ParameterError` — On invalid arguments.

## `sphere` <Badge type="info" text="function" /> {#sphere}

```python
sphere(radius: float = 1.0, subdivisions: int = 3) -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/shapes/_primitives.py#L15">source</a>

Icosahedral sphere centred at the origin.

**Parameters**

- **`radius`** — Sphere radius.  Must be positive.
- **`subdivisions`** — Icosphere subdivision level.  `0` = 20 triangles (raw icosahedron); `3` (default) = 1280 triangles; each level quadruples the count.  `4` is dense enough for most CFD wraps.

**Returns**

- `trimesh.Trimesh` — Watertight icosahedral sphere.

## `torus` <Badge type="info" text="function" /> {#torus}

```python
torus(
    r_major: float = 1.0,
    r_minor: float = 0.25,
    major_sections: int = 32,
    minor_sections: int = 32,
) -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/shapes/_primitives.py#L72">source</a>

Torus centred at the origin, axis of revolution = z.

**Parameters**

- **`r_major`** — Major radius (centre of the tube to centre of the torus).
- **`r_minor`** — Minor radius (tube radius).  Must be smaller than `r_major` for a self-non-intersecting torus.
- **`major_sections`** — Segments around the major circle.
- **`minor_sections`** — Segments around the tube (minor circle).

**Returns**

- `trimesh.Trimesh` — Watertight ring torus (genus-1 surface).
