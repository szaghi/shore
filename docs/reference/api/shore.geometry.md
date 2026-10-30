---
title: shore.geometry
outline: [2, 3]
---

# `shore.geometry`

Geometry backends behind one protocol (ADR 0002).

**Modules**

- [`shore.geometry.base`](./shore.geometry.base) — Geometry protocol
- [`shore.geometry.keys`](./shore.geometry.keys) — Persistent entity identity
- [`shore.geometry.occt`](./shore.geometry.occt) — OCCT backend of the `Geometry` protocol
- [`shore.geometry.stl`](./shore.geometry.stl) — STL backend of the `Geometry` protocol

**Re-exported**

- [`BODY`](./shore.geometry.base#BODY)
- `CAD_SUFFIXES`
- [`AnchoredCaster`](./shore.geometry.base#AnchoredCaster)
- [`EntityKind`](./shore.geometry.base#EntityKind)
- [`EntityRef`](./shore.geometry.base#EntityRef)
- [`Geometry`](./shore.geometry.base#Geometry)
- [`Projection`](./shore.geometry.base#Projection)
- [`StlGeometry`](./shore.geometry.stl#StlGeometry)

## `load_geometry` <Badge type="info" text="function" /> {#load_geometry}

```python
load_geometry(path: str | Path, skip_hygiene: bool = False) -> Geometry
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/geometry/__init__.py#L33">source</a>

Load a geometry file, choosing the backend from its suffix.

**Parameters**

- **`path`** — `.stl` (STL backend) or a CAD file (`.step`/`.stp`/`.iges`/ `.igs`/`.brep`/`.brp`, OCCT backend).
- **`skip_hygiene`** — STL only: bypass mesh cleaning (see `shore.surface.io.load_stl`).

**Raises**

- `MeshInputError` — Unknown suffix, unreadable file, or a CAD file while the OCCT backend (`shore-mesh[cad]`) is not installed.
