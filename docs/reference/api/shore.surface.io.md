---
title: shore.surface.io
outline: [2, 3]
---

# `shore.surface.io`

STL reader: load a raw triangulated surface into a trimesh.Trimesh.

## `load_stl` <Badge type="info" text="function" /> {#load_stl}

```python
load_stl(path: str | Path, skip_hygiene: bool = False) -> trimesh.Trimesh
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/io.py#L16">source</a>

Load an STL file and return a trimesh.Trimesh.

**Parameters**

- **`path`** — Path to the .stl file (ASCII or binary).
- **`skip_hygiene`** — When *True*, bypass mesh hygiene steps (merge_vertices, remove_degenerate_faces, fix_normals).  Preserves legacy behaviour for regression tests.  Default *False*.

**Returns**

- `trimesh.Trimesh` — The loaded (and optionally repaired) surface mesh.

**Raises**

- `MeshInputError` — If *path* does not exist or the loaded geometry is not a single mesh.
