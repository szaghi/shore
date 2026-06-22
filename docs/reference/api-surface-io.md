# shore.surface.io

STL reader: load a raw triangulated surface into a `trimesh.Trimesh`.

## Functions

### `load_stl`

```python
def load_stl(path: str | Path, skip_hygiene: bool = False) -> trimesh.Trimesh
```

Load an STL file and return a hygiene-repaired `Trimesh`.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `path` | `str \| Path` | — | Path to the `.stl` file (ASCII or binary) |
| `skip_hygiene` | `bool` | `False` | When `True`, bypass the mesh-hygiene pass (`merge_vertices`, drop degenerate faces, `fix_normals`) and load the raw geometry. Preserves legacy behaviour for regression tests. |

By default (`skip_hygiene=False`) the loader merges duplicate
vertices, removes degenerate (zero-area) faces, and fixes face
winding via `trimesh.repair.fix_normals`. If the result is still not
watertight it emits a `UserWarning` (it does **not** raise) — ray
casting may then miss rays.

**Returns**

`trimesh.Trimesh` — the loaded (and, unless `skip_hygiene`, repaired) surface mesh.

**Raises**

| Exception | Condition |
|---|---|
| `MeshInputError` | `path` does not exist, or the loaded geometry is not a single `Trimesh` (e.g. a scene with multiple components) |

**Example**

```python
from shore.surface.io import load_stl

mesh = load_stl("fuselage.stl")
print(mesh.vertices.shape)   # (N, 3)
print(mesh.faces.shape)      # (M, 3)
print(mesh.is_watertight)    # True for a clean STL
```

## Notes

- Uses `trimesh.load(..., force="mesh")` to coerce scene files to a single mesh.
- Default loading applies a hygiene pass (merge vertices, drop degenerate faces, fix winding). It does **not** fill holes. If `project_sphere_to_surface` still fails with ray-miss errors on a non-watertight body, try `mesh.fill_holes()` before projecting, or pass `skip_hygiene=True` to inspect the raw geometry.
