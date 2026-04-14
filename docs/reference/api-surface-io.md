# shore.surface.io

STL reader: load a raw triangulated surface into a `trimesh.Trimesh`.

## Functions

### `load_stl`

```python
def load_stl(path: str | Path) -> trimesh.Trimesh
```

Load an STL file and return a watertight-checked `Trimesh`.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `path` | `str \| Path` | Path to the `.stl` file (ASCII or binary) |

**Returns**

`trimesh.Trimesh` — the loaded surface mesh.

**Raises**

| Exception | Condition |
|---|---|
| `FileNotFoundError` | `path` does not exist |
| `ValueError` | The loaded geometry is not a single mesh (e.g. a scene with multiple components) |

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
- The function does not repair the mesh (fill holes, remove degenerate faces). If `project_sphere_to_surface` fails with ray-miss errors, try `mesh.fill_holes()` or `trimesh.repair.fix_winding(mesh)` before projecting.
