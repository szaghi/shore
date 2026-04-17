# shore.volume.hyperbolic

Hyperbolic marching: extrude a structured surface layer into a structured volume.

## Functions

### `march`

```python
def march(
    surface: np.ndarray,
    nk: int,
    ds: float,
    growth: float = 1.15,
    smooth_strength: float = 0.2,
    smooth_iters: int = 2,
    periodic: tuple[bool, bool] = (False, True),
    layer_callback: Callable[[int, float, float], None] | None = None,
) -> np.ndarray
```

Extrude a structured surface layer into a volume by hyperbolic marching.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `surface` | `np.ndarray (ni, nj, 3)` | — | k=0 wall layer |
| `nk` | `int` | — | Total layers including k=0 |
| `ds` | `float` | — | First-layer thickness |
| `growth` | `float` | 1.15 | Geometric growth ratio |
| `smooth_strength` | `float` | 0.2 | Laplacian blend factor [0–1] |
| `smooth_iters` | `int` | 2 | Laplacian sweeps per layer |
| `periodic` | `tuple[bool, bool]` | `(False, True)` | `(periodic_i, periodic_j)` — when `True`, the corresponding direction uses modular indexing for end nodes; when `False`, one-sided differences are used. Default `(False, True)` is the standard O-grid (i open, j periodic). Use `(False, False)` for polar cap blocks. |
| `layer_callback` | `Callable[[int, float, float], None] \| None` | None | Called after each accepted layer with `(k, jmin, step)` |

**Returns**

`numpy.ndarray` of shape `(ni, nj, nk, 3)` — the full structured volume grid.

**Raises**

| Exception | Condition |
|---|---|
| `RuntimeError` | A negative-Jacobian cell is detected. Message includes the layer index and the Jacobian value. |

**Example**

```python
from shore.surface.io        import load_stl
from shore.surface.remesh    import project_sphere_to_surface
from shore.volume.hyperbolic import march
from shore.io.geo             import write_geo

mesh = load_stl("sphere.stl")
k0   = project_sphere_to_surface(mesh, ni=40, nj=60)
grid = march(k0, nk=30, ds=1e-3, growth=1.15)

write_geo("sphere.geo", grid)
```

**Custom progress callback:**

```python
def on_layer(k: int, jmin: float, step: float) -> None:
    print(f"k={k:3d}  jmin={jmin:.3e}  step={step:.3e}")

grid = march(k0, nk=30, ds=1e-3, growth=1.15, layer_callback=on_layer)
```

---

### `_compute_normals` *(internal)*

```python
def _compute_normals(
    layer: np.ndarray,
    periodic: tuple[bool, bool] = (False, True),
) -> np.ndarray
```

Estimate outward unit normals for every node of a structured layer. Central differences in both directions; one-sided at non-periodic boundaries. The `periodic` flag controls which directions wrap around.

---

### `_laplacian_smooth` *(internal)*

```python
def _laplacian_smooth(
    layer: np.ndarray,
    strength: float,
    periodic: tuple[bool, bool] = (False, True),
) -> np.ndarray
```

One step of in-plane Laplacian smoothing. Non-periodic boundary rows/columns are clamped (never moved). The k=0 wall row is always clamped.

::: tip Internal functions
Functions prefixed with `_` are internal implementation details. Their signatures may change without notice. Use `march` for all production code.
:::

## Jacobian checking

The marching loop uses `min_jacobian_two_layers` from [`shore.volume.jacobian`](api-jacobian) to validate each new layer before accepting it. The Jacobian computation has been extracted into a dedicated module so that `shore check`, `shore export`, and user code can all use the same implementation.
