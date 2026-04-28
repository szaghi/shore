# shore.io.geo

Xall geo format reader and writer.

## File format

The geo format is a plain ASCII file containing a single structured block:

```
ni nj nk
x(0,0,0)  x(1,0,0)  x(2,0,0)  ...    (all x values, i-major Fortran order)
y(0,0,0)  y(1,0,0)  ...
z(0,0,0)  ...
```

**Coordinate ordering** is Fortran column-major: the i-index is fastest (innermost), then j, then k. This matches the memory layout of Fortran arrays declared as `real :: x(ni, nj, nk)`.

Six values are written per line for readability. Coordinates are written with 10 significant digits in scientific notation (`1.0000000000e+00`).

**Example** for a 2×2×2 grid:

```
2 2 2
1.0000000000e+00  2.0000000000e+00  1.0000000000e+00  2.0000000000e+00
1.0000000000e+00  2.0000000000e+00  1.0000000000e+00  2.0000000000e+00
0.0000000000e+00  0.0000000000e+00  1.0000000000e+00  1.0000000000e+00
0.0000000000e+00  0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
0.0000000000e+00  0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
1.0000000000e+00  1.0000000000e+00  1.0000000000e+00  1.0000000000e+00
```

## Functions

### `write_geo`

```python
def write_geo(path: str | Path, grid: np.ndarray) -> None
```

Write a single structured block to an ASCII geo file.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `path` | `str \| Path` | Output file path (created or overwritten) |
| `grid` | `np.ndarray (ni, nj, nk, 3)` | Coordinate array; last axis is (x, y, z) |

**Example**

```python
from shore.io.geo import write_geo
import numpy as np

# Create a trivial 4×5×6 grid
ni, nj, nk = 4, 5, 6
i = np.arange(ni, dtype=float)
j = np.arange(nj, dtype=float)
k = np.arange(nk, dtype=float)
ig, jg, kg = np.meshgrid(i, j, k, indexing="ij")
grid = np.stack([ig, jg, kg], axis=-1)

write_geo("grid.geo", grid)
```

---

### `read_geo`

```python
def read_geo(path: str | Path) -> np.ndarray
```

Read a single structured block from an ASCII geo file.

**Parameters**

| Name | Type | Description |
|---|---|---|
| `path` | `str \| Path` | Path to the `.geo` file |

**Returns**

`numpy.ndarray` of shape `(ni, nj, nk, 3)`.

**Raises**

| Exception | Condition |
|---|---|
| `FileNotFoundError` | `path` does not exist |
| `ValueError` | The number of values in the file does not match `3 × ni × nj × nk` |

**Example**

```python
from shore.io.geo import read_geo

grid = read_geo("sphere.geo")
ni, nj, nk, _ = grid.shape
print(f"Grid: {ni}×{nj}×{nk}")
```

## Round-trip guarantee

`write_geo` followed by `read_geo` on the same path is guaranteed to reproduce the original array within floating-point precision (`numpy.testing.assert_allclose` with default tolerances). This is verified by `tests/test_geo.py::test_round_trip`.
