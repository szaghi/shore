---
title: shore.io.geo
outline: [2, 3]
---

# `shore.io.geo`

Xall geo format writer/reader.

Format (ASCII, one block per file)
-----------------------------------
ni nj nk
x(0,0,0)  x(1,0,0)  ...   (all X, i-major, then j, then k)
y(0,0,0)  ...
z(0,0,0)  ...

Coordinate ordering matches Fortran column-major: i is the fastest index.
The grid array passed in must have shape (ni, nj, nk, 3) with the last
axis being (x, y, z).

## `write_geo` <Badge type="info" text="function" /> {#write_geo}

```python
write_geo(path: str | Path, grid: np.ndarray) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/geo.py#L24">source</a>

Write a single structured block to an ASCII geo file.

**Parameters**

- **`path`** — Output file path.
- **`grid`** — Coordinates array, shape (ni, nj, nk, 3).

## `read_geo` <Badge type="info" text="function" /> {#read_geo}

```python
read_geo(path: str | Path) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/geo.py#L47">source</a>

Read a single structured block from an ASCII geo file.

**Parameters**

- **`path`** — Path to the .geo file.

**Returns**

- `numpy.ndarray, shape (ni, nj, nk, 3)`
