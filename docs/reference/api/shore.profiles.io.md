---
title: shore.profiles.io
outline: [2, 3]
---

# `shore.profiles.io`

Selig `.dat` reader and writer for 2D airfoil sections.

The Selig format is the de-facto standard used by the UIUC airfoil
database and most CFD pre-processors:

- An optional first line with the airfoil name (any non-numeric text).
- `N` rows of two whitespace-separated floats: `x y`.
- Coordinates are normalised to chord = 1, ordered upper-TE -&gt; LE
  -&gt; lower-TE.

This module is the single entry / exit point between SHORE's in-memory
`(N, 2)` numpy arrays and on-disk `.dat` files.

## `read_section` <Badge type="info" text="function" /> {#read_section}

```python
read_section(path: str | Path) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/profiles/io.py#L24">source</a>

Load a Selig-format `.dat` airfoil section file.

**Parameters**

- **`path`** — Path to a `.dat` file.  May start with a non-numeric name line (UIUC convention); rows are then `x y` whitespace- separated floats.

**Returns**

- `numpy.ndarray, shape ``(N, 2)``` — Coordinates as read from the file, no reordering.

**Raises**

- `MeshInputError` — If the file is missing, empty, or contains malformed rows.

## `write_section` <Badge type="info" text="function" /> {#write_section}

```python
write_section(path: str | Path, coords: np.ndarray, *, name: str | None = None) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/profiles/io.py#L77">source</a>

Write a 2D section to a Selig-format `.dat` file.

**Parameters**

- **`path`** — Output path.  Parent directory must exist.
- **`coords`** — `(N, 2)` array of section coordinates in Selig order.
- **`name`** — Optional first-line name comment.  Default `None` (no name line written; only coordinate rows).

**Raises**

- `ParameterError` — If `coords` is not `(N, 2)` with `N >= 4`.
