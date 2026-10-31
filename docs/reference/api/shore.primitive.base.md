---
title: shore.primitive.base
outline: [2, 3]
---

# `shore.primitive.base`

Primitive mesh base class and shared utilities.

A *primitive* is a parametric structured grid built from simple
geometric inputs (corners, radii, axis) — no STL, no projection, no
hyperbolic marching.  Primitives are the connective tissue of an
overset assembly: backgrounds, buffers, and the conforming alternative
to Chimera fringes.

Each concrete primitive is a subclass of :class:`PrimitiveMesh` that
populates `self.blocks` with one or more :class:`HexBlock` instances
at construction time.  Downstream operations (Jacobian check, `.geo`
write, VTK export) are inherited.

## `PrimitiveMesh` <Badge type="tip" text="class" /> {#PrimitiveMesh}

```python
class PrimitiveMesh(ABC)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L354">source</a>

Abstract base for parametric structured grids.

Concrete subclasses populate `self.blocks` (a list of
:class:`HexBlock`) and `self.label_default` (the file name to
use when `write_geo(stem)` is called and the writer must pick a
suffix for multi-block primitives).

### `PrimitiveMesh.write_geo` <Badge type="info" text="method" /> {#PrimitiveMesh-write_geo}

```python
write_geo(output: str | Path) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L371">source</a>

Write the primitive to one or more `.geo` files.

Single-block primitives write *output* verbatim (a single file).
Multi-block primitives append `_<label>.geo` to *output*'s stem.

### `PrimitiveMesh.write_vtk` <Badge type="info" text="method" /> {#PrimitiveMesh-write_vtk}

```python
write_vtk(output: str | Path, with_jacobian: bool = True) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L392">source</a>

Write the primitive as a single `.vtm` MultiBlock file.

Single- and multi-block primitives both produce one `.vtm`;
the difference is just the number of children inside it.

## `axis_field_tfi` <Badge type="info" text="function" /> {#axis_field_tfi}

```python
axis_field_tfi(
    n_axis: int,
    n_t1: int,
    n_t2: int,
    start: float,
    stop: float,
    edge_specs: dict[tuple[int, int], Spacing1D],
    axis_label: str,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L167">source</a>

Return `(n_axis, n_t1, n_t2)` axis-coordinate field via TFI.

The 4 parallel edges of an axis are keyed by their endpoint indices in
the two transverse directions:

- `(0, 0)` — both transverse axes at lo (e.g. j_lo_k_lo for axis i)
- `(1, 0)` — first transverse hi, second lo (j_hi_k_lo)
- `(0, 1)` — first lo, second hi              (j_lo_k_hi)
- `(1, 1)` — both hi                          (j_hi_k_hi)

Each `edge_specs` value is a `Spacing1D` describing the node
distribution along that edge.  All 4 keys must be present.

**Parameters**

- **`n_axis`** — Number of nodes along the axis being built.
- **`n_t1, n_t2`** — Number of nodes along the two transverse axes.
- **`start, stop`** — World coordinates of the axis interval.
- **`edge_specs`** — Mapping `{(0,0), (1,0), (0,1), (1,1) → Spacing1D}`.
- **`axis_label`** — Used in error messages (`'i'` / `'j'` / `'k'`).

**Returns**

- `numpy.ndarray, shape ``(n_axis, n_t1, n_t2)```

## `axis_frame` <Badge type="info" text="function" /> {#axis_frame}

```python
axis_frame(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L85">source</a>

Return an orthonormal right-handed frame `(e1, e2, axis)` for *axis*.

*axis* is assumed to be a unit vector.  `e1` and `e2` span the
plane perpendicular to *axis*; the choice of `e1` is deterministic
but arbitrary (whichever world axis is "least parallel" to *axis*).

## `coords_1d` <Badge type="info" text="function" /> {#coords_1d}

```python
coords_1d(
    n: int,
    start: float,
    stop: float,
    law: str = 'uniform',
    beta: float = 3.0,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L108">source</a>

Build *n* sample points on `[start, stop]` using *law*.

*law* is one of `"uniform"`, `"tanh"` (one-sided, dense at
*start*), or `"tanh2"` (two-sided, dense at both ends).

**Raises**

- `ValueError` — If *law* is unrecognised or *n* &lt; 2.

## `resolve_axis` <Badge type="info" text="function" /> {#resolve_axis}

```python
resolve_axis(spec: str | tuple[float, float, float] | np.ndarray) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L39">source</a>

Resolve an axis specifier into a unit vector in R^3.

**Parameters**

- **`spec`** — Either a single-character axis label (`'x'`, `'y'`, `'z'`, case-insensitive) or a 3-component sequence interpreted as a free axis direction (will be normalised).

**Returns**

- `numpy.ndarray, shape (3,)` — Unit vector in the chosen direction.

**Raises**

- `ValueError` — If *spec* is neither a known label nor a 3-vector, or if a 3-vector has zero magnitude.

## `resolve_edge_specs` <Badge type="info" text="function" /> {#resolve_edge_specs}

```python
resolve_edge_specs(
    axis: str,
    edge_pins: dict[str, Spacing1D] | None,
    axis_default: Spacing1D,
) -> dict[tuple[int, int], Spacing1D]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L237">source</a>

Resolve the 4 parallel-edge `Spacing1D`s for one axis.

**Parameters**

- **`axis`** — `'i'` / `'j'` / `'k'`.
- **`edge_pins`** — Optional sparse user override: `{edge_name: Spacing1D}`.  Edges not in the dict (or unmentioned for this axis) fall back to *axis_default*.  Edge names are the canonical 12 in :mod:`shore.mesh.edge`.
- **`axis_default`** — `Spacing1D` used for any edge not pinned in *edge_pins*.

**Returns**

- `dict ``{(0,0), (1,0), (0,1), (1,1) → Spacing1D}```

**Notes**

Edge-name → axis-key mapping (from :mod:`shore.mesh.block`):

- axis `i` (j and k constant): j_lo_k_lo=(0,0), j_hi_k_lo=(1,0),
  j_lo_k_hi=(0,1), j_hi_k_hi=(1,1)  — first key = j, second = k
- axis `j` (i and k constant): i_lo_k_lo=(0,0), i_hi_k_lo=(1,0),
  i_lo_k_hi=(0,1), i_hi_k_hi=(1,1)  — first key = i, second = k
- axis `k` (i and j constant): i_lo_j_lo=(0,0), i_hi_j_lo=(1,0),
  i_lo_j_hi=(0,1), i_hi_j_hi=(1,1)  — first key = i, second = j

## `stamp_edge_specs_onto_block` <Badge type="info" text="function" /> {#stamp_edge_specs_onto_block}

```python
stamp_edge_specs_onto_block(
    block: HexBlock,
    *,
    edge_pins: dict[str, Spacing1D] | None,
) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/primitive/base.py#L301">source</a>

Populate `block.edge(name).spacing` for every name in *edge_pins*.

After construction, downstream readers (e.g. `Mesh._resolve_k_steps`)
may want to learn the actual per-edge schedule used.  This stamps the
user's pins onto the block so they survive past `build()`.
