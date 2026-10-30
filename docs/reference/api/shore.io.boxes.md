---
title: shore.io.boxes
outline: [2, 3]
---

# `shore.io.boxes`

Hexahedral marker boxes for Xall `cc.par`.

Boxes are the last section of an Xall `cc.par` file (after the
edges section).  Each box defines a hexahedral region of the
simulation domain that marks cells inside it as "internal wall" —
i.e. cells lying inside a solid body that the chimera hole-cutter
should treat as donors-only, not active flow cells.

Conceptually, a box is a closed region that contains the body it
represents and itself fits inside the body's near-mesh: the union
of one body's boxes lies in the thin shell between the body's STL
surface and the outer extent of the body's near-mesh, so that a
*background* mesh's cells overlapping the body land inside at least
one box and get hole-cut.

Box types
---------

The Xall `cc.par` reader (`overset/src/scatole.f90`) supports four
types — parallelepiped (1), single cylinder (2, deprecated), coaxial
cylinders (3), and general hexahedron (9).  SHORE only emits **type
9 (hexahedron)**: it is the universal form, and any axis-aligned or
oriented bounding box can be expressed as 8 vertices.  The
parallelepiped / cylinder shorthands are reader-side conveniences
SHORE has no reason to produce.

Vertex ordering
---------------

Type-9 boxes carry 8 vertices in the order required by Xall's
`esaedro` subroutine (`scatole.f90:30-44`).  The convention is
the standard "VTK voxel" / Fortran-column-major ordering, where the
8 corners enumerate (i, j, k) bits in the order ``v_idx = 1 + i +
2\*j + 4\*k``:

```
1 = (i=lo, j=lo, k=lo)
2 = (i=hi, j=lo, k=lo)
3 = (i=lo, j=hi, k=lo)
4 = (i=hi, j=hi, k=lo)
5 = (i=lo, j=lo, k=hi)
6 = (i=hi, j=lo, k=hi)
7 = (i=lo, j=hi, k=hi)
8 = (i=hi, j=hi, k=hi)

```
The `(i, j, k)` directions need not be world-aligned: a box may be
arbitrarily skewed and rotated so long as the connectivity above is
preserved.  :class:`Box` validates this at construction by
reconstructing the six face centres and rejecting hexahedra whose
faces self-intersect or are degenerate.

Block association
-----------------

Every box references a SHORE block label via :attr:`Box.associated_block`.
At `cc.par` write time the label is resolved to the block's 1-based
index in the adjacency JSON's block list; the box inherits that
block's `level / group / priority` for the chimera precedence
rules in `check_parete_interna`.

## `Box` <Badge type="tip" text="dataclass" /> {#Box}

```python
Box(vertices: np.ndarray = field(), associated_block: str = field())
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/boxes.py#L74">source</a>

One Xall `cc.par` type-9 hexahedral marker box.

**Parameters**

- **`vertices`** — `(8, 3)` array of corner coordinates in the canonical order documented above (see module docstring).  Stored as a `float64` numpy array; lists / tuples are accepted and converted at construction.
- **`associated_block`** — Label of the SHORE block this box is associated with.  At `cc.par` write time this is resolved to a 1-based block index; the box inherits the block's `level / group / priority` for chimera precedence.

**Examples**

Build an axis-aligned cube of half-side 0.6 around the origin,
associated with block `"sub0"`:

```python
>>> import numpy as np
>>> half = 0.6
>>> verts = np.array([
...     [-half, -half, -half],   # v1: (i=lo, j=lo, k=lo)
...     [+half, -half, -half],   # v2: (i=hi, j=lo, k=lo)
...     [-half, +half, -half],   # v3: (i=lo, j=hi, k=lo)
...     [+half, +half, -half],   # v4: (i=hi, j=hi, k=lo)
...     [-half, -half, +half],   # v5: (i=lo, j=lo, k=hi)
...     [+half, -half, +half],   # v6: (i=hi, j=lo, k=hi)
...     [-half, +half, +half],   # v7: (i=lo, j=hi, k=hi)
...     [+half, +half, +half],   # v8: (i=hi, j=hi, k=hi)
... ])
>>> box = Box(vertices=verts, associated_block="sub0")
```

### `Box.from_aabb` <Badge type="info" text="classmethod" /> {#Box-from_aabb}

```python
from_aabb(
    vmin: tuple[float, float, float] | np.ndarray,
    vmax: tuple[float, float, float] | np.ndarray,
    associated_block: str,
) -> Box
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/boxes.py#L137">source</a>

Build a Box from an axis-aligned bounding box.

**Parameters**

- **`vmin, vmax`** — Diagonal corners of the AABB in world coordinates. `vmax > vmin` componentwise.
- **`associated_block`** — Block label.

**Returns**

- `Box` — A type-9 hexahedron whose 8 vertices enumerate the corners of the AABB in canonical order.

## `box_to_json` <Badge type="info" text="function" /> {#box_to_json}

```python
box_to_json(box: Box) -> dict[str, Any]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/boxes.py#L299">source</a>

Serialise a :class:`Box` to a JSON-friendly dict.

## `box_from_json` <Badge type="info" text="function" /> {#box_from_json}

```python
box_from_json(payload: dict[str, Any]) -> Box
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/io/boxes.py#L307">source</a>

Parse one `{associated_block, vertices}` dict into a :class:`Box`.

**Raises**

- `ParameterError` — On missing keys, wrong vertex shape, or an invalid hexahedron (forwarded from :class:`Box.__post_init__`).
