---
title: shore.surface.anchor
outline: [2, 3]
---

# `shore.surface.anchor`

Anchor-point resolution strategies for ray-cast projections.

Phase 1: resolve_anchor with centroid / bbox / incenter / user / auto strategies.

## `resolve_anchor` <Badge type="info" text="function" /> {#resolve_anchor}

```python
resolve_anchor(
    mesh: trimesh.Trimesh,
    strategy: str = 'auto',
    user_point: np.ndarray | None = None,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/surface/anchor.py#L20">source</a>

Resolve a validated ray-cast anchor point for *mesh*.

**Parameters**

- **`mesh`** — Triangulated body surface (already hygiene-processed is recommended).
- **`strategy`** — `'centroid'` — volumetric centroid (legacy; may lie outside on non-convex bodies). `'bbox'` — bounding-box centre (always inside a convex body). `'incenter'` — centre of the minimum enclosing sphere (`trimesh.nsphere.minimum_nsphere`); robust for smooth convex bodies. `'user'` — caller-supplied *user_point*; validated before use. `'auto'` — try incenter, validate, fall back to bbox, validate, raise.
- **`user_point`** — Required when *strategy* is `'user'`.

**Returns**

- `np.ndarray, shape (3,)` — Validated anchor point confirmed inside the mesh.

**Raises**

- `ParameterError` — If *strategy* is unknown or *user_point* is absent for `'user'`.
- `MeshInputError` — If the resolved anchor fails the containment check, or if `'auto'` exhausts all fallbacks.
