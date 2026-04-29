# Per-edge spacing

A `HexBlock` has **12 axis-parallel edges** — 4 along each of i, j,
k. By default, all 4 parallel edges of an axis share the same
`Spacing1D` law; the block's interior is then the tensor product of
the three per-axis distributions. **Per-edge spacing** breaks that
symmetry: replace the spacing on a subset of parallel edges and the
block builder switches to a transfinite-interpolation (TFI) blend
between the four edge distributions.

This guide covers when to use the feature, the three example
scripts that demonstrate it, and the precedence rules SHORE uses to
keep multi-block meshes consistent.

For the data model, see [`Spacing1D`](/reference/api-spacing); for
the underlying laws, see [Spacing laws](/algorithm/spacing-laws).

## When to reach for per-edge spacing

The classic case is a **boundary-layer-style request on a parametric
background**: "uniform everywhere except near the body, where I
want clustering on one side." The legacy tensor-product API forces
you to either accept uniform-everywhere or restructure the block.
Per-edge spacing expresses the request as a single edge override,
and TFI takes care of the rest.

Three flavours of per-edge support exist in SHORE:

| Topology / primitive | Per-edge i | Per-edge j | Per-edge k |
|---|---|---|---|
| `BoxMesh` | yes (TFI) | yes (TFI) | yes (TFI) |
| `AnnularCylinderMesh` (sector) | yes (TFI) | yes (TFI) | yes (TFI) |
| `AnnularCylinderMesh` (full annulus) | yes | rejected (j periodic) | yes |
| `FlatCapsCylinder` | rejected (geometric) | rejected (geometric) | yes (TFI) |
| `CubedSphereTopology` | rejected (geometric) | rejected (geometric) | yes (identity-check) |
| `OGridTopology` | rejected (j periodic / polar) | rejected (j periodic) | yes (identity-check) |

Where it's "rejected", the construction is pinned by geometric
constraints (cap / equator C1 matching, periodic-j, polar singular
axes) that override per-edge control would invalidate. Where it's
"yes", it's the right tool.

## Edge naming

Each edge's name encodes the **two constant axes** along it:

- i-axis edges: `j_lo_k_lo`, `j_hi_k_lo`, `j_lo_k_hi`, `j_hi_k_hi`
- j-axis edges: `i_lo_k_lo`, `i_hi_k_lo`, `i_lo_k_hi`, `i_hi_k_hi`
- k-axis edges: `i_lo_j_lo`, `i_hi_j_lo`, `i_lo_j_hi`, `i_hi_j_hi`

So `j_lo_k_lo` runs along i with j=0, k=0 — it's the bottom-front
i-edge of the block.

## Recipe 1 — per-edge i / j on `BoxMesh`

The full example is at `examples/03-per-edge-spacing/per_edge_ij_box.py`. The
mechanism: pass `edge_pins={edge_name: Spacing1D, ...}` to
`BoxMesh.build`. Any edge of an axis that's listed in `edge_pins`
overrides that axis's default; once any edge of an axis is pinned,
the block builder switches that axis to TFI between the 4 parallel
edges.

```python
from shore.primitive import BoxMesh
from shore.mesh.spacing import Spacing1D

# Tensor-product baseline: uniform everywhere
baseline = BoxMesh.build(
    vmin=(0.0, 0.0, 0.0), vmax=(10.0, 5.0, 3.0),
    ni=40, nj=30, nk=20,
)

# Per-edge i-pin: cluster the bottom-front i-edge near the wall
pinned = BoxMesh.build(
    vmin=(0.0, 0.0, 0.0), vmax=(10.0, 5.0, 3.0),
    ni=40, nj=30, nk=20,
    edge_pins={"j_lo_k_lo": Spacing1D(law="tanh", beta=4.0)},
)
```

In the result, the `(j=0, k=0)` i-row is densely packed near
`vmin[0]`; opposite corners stay uniform; interior rows blend the
four. See the example script for a print-out that compares
coordinates at the pinned corner, opposite corner, and midpoint.

The same mechanism works for j-edges (override e.g.
`i_lo_k_lo`) and k-edges (override e.g. `i_lo_j_lo`). Mix axes
freely.

## Recipe 2 — per-edge k on `CubedSphereTopology`

The full example is at `examples/03-per-edge-spacing/per_edge_k_spacing_cubed_sphere.py`.

The cubed-sphere mesh has 6 blocks with k-coupled seams (cap /
equator DIRICHLET pins at the i-boundary, sub-sub SHARED views at
the j-boundary), so all 6 blocks must use a single k-schedule.
SHORE detects per-edge overrides by **identity check**: replace the
4 k-edges of every block with the same `Spacing1D` instance, and
`Mesh.march()` reads that schedule instead of its `(ds, growth)`
arguments.

```python
from shore.mesh import CubedSphereTopology, Mesh
from shore.mesh.spacing import Spacing1D

mesh = Mesh.from_stl("sphere.stl",
                     topology=CubedSphereTopology(),
                     ni=40, nj=60, nk=40)

# One Spacing1D, shared across every k-edge of every block
sp_k = Spacing1D(law="tanh", beta=4.0, total_thickness=2.0)
for blk in mesh.blocks:
    for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
        blk.edge(name).spacing = sp_k

mesh.march()                  # tanh schedule used; ds/growth defaults ignored
mesh.write_geo("sphere")
```

### Precedence rules

`Mesh.march()` resolves the k-step schedule by precedence:

1. The `steps=` argument (highest).
2. Block 0's k-edges, if all 4 carry a user-overridden `Spacing1D`
   that agrees, **and** every other block's k-edges agree.
3. Block 0's `spacing_k` constructor argument.
4. The `(ds, growth)` defaults (lowest).

Disagreement at level (2) is not silent. Two cases emit
`UserWarning` and fall back to the next level:

- **Partial override** — a block's 4 k-edges are *some* user
  instance and *some* the constructor default. The block falls back
  to its default.
- **Inter-block disagreement** — block 0 wants schedule A but block
  N wants B. Block 0 wins (the canonical choice).

These warnings catch the easy mistake — "I overrode the equatorial
sub-blocks but forgot the caps" — without breaking the build.

## Recipe 3 — per-edge k on `OGridTopology`

The full example is at `examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py`.

Identical mechanism to the cubed sphere, single block. Set the
4 k-edges of `mesh.blocks[0]` to a shared `Spacing1D` and
`Mesh.march()` picks it up:

```python
from shore.mesh import Mesh, OGridTopology
from shore.mesh.spacing import Spacing1D

mesh = Mesh.from_stl("sphere.stl", topology=OGridTopology(),
                     ni=40, nj=60, nk=40)

sp_k = Spacing1D(law="tanh2", beta=3.0, total_thickness=2.0)
for name in ("i_lo_j_lo", "i_hi_j_lo", "i_lo_j_hi", "i_hi_j_hi"):
    mesh.blocks[0].edge(name).spacing = sp_k

mesh.march()
```

Per-edge i (the meridional direction) is rejected on `OGridTopology`
because the j-direction is periodic — the four parallel i-edges
sit at *adjacent* azimuths under the j-wrap rather than at distinct
transverse corners. Use the per-axis `i_spacing` / `i_beta`
arguments on `from_stl()` instead.

## What about `first_step` / `last_step`?

[`Edge`](/reference/api-spacing#edge) has two extra fields,
`first_step` and `last_step`, that request **C1 endpoint matching**:
the first / last interval along the edge takes the chosen physical
length, with the interior shaped by `edge.spacing`.

This is the kernel that `CubedSphereTopology` uses internally to
pin every meridian's first cell to the cap-edge first cell at the
seam vertex (the **per-j seam pinning** discussed in
[k=0 mesh construction](/algorithm/k0-construction)). The
`Spacing1D(law="pinned", first_step=..., last_step=...)` law
exposes that machinery as user-facing API: any topology or
primitive that wants C1 endpoint matching at multi-block edges can
read `edge.first_step` / `edge.last_step` and apply it.

## See also

- [`Spacing1D`](/reference/api-spacing) — per-edge dataclass.
- [`shore.volume.spacing`](/reference/api-volume-spacing) —
  low-level laws.
- [Spacing laws](/algorithm/spacing-laws) — algorithmic overview.
- [`BoxMesh`, `AnnularCylinderMesh`, `FlatCapsCylinder`](/reference/api-primitives) — primitives that accept `edge_pins`.
- Example scripts: `examples/03-per-edge-spacing/per_edge_ij_box.py`,
  `examples/03-per-edge-spacing/per_edge_k_spacing_cubed_sphere.py`,
  `examples/03-per-edge-spacing/per_edge_k_spacing_ogrid.py`.
