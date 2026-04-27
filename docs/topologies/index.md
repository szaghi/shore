# Topologies

A **topology** in SHORE is the way the structured mesh is decomposed
into hexahedral blocks and the rules for stitching those blocks
together at their shared faces, edges, and vertices. The
[algorithm](/algorithm/) section describes the building blocks
(projection, k=0 construction, hyperbolic marching, seam-aware normals,
spacing laws); this section catalogues the **assemblies** built from
them.

## Catalogue

| Topology | Status | Blocks | Use case |
|---|---|---|---|
| [O-grid](o-grid) | first-class | 1 (single, j-periodic) | `shore mesh --topology ogrid`; Chimera background, downstream tools that expect a single periodic-j block |
| [Cubed sphere (6-block)](cubed-sphere) | **production default** | 6 (4 sub-blocks + 2 caps) | `shore mesh` (default); clean polar coverage with no singular axis |
| [Yin-Yang](yin-yang) | Python-API only | 2 (overlapping) | Specialised; not wired to the CLI |

Future topologies that fit naturally into this framework:

- **C-grid** for airfoils and wings with sharp trailing edges
- **H-grid** for bluff bodies and internal flows
- **Multi-body Chimera assemblies** for compound geometries

Each new topology gets one page here; the algorithm pages stay shared.

## How to choose

The CLI defaults to the **cubed-sphere topology** because it has no
singular axis and produces conforming SHARED seams across all 6 blocks
— the right choice for production near-body meshing. The two other
topologies cover specific use cases:

1. The **O-grid** is the simpler, single-block sphere covering. Pick
   it (`--topology ogrid`) when you need a single periodic-j block —
   typically as a Chimera background, or as input to a downstream
   tool that doesn't grok multi-block. Note the polar singular axis:
   `i = 0` and `i = ni - 1` are degenerate point rows. The
   cubed-sphere topology is structurally independent of the O-grid
   (it does not wrap an O-grid in caps); the two share only the
   underlying spacing-law primitives.
2. The **Yin-Yang** is in the Python API for users who want a
   pole-singularity-free single-block alternative; it requires a
   downstream Chimera-capable solver to handle the overlap.

For new geometries the workflow is:

1. Confirm the body is **star-shaped** from the chosen anchor (the
   default centroid). Concave geometries with re-entrant features must
   be decomposed first.
2. Pick `theta_cap` and `(ni, nj)` so they land in the cap-equator
   compatibility window (see
   [k=0 mesh construction](/algorithm/k0-construction#compatibility-window-and-fallback)).
3. Run `shore mesh body.stl --ni N --nj M --nk K --theta-cap-deg ANGLE`.

## See also

- [Pipeline overview](/concepts/pipeline) — how a topology fits into
  the overall STL → mesh flow
- [Algorithm overview](/algorithm/) — the building blocks
- [Chimera context](/concepts/chimera) — why these topologies are
  designed for overset workflows
