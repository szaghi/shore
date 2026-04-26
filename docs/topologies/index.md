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
| [O-grid](o-grid) | conceptual primitive | 1 (single) | Pedagogy + downstream tools that expect a single periodic-j block |
| [Caps (6-block)](caps) | **production default** | 6 (4 sub-blocks + 2 caps) | The `shore mesh` CLI; clean polar coverage with no singular axis |
| [Yin-Yang](yin-yang) | Python-API only | 2 (overlapping) | Specialised; not wired to the CLI |

Future topologies that fit naturally into this framework:

- **C-grid** for airfoils and wings with sharp trailing edges
- **H-grid** for bluff bodies and internal flows
- **Multi-body Chimera assemblies** for compound geometries

Each new topology gets one page here; the algorithm pages stay shared.

## How to choose

The CLI builds the **caps topology** unconditionally. The other
topologies are documented because:

1. The **O-grid** is the conceptual primitive that the caps topology is
   built from — understanding it explains the equator and the j-periodic
   wrap.
2. The **Yin-Yang** is still in the Python API for users who want a
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
