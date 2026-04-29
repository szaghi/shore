# Algorithm

SHORE turns a triangulated body surface into a structured volume mesh
through a fixed pipeline of well-defined stages. This section describes
each stage in algorithmic detail, independent of which **topology** the
stages are arranged into. For the topology view, see
[Topologies](/topologies/).

![Pipeline overview](/figures/pipeline.png)

The four algorithm pages, in order of execution:

1. [**Surface projection**](projection) — ray-cast the body STL from a
   single anchor point to obtain a structured `(ni, nj, 3)` array of
   k=0 surface nodes. Pure geometric kernel; topology-agnostic.
2. [**k=0 mesh construction**](k0-construction) — combine the projected
   nodes into the 6-block cubed-sphere topology: gnomonic flat-square cap k=0
   seed + great-circle meridian equator with C1 seam-step pinning.
3. [**Hyperbolic marching**](hyperbolic-marching) — extrude all 6 blocks
   simultaneously layer by layer along outward normals, with Laplacian
   smoothing and a per-layer Jacobian guard.
4. [**Seam-aware normals**](seam-aware-normals) — the centred-difference
   trick that keeps spacing continuous across multi-block seams during
   marching.

A fifth, transverse page documents the
[**spacing laws**](spacing-laws) (`uniform`, `tanh`, `tanh2`, `geometric`)
that control node distribution along i, j, and k.

Two auxiliary algorithms, branched off the main pipeline:

- [**Marker boxes**](boxes) — voxel-fill + greedy-merge generation of
  `cc.par` hole-cutting boxes from a body STL.
- [**C-grid construction**](c-grid) — TE detection, span slicing, and
  C-curve build that feed into the
  [C-grid topology](/topologies/c-grid).

::: tip Read in pipeline order
The pages are designed to be read top-to-bottom in the order above.
Each page assumes the previous ones as background and links to them
where needed.
:::
