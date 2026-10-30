---
title: shore.topology.sphere_cap
outline: [2, 3]
---

# `shore.topology.sphere_cap`

Polar cap block generation for sphere O-grids.

Algorithm overview
------------------
A single-block sphere O-grid (ni, nj, nk, 3) has a periodic j-direction and
open i-boundaries near the poles.  This module closes both polar holes by:

1. Splitting the equatorial block into 4 conforming sub-blocks along j.
2. Building a gnomonic seed face for each polar cap from the 4 corner nodes
   that meet at each pole.
3. Marching all 6 blocks simultaneously, layer by layer, with the equatorial
   sub-blocks driving the shared faces of the cap blocks at every layer.

The result is a list of 6 arrays in order:

```
[sub0, sub1, sub2, sub3, cap_north, cap_south]

```
where `sub_q` has shape `(ni, nj//4+1, nk, 3)` and each cap has shape
`(nc, nc, nk, 3)` with `nc = nj//4 + 1`.

Block adjacency (shared faces)
-------------------------------
- sub0 ↔ sub1  (at j = nj//4)
- sub1 ↔ sub2  (at j = nj//2)
- sub2 ↔ sub3  (at j = 3\*nj//4)
- sub3 ↔ sub0  (at j = 0 / nj, periodic seam)
- sub0, sub1, sub2, sub3 ↔ cap_north  (at i = ni-1)
- sub0, sub1, sub2, sub3 ↔ cap_south  (at i = 0)
