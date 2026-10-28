# 0032 — STL projection: each query point computed once

- **Status:** Accepted
- **Date:** 2026-10-28

## Context

After ADR 0031, tracing a junction loop on the sliver-STL wing still took
5-8 s. Nearly all of it is trimesh's `closest_point`: on sliver triangles
the broad phase returns most faces for every point. User decision of
2026-10-28: make it faster, with every result bit-identical.

## Investigation (2026-10-28; T1)

**The answer cannot change, so the candidate sets cannot change either.**
trimesh picks the closest of the candidate faces and breaks ties by
candidate order. Ties are common on crease edges, where they decide the
normal. A tighter candidate set could therefore change normals at creases.

**But most queries are repeats.** On the 1000-point blunt wing, 61 % of the
7898 queried points are exact repeats (4832): polish re-projects nodes that
have already converged, and the tracer revisits points.

**Each point's answer is independent of its batch.** For each query point
trimesh computes candidates, distances and the tie correction on that
point's rows alone. So remembering each point's answer, keyed by its exact
bytes, returns the same arrays bit for bit.

| Loop (NACA 0012, 2 stations) | before | after | result |
|---|---|---|---|
| 1000-point sections, blunt TE | 7.7 s | 4.0 s | bit-identical |
| 1000-point sections, sharp TE | 5.4 s | 3.5 s | bit-identical |
| 200-point sections | 1.3 s | 0.7 s | bit-identical |

The full suite drops from 356 s to 336 s.

## Decision

- `StlGeometry.project` asks trimesh only for the points it has not seen
  (`_closest_points`) and remembers each answer (closest point, distance,
  triangle) under the point's bytes.
- The cache holds at most `StlGeometry._CACHE_SIZE` (200 000) points, about
  20 MB, and is cleared when full.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Tighter candidate faces (a better broad phase) | Can change the tie-break at crease edges, so it would not be bit-identical. |
| An LRU cache | More bookkeeping for no measured gain; clearing a full cache keeps the bound. |

## Consequences

- **Faster:** junction loops and marches on STL bodies, about 2 × on
  sliver STLs.
- **Unchanged:** every result, bit for bit.
- **Memory:** up to about 20 MB per STL geometry.

## Verification

- **Unit tests** (`tests/test_geometry.py::TestProjectCache`):
  - batches, shuffled subsets and repeats give bit-identical points,
    distances and normals, with the cache cleared or not;
  - the cache stays within its bound.
- **Full suite:** passes unchanged.
