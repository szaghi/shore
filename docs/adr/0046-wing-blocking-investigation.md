# 0046 — A wing on a plane with the blocking tools (investigation)

- **Status:** Investigation (no code change; options for a decision)
- **Date:** 2026-10-29

## Context

Phase 7's acceptance is a wing–body blocking built interactively and
replayed from its journal after a CAD edit. ADR 0042–0045 provide the
blocking commands: create, split, O-grid, delete, associate, move vertex,
premesh. Before choosing the next slice, this checks how far they get on
the simplest wing case, with no code change.

**Run unattended (night of 2026-10-29).** The scripts are in the session
scratchpad `wing46/`.

## Case

- **Geometry:** a NACA 0012 wing with blunt TE, chord 1, span 2 along z,
  root on the plane z = 0. It is an STL with end caps (`extrude_section`).
- **Blocking:** box [-3, 4] × [-3, 3] × [0, 5] with 41 × 21 × 21 nodes, cut
  at x = -0.2 and 1.2, y = ±0.3, z = 2.3. The block around the wing is the
  centre block.

**Metric.** Wall-line angle: the angle between a wall node's first grid
segment and the STL normal at its closest point. The root row (z = 0) is
left out; there the nearest triangle is the root cap, which gives an
artifact of 90°.

## Results (T1)

| Route | Blocks | Min scaled Jacobian | Wall angle: median / p95 / max |
|---|---|---|---|
| **A.** O-grid the centre block (open at the root), associate the core's faces, delete the core | 22 | **−0.10** (LE, TE, tip layers) | — |
| **B.** Associate the centre's faces, delete it, O-grid the hole (factor 0.5) | 22 | 0.44 | 35° / 52° / 55° |
| B with 57 × 31 × 31 nodes | 22 | 0.39 | 36° / 50° / 55° |
| B with factor 0.3 | 22 | −0.09 (corner layers) | — |
| B with the hole's corners moved onto the wing (LE shoulders, TE corners) | 22 | ~~−1 (collapsed)~~ 0.26 / 0.11 (shoulders at x = 0.04 / 0.10), see the correction | 37° / 55° / 71°; 39° / 56° / 84° |

## Findings

1. **Closest-point projection from inside collapses.**
   - Route A's core lies inside the wing's envelope (core y half-width
     0.15 against thickness 0.06 at most). Its faces project onto the
     nearest upper or lower surface, never around the LE or TE.
   - This is ADR 0043's failure of the centre-block route when the box is
     not larger than the body, made worse by a thin body.
2. **Projection from outside is valid but skews the wall lines.**
   - In route B the hole's faces lie outside the wing and project
     smoothly.
   - But a flat face's closest points crowd at the LE and TE (the face is
     longer than the chord). The lines from wall to layer then follow the
     box, not the normal: median 35°.
   - Refining does not help; it is the parametrization, not the
     resolution.
3. **~~Placing vertices on the wing needs edges that do not cut through it.~~** Wrong; see the correction below.
   - A straight edge between an upper and a lower LE shoulder crosses the
     wing's interior. Projecting it folds onto one side, so the LE face
     collapses.
   - ICEM avoids this with curve associations (section curves) and
     smoothing. Edges projected by closest point from a straight chord
     cannot wrap a nose.
4. **A sharp-ish TE is the weakest spot** in every route. Faces facing the
   TE project to its edge line.

So the blocking commands make valid wing topologies (route B), but premesh
by TFI and closest-point projection gives wall quality far below SHORE's
marched C-grids, whose lines leave the wall nearly normal.

## Correction (2026-10-29, ADR 0047)

- **Finding 3 came from a script bug.** The moves were applied after the
  O-grid, to the hole faces looked up by their old `(block, face)` names.
  After the O-grid those names give the set blocks' inner copies, which
  carry no association. The moved vertices were therefore the copies, not
  the wall.
- **Moved correctly** (the vertices of the associated faces), the corners
  land on the wing and premesh stays valid, though worse than without the
  moves:
  - minimum scaled Jacobian 0.26 with the LE shoulders at x = 0.04, 0.11
    at x = 0.10;
  - wall spacing ratio 3.6 / 9.7;
  - median wall angle 37° / 39°.
- **The other numbers stand.** Re-measured on a wing STL without the root
  cap (the root lies on the plane), route B gives median 35.9°, p95
  53.3°: the same.
- **Option 2's premise does not hold:** edges between vertices on the wing
  do not fold. ADR 0047 tried that option anyway and rejects it.

## Options (for a decision)

| Option | What | Reuses | Cost |
|---|---|---|---|
| **1. Smooth after premesh** | Wall faces smoothed on the geometry (Laplacian or elliptic with re-projection), then the blocks' interiors by Winslow with orthogonality control at walls (Thomas–Middlecoff) | `collar.py`'s TFI + Winslow interior; ADR 0007's elliptic smoother (CAD (u, v) only, not STL) | Medium: a smoother for STL walls, orthogonality forcing; spec in quality metrics |
| **2. Edges that wrap the body** | Initialise an associated edge by projecting the layer's outer edge (from outside) instead of the straight chord; vertices placed on the body then work | the O-grid's correspondence of layer and hole | Small; expected (untested, T0) to fix finding 3, not 2 |
| **3. Hybrid: march the wall layer** | Use the blocking for the far field and SHORE's hyperbolic marching for the wall layer: premesh the hole's faces as a surface grid, march it, join it to the layer blocks | `march`, `march_surface` (ADR 0005, 0008) | Medium to large; best wall quality (SHORE's strength) |
| **4. Curve associations on STL** | Section, LE and TE curves extracted from the STL (feature edges, plane cuts) for edges to follow | `surface/cgrid.py`'s TE detection and section slicing | Medium; ICEM-like, but CAD covers it already |

**Recommendation:** option 2 first (small, and it makes vertex placement
usable), then option 1 (smoothing). Together they cover ICEM's workflow on
STL and CAD. Option 3 is the stronger long-term route for viscous walls,
and deserves its own ADR.

## Verification

None: an investigation, with T1 numbers from the scripts in `wing46/`
(`hole.py`, `variants.py`, `measure.py`).
