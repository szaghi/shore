# 0023 — Surface marching at creases

- **Status:** Proposed
- **Date:** 2026-10-27

## Context

ADR 0022 left the wing-plate collar folding, and blamed ADR 0008's surface
marcher at a clustered leading edge. The user asked for that to be fixed,
with a valid wing-plate collar as the goal, and with smooth marches
unchanged bit for bit.

## Investigation (2026-10-27, prototypes in the session scratchpad; T1)

**The premise was wrong.**
- **Not the leading-edge radius:** ellipses with NACA 0012's LE radius
  march without flips, analytic and STL, uniform and curvature-clustered,
  with steps up to 9.5 × the front spacing.
- **The inside-out wing:** `shore.shapes.extrude_section` wrote inside-out
  STLs, with inward side faces and outward caps. The loop then ran
  backwards, and the plate march went into the wing at the LE. This was
  fixed separately: outward, consistent winding.
- **Two tracer bugs:** the same work exposed the tracer bugs on fine
  sharp trailing edges, and `load_stl`'s missing `networkx` dependency.
  Both were fixed separately.

**The real defect, on a correctly wound wing.** Wing-side marches flip
cells near the blunt TE, and the nodes drift chordwise by up to 0.13 chord
while marching up a prismatic wing, where the march should be a pure
translation. The exact analytic prism of the section drifts by the same
amount, so it's the scheme, not the STL.

The mechanism, at a corner of the front on an edge of the surface (the TE
crease):
- the node's normal is one face's normal and leans along the
  central-difference tangent `t`;
- the scheme assumes `n ⊥ t`, so the area row divides by `|n × t| < |t|`,
  and the marching distances vary from node to node;
- the orthogonality row `t·dr + dr0·D(dr) = 0` turns those varying
  distances into tangential motion.

Sub-stepping (×4, ×16) and lower dissipation don't remove the drift:
- with uniform spacing it stays at 0.06-0.15;
- with curvature it reaches up to 0.94.

**The fix.** Where `|n · t̂| > 0.1` (`MarchParams.crease_tolerance`), the
marching frame uses the normal made orthogonal to the tangent.

| Blunt NACA 0012 × plate (n = 96) | wall flips before → after | B drift before → after |
|---|---|---|
| analytic prism, uniform | 11 → 0 | 0.130 → 0.060 |
| analytic prism, curvature 1 | 0 → 0 | 0.011 → 0.010 |
| STL, uniform | 10 → 0 | 0.131 → 0.060 |
| STL, curvature 1 | 73 → 0 | 0.217 → 0.011 |

- **Smooth fronts stay untouched.** The threshold is never crossed on the
  13 golden marches: ADR 0008's seven cases, the elliptic junction of
  ADR 0020 (loop and half curve, both bodies) and a faceted STL sphere ×
  cylinder. All are bit-identical.
- **The blunt-wing collar:**
  - before (correctly wound wing): 155-180 folded cells, uniform;
  - after: **7** (uniform) and **6** (curvature) at 1.9 × h, **none** in
    the wall layer.

  The remaining folds sit in the volume, at the junction row next to the TE
  corners. There the base (0.0025) is much thinner than the collar, and the
  two 90° corners act like a sharp edge: ADR 0022's limit, the split
  collar's job.
- **Sharp TE:** the plate march folds at the 164° turn (1 cell uniform,
  35-50 with curvature). Out of scope, as in ADR 0022.

## Decision

User decisions of 2026-10-27: fix first, then this ADR; smooth marches
bit-identical when the fix is not triggered.

- `MarchParams.crease_tolerance = 0.1`: in `_layer`, nodes with
  `|n · t̂|` above it use `n - (n · t̂) t̂`, normalised, for the whole
  marching frame.
- The wing-plate collar of a blunt TE is valid at the wall. Its volume
  next to the TE stays with the split collar, together with sharp
  trailing edges.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Sub-stepping | Drift stays (0.12-0.86); flips reduced, not removed. |
| Lower dissipation / area smoothing | Drift up to 0.94 with curvature. |
| Always orthogonalise | Same result here, but changes every march bit-wise (saved projects). |
| Normal from the front's own geometry (cross of t and the marching direction) | Needs the marching direction, which the frame defines: circular. |

## Consequences

- Surface marches from corner loops, the collar walls of struts and blunt
  trailing edges, no longer drift or flip at the creases.
- The wing-plate collar's remaining folds are volume folds at the TE
  corners. The split collar ADR is next for those and for sharp TEs.

## Verification

### Acceptance spec (frozen 2026-10-27, before implementation)

**Spec.** `tests/spec/test_march_creases.py` runs in the main CI job. The
golden file `tests/spec/data/march_default.npz` was written by the marcher
before this ADR and is hashed in `SHA256SUMS`.

**Verdict.** It must report 22 passed, 0 skipped, 0 xfailed. At freeze it
reports 15 passed and 7 xfailed. The 15 are regression guards: the golden
marches, and the prism cases with curvature, which already held.
Validated with the prototype behind the interface: 22 passed.

| Case | Required |
|---|---|
| C1 (× 4) | no flipped surface cells on either march: prism and STL, uniform and curvature 1 |
| C2 (× 4) | B's in-section drift ≤ 0.08 uniform (0.060), ≤ 0.02 with curvature (0.011) |
| C3 | the blunt STL wing collar's wall layer has no folded cells |
| C4 (× 13) | golden marches within 1e-12 |

**Out of scope:**
- sharp trailing edges;
- the collar volume at corners;
- marching speed on sliver STL triangles.

### Review questions

1. **The threshold.** Is 0.1 (about 6° of lean) right, given STL facets on
   coarse meshes? It was never crossed on the golden set.
2. **Next step:** the split collar (sharp TE, the volume at blunt TE
   corners, concave corners, thick collars), or a 0.3.0 release first?
