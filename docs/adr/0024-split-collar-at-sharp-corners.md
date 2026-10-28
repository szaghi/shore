# 0024 — The split collar at sharp corners

- **Status:** Accepted
- **Date:** 2026-10-27

## Context

The single collar block (ADR 0010) folds at the corners of its junction
loop (ADR 0021, 0022). The user asked for a split, multi-block collar,
covering four failures:
- sharp corners;
- concave corners;
- thick collars at convex corners;
- the blunt trailing edge.

The user also asked for prototypes of both conforming and overlapping
block layouts. The acceptance criteria are:
- no folded cells;
- consistent seams;
- smooth loops unchanged;
- quality bounds.

## Investigation (2026-10-27, prototypes in the session scratchpad; T1)

**Where the single block folds.** n_loop 96, surface steps
`0.005 * 1.12**j`, wall-normal `5e-5 * 1.2**k`; thickness in units of the
first surface step *h* (0.005). Each entry is folded cells / folded cells
in the wall layer.

| Junction | 1.9 × h | 11.8 × h | 11.8 × h, curvature 1 |
|---|---|---|---|
| rectangle prism, convex 90° | 0 | 20 / 0 | 52 / 0 |
| L-shape, one concave 90° | 4 / 0 | 126 / 0 | 291 / 0 |
| blunt NACA 0012 STL wing | 7 / 0 | 27 / 0 | 71 / 0 |
| sharp NACA 0012 STL wing (164°) | 12 / **2** | 27 / **2** | 7534 / **225** |
| biconvex lens (two 157° tips) | 40 / **4** | 60 / **4** | 2885 / **55** |

**Two different failures.**
1. **Sharp corners fold the wall.** The surface march on A cannot fan
   through a 164° turn. This one needs a different topology.
2. **Every other fold is in the volume.** The wall layer is valid. The
   cause is the construction: the volume is built per loop station, as 2-D
   (i, k) sheets that know nothing of their neighbours. It is not caused by
   the topology:
   - **The folds grow with refinement.** At 4.7 × h the rectangle has 0, 0
     and 8 folds at n_loop 96, 192 and 384; the L-shape has 10, 271 and
     1014. They scale with the thickness over the loop spacing.
   - **A split cannot take a convex corner out of a block.** A corner block
     would have a wall collapsed onto the edge, which Xnavis rejects
     (ADR 0022).
   - **Three volume fixes failed:**
     - corner-aware normals at corner stations: blunt wing 7 → 1, no effect
       elsewhere;
     - a full 3-D Winslow solve with the outer boundary fixed: rectangle at
       11.8 × h 20 → 12, about 300 s;
     - redistributing the outer surface in j, row by row: worse, because it
       breaks the alignment in i.

   These go to ADR 0025 (user decision of 2026-10-27).

**Overlapping blocks are ruled out by argument, not measured** (user
decision of 2026-10-27). A block extended past a corner has no wall of B
beyond the edge. SHORE's whole-face boundary conditions cannot make a face
partly wall and partly free. Overlap belongs between the collar and the
body grids, not inside the collar.

**The C-collar.** The design takes the C-grid's wake cut to the collar:
- **Cut.** The loop is cut at the sharp corner. The piece is extended at
  both ends by a wake line on A, along the bisector of the corner's two
  loop directions.
- **Sheet.** Over the wake, the wall is a virtual sheet: the corner's
  B-line translated along the wake. The two sides of the sheet are joined
  by a SHARED seam.
- **Small turns.** The front on A now turns by half the wedge angle (8° on
  the wing), not by 164°.

Measured results:

| Sharp junction, 1.9 × h | single block (folded) | C-collar (min scaled Jacobian) |
|---|---|---|
| NACA 0012 STL wing, uniform | 12 | 0.155 |
| NACA 0012 STL wing, curvature 1 | 5500 | 0.228 |
| lens, uniform | 40 | 0.462 |
| lens, curvature 1 | 971 | 0.447 |
| lens, uniform, 4.7 × h | — | 0.460 |
| rounded wedges, turn 120-160°, uniform | 2-24 | 0.439-0.471 |

Two things were needed on the way:
- **B's march between two planes.** The march along B must hold the
  piece's ends on the edge. Free ends drift off the edge. Mirror ghosts
  (ADR 0012) work, but the lens's pieces end on two different sheets, so
  the marcher takes one plane per end (private `_Ends`).
- **Sharp creases in the marcher.** At a sharp crease, the central
  difference across the fold is about sin 8° of the spacing. The
  prescribed area then gives an outsized step: the node overshoots and
  drags its neighbours.
  - This is also the reason why the existing *periodic* march of the wing
    flips with curvature: 66 cells, drift 0.98.
  - Lengthening the tangent to the node's arc spacing where it is shorter
    than half of it fixes it: 0 flips, drift 0.011.
  - The blunt trailing edge is not affected: there the difference is 0.7
    of the spacing.

**Threshold.** On rounded wedges, the single block is valid up to a 110°
turn with uniform spacing (min 0.021) and folds from 120°. With curvature 1
it folds from 100°. Sharp corners are those turning more than 105°.

**Limit.** With curvature 1, the tip spacing is about 0.001, against a
collar of 0.0093 (1.9 × h). The C-collar is valid at that thickness. At
4.7 × h it folds in the volume (164-187 cells, wall valid). On wedges
with a 0.1 rounded back, curvature already folds it at 1.9 × h (32-516),
next to the corner station. All of this is the volume construction of
ADR 0025.

## Decision

User decisions of 2026-10-27:
- two ADRs, this one first;
- a new `split_collar`; `build_collar` unchanged;
- acceptance: no folds, seams, unchanged, quality bounds.

**`shore.collar.split_collar(body_a, body_b, seed, *, n_loop, steps_a,
steps_b, steps_k, label, params, loop_spacing, loop_curvature,
wake_length) -> list[HexBlock]`.**

- **Sharp corners.** They are the tracer's corners (ADR 0021) that turn by
  more than `CollarParams.sharp_deg` (105°). The private
  `intersect._loop` returns the corner node indices.
- **Without sharp corners,** the result is `build_collar`'s block, bit for
  bit.
- **Wake.** At each sharp corner a wake line on A runs along the bisector
  `t_in - t_out`. Its cells start at the corner's mean loop spacing and
  grow by at most `CollarParams.wake_ratio` (1.2) to `wake_length`
  (default `sum(steps_a)`). The points are projected onto A.
- **Sheet.** The plane through B's edge (`n1 × n2` of B's normals on both
  sides) along the wake. B is marched from each piece with its ends held
  on the two sheets (`_march_between`). The corner's B-line is the mean of
  the piece ends that meet there. The sheet is that B-line translated
  along the wake.
- **Volume.** Each piece (wake, piece, wake) is extruded as one open block
  (ADR 0010's construction, j not periodic). It is then split at the
  junction row and at the two corners, into six blocks
  `f"{label}_{p}_{part}_{side}"`.
- **Boundary conditions.** The sheet blocks' `k_lo` faces are SHARED
  across the cut, with AxisMap `((i, False), (j, True))`. Every other
  `k_lo` is WALL. The seams inside a piece are SHARED (`split_blocks`).
- **Marcher.** At a crease where `|t|` is less than
  `MarchParams.sharp_crease` (0.5) of the arc spacing, `t` is lengthened
  to the arc spacing. Other marches are unchanged.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Overlapping (chimera) blocks per segment | No wall of B past a corner; whole-face BCs (argument). |
| Corner fan in one block | Collapsed wall cells (ADR 0022). |
| Corners from the loop's polyline turns | A smooth LE turns 82° between nodes at n = 96 (blunt wing); the tracer knows the real corners. |
| Free ends of B's march at the sharp edge | The end leaves the edge: 5 flips uniform, worse with curvature. |
| Linear-extrapolation ghosts at the ends | 42-79 flips with curvature. |
| Corner-aware normals at the cut corners (volume) | No effect on the curvature folds. |
| One ADR for sharp, concave, thick and blunt corners | The last three are volume construction; they share no mechanism with the cut (user decision). |

## Consequences

- **Sharp junction corners now have valid collars:** sharp trailing edges,
  lens tips and wedges, as thin collars at sharp corners. The result is 6
  blocks per piece, written to `cc.par` like any SHARED topology.
- **Curvature clustering at sharp corners works with thin collars only.**
  The volume folds as the thickness grows past about 10 × the corner
  spacing. That limit belongs to ADR 0025.
- **The surface marcher no longer flips at sharp creases.** Single collars
  at sharp corners change: they folded before, and the golden collars and
  marches are unchanged.
- **Not covered:** half-collars with sharp corners, sharp corners on A's
  edges (they raise), and the project/CLI layer, which is a later step.

## Verification

### Acceptance spec (frozen 2026-10-27, before implementation)

**Spec.** `tests/spec/test_split_collar.py` runs in the main CI job. The
golden file `tests/spec/data/collar_default.npz` was written by
`build_collar` before this ADR and is hashed in `SHA256SUMS`.

**Verdict.** It must report 33 passed, 0 skipped, 0 xfailed. At freeze it
reports 5 passed and 28 xfailed. The 5 are regression guards: the golden
collars and the uniform marches at sharp edges. Validated with the
prototype behind the interface (a copy of the spec without `@PENDING`): 33
passed.

| Case | Required |
|---|---|
| S1 (× 9) | min scaled Jacobian ≥ 0.1 in every block: wing and lens, uniform and curvature 1 (0.155-0.462); lens at 4.7 × h (0.460); wedges turning 120-160° (≥ 0.439) |
| S2 (× 5) | labels, face BCs and maps; every SHARED pair equal node for node; `.adjacency.json` + `cc.par` written |
| S3 (× 2) | WALL faces on their body (1e-9 × diagonal, STL 1e-6); the junction rows on both bodies |
| S4 (× 4) | wake lines on A, cell ratio in [1, 1.2], length `wake_length` (default and 0.1); sheets a translation of the corner's B-line |
| S5 (× 4) | periodic marches of B at sharp edges: 0 flips, drift ≤ 0.02, uniform (guards) and curvature 1 (0.011) |
| S6 (× 5) | golden `build_collar` blocks within 1e-12 (rectangle with curvature, blunt STL wing, half ellipse); `split_collar` gives them without sharp corners |
| S7 (× 4) | `ParameterError` for a bad `wake_length` and a sharp corner on A's edge |

**Out of scope:** thicker collars at corners, concave corners and the blunt
trailing edge's volume (ADR 0025); half-collars with sharp corners; sharp
corners on A's edges; project and CLI.

### Review (2026-10-27, recommendations taken while the user was away)

The user asked for the night's work to go ahead on the recommendations:
1. The threshold stays at 105° (`CollarParams.sharp_deg`).
2. The wake length defaults to the collar's width on A (`sum(steps_a)`).
3. Next: the project and CLI layers for `split_collar`, then ADR 0025.

## Acceptance (2026-10-27): split-collar spec met

**Kernel.** As decided. `intersect._loop` returns the corner nodes of
`_resample` (`with_fixed`). `march._march_between` runs the march with an
`_Ends` pair. `march_surface(plane=)` keeps its path: `_Mirror` answers
`end` and `matrices` itself, and the march body is shared (`_march`).
`build_collar`'s closed and half collars go through `_closed_collar` and
`_wall`, bit-identical to before.

**Verdict for `tests/spec/test_split_collar.py`: 33 passed, 0 skipped,
0 xfailed** (T2; `scripts/check_spec.py` reports all files ok). All 28
`@PENDING` lines are removed.

**Unit tests:**
- `tests/test_surface_march.py`:
  - on a front folded round a 16° edge, every node marches the step
    within 1e-4; without the fix the edge node marches 3 % too far;
  - the same plane at both ends is bit-identical to `plane=`;
  - each end stays on its own plane.
- `tests/test_collar.py`:
  - wake grading;
  - `_loop`'s corner nodes on the rectangle;
  - a bad `wake_length`.
