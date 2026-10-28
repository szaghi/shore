# shore collar

Build the collar block around the junction of two bodies: a Chimera volume
block whose wall lies on both bodies along their intersection curve
(ADR 0010). With `--plane`, the half-collar of a half model on that
symmetry plane (ADR 0013).

## Synopsis

```
shore collar [OPTIONS] BODY_A BODY_B
```

## Arguments

| Argument | Description |
|---|---|
| `BODY_A` | Body A (STL or STEP/IGES/BREP); the collar's surface layers `--m-a` march on it. |
| `BODY_B` | Body B, meeting body A along a closed (or, with `--plane`, open) curve. |

## Options

| Option | Default | Description |
|---|---|---|
| `--seed X Y Z` | required | A point near the junction; selects the intersection curve. |
| `--n-loop INT` | required | Nodes along the junction curve. |
| `--ds-surface FLOAT` | required | First surface step away from the junction. |
| `--growth-surface FLOAT` | required | Surface step growth ratio. |
| `--m-a INT`, `--m-b INT` | required | Surface layers on body A and on body B. |
| `--ds-wall FLOAT` | required | First wall-normal step. |
| `--growth-wall FLOAT` | required | Wall-normal growth ratio. |
| `--nk INT` | required | Wall-normal layers. |
| `--label TEXT` | `collar` | Block label. |
| `--plane PX PY PZ NX NY NZ` | | Symmetry plane (point, normal): build the half-collar on the seed's side. Its j-faces are `SYMMETRY` (cc.par code 2). |
| `--loop-law TEXT` | `uniform` | Node law along the junction curve: `uniform` or `tanh2` (`tanh` on a half-collar's open curve too). |
| `--loop-beta FLOAT` | `3.0` | Clustering strength of `tanh` / `tanh2`. |
| `--loop-cluster-at FLOAT` | | Where `tanh2` clusters, as a fraction of the curve in [0, 1] (from the seed round the loop, or from the first plane end). Default: both ends, i.e. at the seed on a loop. |
| `--loop-curvature FLOAT` | `0.0` | Cluster nodes where the curve bends (ADR 0020). Try `1` for rounded sections. |
| `--split` | off | Cut the collar at sharp corners (> 105 degrees) into blocks joined across wake sheets (ADR 0024). Closed loops only (not with `--plane`). |
| `--wake-length FLOAT` | collar width on A | With `--split`: length of the wake lines behind the sharp corners. |
| `--collar-volume TEXT` | `smoothed` | Outer boundary of the volume: `smoothed` (station offsets smoothed along the loop near corners, ADR 0025-0026) or `stations` (per loop station, ADR 0010: grids of earlier versions). |
| `--name-a`, `--name-b TEXT` | file stems | Geometry names in the project. |
| `-o, --output PATH` | | Output stem for the `.geo` file. |
| `--adjacency-json PATH` | | Also write the `.adjacency.json` sidecar. |
| `--project DIR` | | Add the collar to the project in `DIR` (created if missing) and save it. |
| `--component TEXT` | `collar` | Component of `--project` to build. |

At least one of `-o` and `--project` is required.

## Node spacing along the junction curve

By default the nodes are equally spaced in arc length. Around a thick,
rounded section (a strut, a wing with a rounded trailing edge) that leaves
too few nodes at the leading and trailing ends, and the collar degrades or
folds there. `--loop-curvature A` clusters nodes where the curve bends.
Nodes equidistribute `1 + A * kappa * L / (2 pi)`, where `kappa` is the
curvature and `L` the curve length, and adjacent cells stay within 1.2x of
each other. A circle stays uniform at any weight. On a plate pierced by a
1 x 0.2 elliptic post with 128 nodes:

- `A = 0` (uniform): the ends turn 43 degrees per cell, and the collar's
  minimum scaled Jacobian is 0.008;
- `A = 1`: the ends turn 5.5 degrees per cell, and the minimum scaled
  Jacobian is 0.23.

`--loop-law` applies a `Spacing1D` law on top. With `A = 0` it holds
exactly by arc length; with `A > 0` it acts on the curvature-distributed
coordinate. On a half-collar the curvature is measured together with its
mirror image, so the ends on the plane are clustered too.

Junction curves with corners (a sharp or blunt trailing edge through a
plate, a prismatic strut) are traced: every corner, where the curve turns
by 30 degrees or more, becomes a node, and each segment between corners
(and the seed) gets cells in proportion to its length and the law by
itself (ADR 0021). With `--loop-curvature`, one graded spacing runs over the
whole curve, so adjacent cells stay within 1.2x across the corners too
(ADR 0022). The collar block is valid at convex corners of up to 90 degrees
only while it stays thin: up to about 4.7 junction spacings with a uniform
loop, about 1.9 with `--loop-curvature 1`. Thicker collars and concave
corners folded with the per-station outer boundary of earlier versions
(`--collar-volume stations`). The default `smoothed` boundary (ADR 0025,
0026) builds each station's outer line on
its 3-D offset line and, near corners and bends tighter than the collar is
thick, smooths the station offsets along the loop over half the collar's
thickness, so neighbouring stations agree: it fans round
convex corners and stops crossing at concave ones. On the rectangular strut
at 11.8 junction spacings it turns 20 folded cells (52 with curvature) into
none; on 49 test junctions the folded cells outside a folded wall go from
2367 to 72. Concave corners at about 5 junction spacings and thicker, and a
blunt trailing edge's base, still fold a few cells.

Sharp corners, those turning by more than 105 degrees (a sharp trailing
edge, the tips of a lens), fold even the collar's wall: the surface grid
on the plate cannot turn round them. `--split` (Python:
`shore.collar.split_collar`, `Project.build_collar(split=True)`) cuts the
loop there (ADR 0024). A wake line runs
on the plate behind each sharp corner, along its bisector, for
`wake_length` (default: the collar's width on the plate), with cells
growing by at most 1.2 from the corner spacing. Over the wake, the wall is
a virtual sheet: the edge's grid line on the other body, translated along
the wake. The two sides of the sheet are joined by a shared seam, as in a
C-grid. Each piece of the loop between sharp corners becomes six blocks
(the wake, the piece and the other wake, below and above the junction
row). On a NACA 0012 wing with a sharp trailing edge, the single block
folds 12 cells (5500 with `--loop-curvature 1`); the split collar's
smallest scaled Jacobian is 0.155 (0.228). With curvature clustering at the
corner the split collar must stay thin: about 10 corner spacings. With
`-o` each block goes to `<stem>_<label>.geo`, labels
`<label>_<piece>_<wake_lo|body|wake_hi>_<b|a>`; without sharp corners
`--split` gives the single block.

## Projects

With `--project`, the project is loaded (its journal replays), the collar
is built in `--component`, and the project is saved back. A geometry name
already loaded from a file with the same content is reused, so earlier
components built on it survive. A file whose stem matches a geometry
loaded from a *different* file is an error: name it with `--name-a` or
`--name-b`.

```bash
shore mesh fuselage.stl --topology ogrid -o out/fus --project proj
shore collar fuselage.stl wing.stl --seed 1.2 0.0 0.3 --n-loop 96 \
    --ds-surface 0.005 --growth-surface 1.12 --m-a 20 --m-b 20 \
    --ds-wall 1e-5 --growth-wall 1.15 --nk 40 --project proj
shore xall proj --out run --options xall.toml
```
