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
itself (ADR 0021). The collar block itself still folds at such corners:
a collar that wraps a trailing edge is a later step.

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
