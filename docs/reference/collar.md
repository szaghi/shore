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
| `--name-a`, `--name-b TEXT` | file stems | Geometry names in the project. |
| `-o, --output PATH` | | Output stem for the `.geo` file. |
| `--adjacency-json PATH` | | Also write the `.adjacency.json` sidecar. |
| `--project DIR` | | Add the collar to the project in `DIR` (created if missing) and save it. |
| `--component TEXT` | `collar` | Component of `--project` to build. |

At least one of `-o` and `--project` is required.

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
