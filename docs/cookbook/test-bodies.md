# Make a test body

You need a body to mesh before you have the real CAD: an analytic shape, an
airfoil section or a wing built from one. SHORE generates all three, as STL
files from the CLI or as `trimesh.Trimesh` / `numpy` arrays from Python.

The Python code on this page is
[`examples/cookbook/bodies.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/bodies.py);
it runs in under a second.

## Write an analytic body as STL

::: code-group
```bash [CLI]
shore body sphere --radius 1.0 --subdivisions 3 -o sphere.stl
shore body capsule --radius 0.5 --height 2.0 --sections 24 -o capsule.stl
shore body box --extents 2 1 0.5 -o box.stl
```
<<< @/../examples/cookbook/bodies.py#analytic{python} [Python]
:::

```text
torus: 2048 triangles, watertight=True
```

- Every body is centred at the origin; the cylinder, capsule and torus have
  their axis along z.
- `sphere` is an icosphere: each subdivision level multiplies the triangle
  count by 4 (level 3 = 1280 triangles). Its radius is exact only at the
  vertices, so a fine level keeps the faceting error small against your
  first cell.
- `capsule` (`height` = length of the cylindrical middle) and `sphere` are
  star-shaped from the centroid, so they suit the `cubed_sphere`, `ogrid`
  and `oh` topologies. The torus is not star-shaped: use it to test
  geometry and association code, not the ray-cast projection.

## Generate a NACA section

::: code-group
```bash [CLI]
shore profile naca 0012 -o naca0012.dat
shore profile naca 2412 --open-te -o naca2412_blunt.dat
```
<<< @/../examples/cookbook/bodies.py#naca{python} [Python]
:::

```text
NACA 0012: 159 points; NACA 2412 blunt TE gap 2.52e-03 chord
```

- 4-digit (`0012`, `2412`) and 5-digit (`23012`, reflexed `23112`)
  designations are supported; anything else raises `ParameterError`.
- `n_pts` (`--n-pts`, default 80) is the number of cosine-spaced chord
  samples per side; the section has `2 * n_pts - 1` points, ordered upper
  TE → LE → lower TE (Selig order).
- The default closes the trailing edge exactly (sharp TE).
  `closed_te=False` (`--open-te`) keeps the standard polynomial's residual
  thickness, a blunt TE about 0.25 % chord thick, which the C-grid
  topology detects and meshes with its 4-block blunt variant.
- The `.dat` file is plain `x y` rows with an optional name line, the UIUC
  airfoil database format; `read_section` reads such files.

## Extrude a wing with taper, sweep and twist

::: code-group
```bash [CLI]
shore body extrude --section naca0012.dat --span 1.0 --n-stations 8 -o wing_straight.stl
shore body extrude --section naca0012.dat --span 2.0 --n-stations 12 \
    --taper 0.6 --sweep-deg 25 --twist-deg -3 --span-axis y --end-caps -o wing_swept.stl
```
<<< @/../examples/cookbook/bodies.py#wing{python} [Python]
:::

```text
swept wing: extents [1.532 2.    0.12 ], watertight=True
```

| Parameter | Default | Meaning |
|---|---|---|
| `span` | 1.0 | total length along `span_axis` |
| `n_stations` | 8 | span stations including root and tip (>= 2) |
| `taper` | 1.0 | tip-to-root chord ratio |
| `sweep_deg` | 0.0 | leading-edge sweep, positive = sweepback |
| `twist_deg` | 0.0 | tip twist, linear from 0 at the root, positive = nose-up |
| `span_axis` | `z` | `x`, `y` or `z` |
| `end_caps` | `False` | centroid-fan caps at both ends, watertight |

- Leave the ends open (the default) for the `cgrid` and `ch` topologies:
  they cut the wing into span stations and need the open ends to slice
  cleanly. Close them (`end_caps=True`) for anything that needs a
  watertight body: `shore boxes-from-stl` (voxel fill) refuses an open
  STL, and loading one warns.
- The section is stacked as given: the stations are straight copies, so
  the spanwise faceting is linear between stations.

## See also

- [Airfoil profiles](/concepts/airfoil-profiles): the NACA equations and the
  sharp/blunt TE choice.
- [C-grid topology](/concepts/topologies/c-grid): meshing the extruded wing.
- CLI: [`shore body`](/reference/body), [`shore profile`](/reference/profile).
- API: [`shore.shapes`](/reference/api/shore.shapes),
  [`shore.profiles`](/reference/api/shore.profiles).
