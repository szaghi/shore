# shore body

Generate 3D body STL files: analytic primitives (sphere, cylinder,
torus, box, capsule) and arbitrary section extrusions (any 2D Selig
profile → tapered / swept / twisted wing).

`shore body` is a sub-app with six commands:

```
shore body sphere    [OPTIONS]
shore body cylinder  [OPTIONS]
shore body torus     [OPTIONS]
shore body box       [OPTIONS]
shore body capsule   [OPTIONS]
shore body extrude   [OPTIONS]
```

All commands write a single `.stl` file. The primitives wrap
`trimesh.creation` and produce watertight meshes. The `extrude`
subcommand is the bridge between [`shore profile`](profile) and 3D
geometry: feed any Selig `.dat` (yours or from the UIUC database) and
get a wing STL with optional taper, sweep, and twist.

## `shore body sphere`

Icosahedral sphere centred at the origin.

```
shore body sphere -o FILE.stl [--radius FLOAT] [--subdivisions INTEGER]
```

| Option | Default | Description |
|---|---|---|
| `-o` / `--output PATH` | — (required) | Output `.stl` file |
| `--radius FLOAT` | 1.0 | Sphere radius |
| `--subdivisions INTEGER` | 3 | Icosphere subdivision level. 0 = 20 triangles; each level ×4 triangle count. 4 is dense enough for most CFD wraps. |

```bash
shore body sphere --radius 1.5 --subdivisions 4 -o sphere.stl
```

## `shore body cylinder`

Right circular cylinder with circular caps, axis along z.

```
shore body cylinder -o FILE.stl [--radius FLOAT] [--height FLOAT] [--sections INTEGER]
```

| Option | Default | Description |
|---|---|---|
| `--radius FLOAT` | 1.0 | Cylinder radius |
| `--height FLOAT` | 1.0 | Total length along z |
| `--sections INTEGER` | 32 | Segments around the circumference |

```bash
shore body cylinder --radius 0.3 --height 2.0 -o cyl.stl
```

## `shore body torus`

Ring torus, axis of revolution = z. Genus-1 watertight surface.

```
shore body torus -o FILE.stl [--r-major FLOAT] [--r-minor FLOAT] \
                 [--major-sections INTEGER] [--minor-sections INTEGER]
```

| Option | Default | Description |
|---|---|---|
| `--r-major FLOAT` | 1.0 | Major radius (centre of tube to centre of torus) |
| `--r-minor FLOAT` | 0.25 | Minor radius (tube radius). Must be smaller than `r-major`. |
| `--major-sections INTEGER` | 32 | Segments around the major circle |
| `--minor-sections INTEGER` | 32 | Segments around the tube |

```bash
shore body torus --r-major 2.0 --r-minor 0.5 -o torus.stl
```

## `shore body box`

Axis-aligned rectangular box centred at the origin.

```
shore body box -o FILE.stl [--extents FLOAT FLOAT FLOAT]
```

| Option | Default | Description |
|---|---|---|
| `--extents FLOAT FLOAT FLOAT` | `1 1 1` | Side lengths `(dx dy dz)` |

```bash
shore body box --extents 1 2 3 -o box.stl
```

## `shore body capsule`

Cylinder with hemispherical end caps, axis along z. Total length
along z is `height + 2 * radius`.

```
shore body capsule -o FILE.stl [--radius FLOAT] [--height FLOAT] [--sections INTEGER]
```

| Option | Default | Description |
|---|---|---|
| `--radius FLOAT` | 0.5 | Cap and cylinder radius |
| `--height FLOAT` | 1.0 | Cylindrical middle length (caps add `2 * radius` more) |
| `--sections INTEGER` | 16 | Segments per circumferential / meridional axis |

## `shore body extrude`

Extrude a 2D Selig section file into a 3D wing STL. Supports
constant-section, tapered, swept, and twisted wings; any
combination.

```
shore body extrude --section FILE.dat -o FILE.stl [OPTIONS]
```

### Options

| Option | Default | Description |
|---|---|---|
| `--section PATH` | — (required) | Input Selig `.dat` (e.g. from `shore profile naca`) |
| `-o` / `--output PATH` | — (required) | Output `.stl` file |
| `--span FLOAT` | 1.0 | Total wing span along the span axis |
| `--n-stations INTEGER` | 8 | Span stations including root and tip |
| `--taper FLOAT` | 1.0 | Tip-to-root chord ratio (1.0 = no taper) |
| `--sweep-deg FLOAT` | 0.0 | Sweep angle of the leading edge in degrees (positive = sweepback) |
| `--twist-deg FLOAT` | 0.0 | Geometric twist at the tip in degrees (linearly interpolated from 0 at root). Positive = nose-up. |
| `--span-axis TEXT` | `z` | Span axis: `x`, `y`, or `z` |
| `--end-caps` / `--open-ends` | `--open-ends` | Add centroid-fan caps at the span ends (watertight). Default is open ends — the [C-grid topology](/topologies/c-grid) needs them open. |

### Examples

```bash
# Constant-section straight wing (open ends, ready for shore mesh --topology cgrid):
shore profile naca 0012 -o naca0012.dat
shore body extrude --section naca0012.dat --span 1.0 --n-stations 8 -o wing.stl

# Tapered swept twisted wing with watertight ends:
shore body extrude --section naca0012.dat --span 2.0 --n-stations 12 \
    --taper 0.6 --sweep-deg 25 --twist-deg -3 --end-caps -o wing.stl

# UIUC database section (download a .dat file first):
shore body extrude --section eppler395.dat --span 1.5 --n-stations 10 -o e395.stl

# Span along x instead of z:
shore body extrude --section naca0012.dat --span-axis x --span 1.0 -o wing.stl
```

### Geometry conventions

- The 2D section lives in the plane of the two non-span axes.
  By default (`--span-axis z`) the chord direction is `+x` and
  thickness is along `±y`.
- Stations are equally spaced from `0` to `--span` along the
  chosen axis.
- Taper scales each station's section about its leading edge.
- Sweep shifts each station along the chord axis by
  `tan(sweep_deg) * z`.
- Twist rotates each station about its leading edge by an angle
  linearly interpolated from `0` at the root to `--twist-deg` at
  the tip.
- For the [C-grid topology](/topologies/c-grid), default to
  `--open-ends` so the per-station slicer doesn't trip over end-cap
  triangles.

## See also

- [`shore profile`](profile) — generate 2D NACA sections.
- [C-grid topology](/topologies/c-grid) — the natural meshing target for extruded wings.
- [Cubed-sphere](/topologies/cubed-sphere) and [O-grid](/topologies/o-grid) — meshing targets for sphere / cylinder / torus / capsule bodies.
