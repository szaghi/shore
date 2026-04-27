# shore primitive

Build a parametric structured-grid primitive: a Cartesian box, an
annular cylinder, or a 5-block flat-caps cylinder.

`shore primitive` is a sub-app with three commands:

```
shore primitive box [OPTIONS]
shore primitive annulus [OPTIONS]
shore primitive flat-caps [OPTIONS]
```

Every command accepts the same `-c / --config` flag for TOML-driven
invocation and follows the same precedence rule as `shore mesh`:
**explicit CLI flag → config-file value → built-in default**.

See [Primitives](/topologies/primitives) for the conceptual overview
and use cases.

## `shore primitive box`

Cartesian axis-aligned box, single block.

### Synopsis

```
shore primitive box [OPTIONS]
```

### Options

| Option | Default | Description |
|---|---|---|
| `-c` / `--config PATH` | — | TOML primitive config (`type = "box"`) |
| `-o` / `--output PATH` | `box.geo` | Output `.geo` file |
| `--min FLOAT FLOAT FLOAT` | `0 0 0` | Minimum corner of the box |
| `--max FLOAT FLOAT FLOAT` | `1 1 1` | Maximum corner; each component must strictly exceed the corresponding `--min` component |
| `--ni INTEGER` | 40 | Node count along i (x-axis) |
| `--nj INTEGER` | 40 | Node count along j (y-axis) |
| `--nk INTEGER` | 30 | Node count along k (z-axis) |
| `--i-spacing TEXT` | `uniform` | Spacing along i: `uniform` / `tanh` / `tanh2` |
| `--j-spacing TEXT` | `uniform` | Spacing along j |
| `--k-spacing TEXT` | `uniform` | Spacing along k |
| `--i-beta FLOAT` | 3.0 | Clustering strength for tanh / tanh2 along i |
| `--j-beta FLOAT` | 3.0 | Clustering strength for tanh / tanh2 along j |
| `--k-beta FLOAT` | 3.0 | Clustering strength for tanh / tanh2 along k |

### Examples

```bash
# Cubic background, uniform spacing
shore primitive box --min 0 0 0 --max 10 10 10 --ni 40 --nj 40 --nk 40 -o bg.geo

# Boundary-layer-style cluster near the inflow plane
shore primitive box --min 0 0 0 --max 10 5 3 \
                    --i-spacing tanh --i-beta 4.0 \
                    --ni 80 --nj 30 --nk 20 -o channel.geo
```

## `shore primitive annulus`

Annular cylinder, optionally a partial sector. Single block.

### Synopsis

```
shore primitive annulus [OPTIONS]
```

### Options

| Option | Default | Description |
|---|---|---|
| `-c` / `--config PATH` | — | TOML primitive config (`type = "annulus"`) |
| `-o` / `--output PATH` | `annulus.geo` | Output `.geo` file |
| `--r-in FLOAT` | 0.5 | Inner radius |
| `--r-out FLOAT` | 2.0 | Outer radius (must exceed `--r-in`) |
| `--z0 FLOAT` | 0.0 | Axial start along `--axis` |
| `--z1 FLOAT` | 1.0 | Axial end along `--axis` |
| `--axis TEXT` | `z` | Cylinder axis: `x` / `y` / `z` or `"vx vy vz"` |
| `--center FLOAT FLOAT FLOAT` | `0 0 0` | World-space anchor of the axis |
| `--theta-start FLOAT` | 0.0 | Sector start in radians |
| `--theta-end FLOAT` | — | Sector end in radians; omit for a full annulus |
| `--ni INTEGER` | 20 | Radial node count |
| `--nj INTEGER` | 64 | Azimuthal node count |
| `--nk INTEGER` | 30 | Axial node count |
| `--i-spacing` / `--j-spacing` / `--k-spacing` | `uniform` | Per-axis spacing law |
| `--i-beta` / `--j-beta` / `--k-beta` | 3.0 | Per-axis clustering strength |

::: warning Full annulus → uniform j-spacing
For a full annulus the azimuthal direction is periodic and the
distribution must be uniform. `--j-spacing tanh` / `tanh2` is rejected;
use a partial sector instead.
:::

### Examples

```bash
# Full annulus around a near-body sphere mesh
shore primitive annulus --r-in 1.5 --r-out 5.0 --z0 -2 --z1 2 \
                        --ni 30 --nj 64 --nk 40 -o buffer.geo

# Half-pipe sector with non-uniform azimuthal spacing
shore primitive annulus --r-in 1.5 --r-out 5.0 --z0 -2 --z1 2 \
                        --theta-start 0 --theta-end 3.14159 \
                        --j-spacing tanh2 --j-beta 3.0 \
                        --ni 24 --nj 32 --nk 40 -o halfpipe.geo

# Cylinder along x-axis instead of z
shore primitive annulus --axis x --r-in 0.5 --r-out 2 --z0 0 --z1 4 \
                        -o axis_x.geo
```

## `shore primitive flat-caps`

5-block flat-caps cylinder. The disk interior at every axial slice is
decomposed into a central square + 4 trapezoidal sub-blocks. Writes
five `.geo` files.

### Synopsis

```
shore primitive flat-caps [OPTIONS]
```

### Options

| Option | Default | Description |
|---|---|---|
| `-c` / `--config PATH` | — | TOML primitive config (`type = "flat_caps"`) |
| `-o` / `--output PATH` | `flat_caps` | Output stem; writes `<stem>_center.geo` + `<stem>_sub_e.geo` / `_sub_n.geo` / `_sub_w.geo` / `_sub_s.geo` |
| `--r-out FLOAT` | 2.0 | Outer radius |
| `--z0` / `--z1 FLOAT` | 0 / 1 | Axial extent along `--axis` |
| `--axis TEXT` | `z` | Cylinder axis: `x` / `y` / `z` or `"vx vy vz"` |
| `--center FLOAT FLOAT FLOAT` | `0 0 0` | World-space anchor of the axis |
| `--inner-frac FLOAT` | 0.5 | Central square corner radius / outer radius. Must lie in $(0, 1)$. |
| `--ni INTEGER` | 20 | Radial node count for each sub-block (square edge → outer circle) |
| `--nc INTEGER` | 20 | Side node count of the central square (shared between center and sub-blocks) |
| `--nk INTEGER` | 30 | Axial node count |
| `--i-spacing` / `--c-spacing` / `--k-spacing` | `uniform` | Per-axis spacing law |
| `--i-beta` / `--c-beta` / `--k-beta` | 3.0 | Per-axis clustering strength |

The `--c-spacing` controls the distribution along the central-square
sides and the matching arc on the outer circle (it is shared between
the center block and every sub-block to guarantee SHARED conformity
at every seam).

The sub-block radial direction has its **first cell pinned per-j** to
the central square's local tangential cell at column j, giving
aspect-ratio-1 cells at the central / sub-block seam regardless of
`--c-spacing`. When the pin would produce a non-monotone distribution
(extreme `inner_frac` or mismatched `ni` / `nc`), SHORE falls back to
the unpinned distribution and emits a `UserWarning`.

### Examples

```bash
# 5-block conforming cylindrical background
shore primitive flat-caps --r-out 5.0 --z0 -2 --z1 2 \
                          --ni 24 --nc 32 --nk 40 -o bg

# x-aligned, with axial clustering at the inflow plane
shore primitive flat-caps --axis x --r-out 5.0 --z0 0 --z1 10 \
                          --k-spacing tanh --k-beta 4.0 \
                          --ni 24 --nc 32 --nk 50 -o tunnel
```

## TOML schema (primitive configs)

Primitive configs are identified by `kind = "primitive"` at the top
level. Mesh and primitive configs share `version = 1` but otherwise
have disjoint schemas — feeding a mesh config to `shore primitive` (or
vice versa) is rejected with a clear error.

```toml
version = 1
kind    = "primitive"

[geometry]
type   = "box"            # box | annulus | flat_caps

# Box-only:
vmin = [0.0, 0.0, 0.0]
vmax = [10.0, 5.0, 3.0]

# Annulus / flat_caps shared:
# axis   = "z"
# center = [0.0, 0.0, 0.0]
# z0     = 0.0
# z1     = 4.0

# Annulus-only:
# r_in        = 0.5
# r_out       = 2.0
# theta_start = 0.0
# theta_end   = 3.14159   # omit for full annulus

# Flat-caps-only:
# r_out, z0, z1 from above
# inner_frac = 0.5            # central square corner radius / r_out (in (0, 1))

[grid]
ni = 40
nj = 40                   # for flat_caps this maps to --nc
nk = 30

[spacing]
i_law = "uniform"
j_law = "uniform"         # for flat_caps this maps to --c-spacing
k_law = "uniform"
i_beta = 3.0
j_beta = 3.0
k_beta = 3.0

# output = "background.geo"   # optional
```

Unknown keys raise `ParameterError`. Run `shore config validate`
(planned for primitive configs) or just invoke the subcommand — the
loader runs first.

## See also

- [Primitives overview](/topologies/primitives)
- [`shore mesh`](mesh) — body-fitted O-grid (the primary near-body workflow)
- [Config files](/reference/config) — TOML schema and precedence rules
