# Primitives

Primitives are parametric structured grids built from simple
geometric inputs — boxes, annular cylinders, flat-caps cylinders. No
STL, no projection, no hyperbolic marching: each primitive is a tensor
product of three 1D distributions.

They exist for one reason: an overset/Chimera assembly is not just the
near-body mesh produced by `shore mesh`. You also need a **background**
mesh that fills the whole domain, and often one or more **buffer**
meshes that bridge the resolution jump between near-body and
background. Primitives are that connective tissue.

::: tip Topology vs. role
"Background", "buffer", and "near-body" are *roles* in the assembly,
not types. The same Cartesian box can be a background in one
assembly and a buffer in another — the role is determined by where
the user places it relative to the other meshes, not by the primitive
itself.
:::

## Available primitives

### `box` — single-block Cartesian background

A rectangular axis-aligned box, single hexahedral block, shape
`(ni, nj, nk, 3)`. Each axis has independent spacing
(`uniform` / `tanh` / `tanh2`). No projection.

```bash
shore primitive box --min 0 0 0 --max 10 5 3 --ni 60 --nj 30 --nk 20 -o background.geo
```

Use case: outermost background of a Chimera assembly; simple buffer
when the geometry permits axis-aligned cells.

### `annulus` — annular cylinder, optionally a partial sector

A cylindrical shell bounded by `r_in`, `r_out`, `z0`, `z1`. The
azimuthal direction is **periodic** when the annulus is full
($2\pi$); for partial sectors `theta_start` / `theta_end` are explicit
and the j direction is not periodic.

```bash
# full annulus around a sphere
shore primitive annulus --r-in 1.5 --r-out 5.0 --z0 -2 --z1 2 \
                        --ni 30 --nj 64 --nk 40 -o buffer.geo

# half-pipe (180° sector)
shore primitive annulus --r-in 1.5 --r-out 5.0 --z0 -2 --z1 2 \
                        --theta-start 0 --theta-end 3.14159 \
                        --ni 30 --nj 32 --nk 40 -o halfpipe.geo
```

::: warning Periodic-j requires uniform spacing
A full annulus drops the wrap-around node; the azimuthal distribution
must be uniform so successive cells are equal. For non-uniform
azimuthal spacing, use a partial sector with explicit
`theta_start` / `theta_end`.
:::

Use case: buffer between a near-body sphere/cylinder mesh and the
outer Cartesian background.

### `flat-caps` — 5-block cylinder, no central hole

Same outer geometry as `annulus`, but the disk interior at every
axial slice is **conformally meshed** by a central square block plus
4 trapezoidal sub-blocks. This eliminates the Chimera fringe a hollow
cylinder would otherwise need at its inner boundary.

```bash
shore primitive flat-caps --r-out 5.0 --z0 -2 --z1 2 \
                          --ni 24 --nc 32 --nk 40 -o background
# writes 5 files: background_center.geo, background_sub_e.geo, ...
```

Block layout in the (e1, e2) plane (rotated 45°, "diamond" orientation):

```
              sub_n
           ╱       ╲
        ╱──┐     ┌──╲
        │  │     │  │
        │  │ ctr │  │
  sub_w │  │     │  │ sub_e
        │  │     │  │
        ╲──┘     └──╱
           ╲       ╱
              sub_s
```

The 4 corners of the central square sit at radius
`r_inner = inner_frac * r_out` (default `inner_frac = 0.5`) at
azimuths $\pi/4, 3\pi/4, 5\pi/4, 7\pi/4$; **straight radial edges**
connect each corner to the matching point on the outer circle. Each
sub-block fills the wedge between two adjacent radial edges, with one
side on the central square and the opposite side on a 90° arc of the
outer circle. SHARED-memory conformity at every seam: the central
square's edges are bit-exact equal to each sub-block's inner edge,
and the radial line from a square corner to the circle is shared by
two adjacent sub-blocks.

The size of the central square is controlled by `--inner-frac`
(default 0.5):

| Value | Effect |
|---|---|
| Small (e.g. 0.2) | Tiny central square, large sub-blocks. Useful when most of the disk should be radially structured. |
| 0.5 (default) | Central square corners at half-radius — visually balanced, matches the cubed-sphere-topology look. |
| Large (e.g. 0.8) | Big central square, thin sub-blocks. Useful when most of the disk should be Cartesian. |

#### C1 first-step continuity at the central / sub-block seam

The sub-block radial direction is automatically pinned so that the
first radial cell at each j-column equals the central square's local
tangential cell at that column. This gives **aspect-ratio-1 cells at
the seam** on both sides — even when `--c-spacing` clusters the central
distribution non-uniformly. Per-j pinning, same discipline as the
cubed-sphere topology's per-meridian C1 pin in `surface/equator.py`. When
the chosen `ni`, `nc`, and `inner_frac` make the pin geometrically
incompatible with the chosen `--i-spacing` (first step would exceed
or be much smaller than the natural radial cell size), SHORE falls
back to the unpinned distribution and emits a `UserWarning`.

Use case: cylindrical background mesh when you want to avoid the
cost of Chimera fringes — high-order schemes that lose order at
overset interpolations, conservation-critical solvers, or GPU
deployments where the irregular memory access pattern of donor
lookups is a performance ceiling. The trade-off is loss of
flexibility: a flat-caps background must be *axis-aligned* with
whatever it surrounds, whereas an `annulus + Cartesian background`
overset stack imposes no such constraint.

## Axis customisation

`annulus` and `flat-caps` accept either a symbolic axis label or an
explicit 3-vector:

| Form | Meaning |
|---|---|
| `--axis x` | World x axis |
| `--axis y` | World y axis |
| `--axis z` | World z axis (default) |
| `--axis "1.0 1.0 0.0"` | Arbitrary direction; the primitive normalises it. |

In TOML, the same key accepts either a string or a 3-component array:

```toml
axis = "z"                   # symbolic
axis = [1.0, 1.0, 0.0]       # vector
```

The local frame `(e1, e2, axis)` is right-handed; `e1` is whichever
world axis is most perpendicular to *axis*. This choice is
deterministic but otherwise arbitrary — if you depend on a specific
azimuthal orientation, use `--center` to anchor the cylinder and
post-rotate as needed.

## Config-file workflow

Primitives use the same TOML config flow as `shore mesh`. The schema
is identified by `kind = "primitive"` at the top level so a single
file format unambiguously routes to the right loader:

```bash
shore config template primitive box      > box.toml       # not yet implemented as CLI
# Or hand-write a config:
```

```toml
version = 1
kind    = "primitive"

[geometry]
type   = "annulus"
r_in   = 1.5
r_out  = 5.0
z0     = -2.0
z1     =  2.0
axis   = "z"
center = [0.0, 0.0, 0.0]

[grid]
ni = 30
nj = 64
nk = 40

[spacing]
i_law  = "tanh"
i_beta = 4.0
j_law  = "uniform"
k_law  = "uniform"

output = "buffer.geo"
```

Run with:

```bash
shore primitive annulus -c buffer.toml
```

CLI flags override config values exactly as for `shore mesh`. See
[Config files](/reference/config) for the full precedence rule.

## Limitations

- **No mesh assembly logic.** Primitives are produced standalone;
  combining them into a Chimera assembly (hole-cutting, fringe
  identification, donor search) is its own project, not part of
  primitive generation.
- **`flat-caps` is axis-aligned with the cylinder axis only.** If
  your near-body mesh axis differs from the desired flat-caps
  alignment, use the `annulus + box` overset route instead.
- **No periodic flag on `box`.** Backgrounds rarely need it;
  periodic channels are a separate use case.
- **No multi-axis annulus.** One axis per cylinder.

## See also

- [`shore primitive` reference](/reference/primitive) — full CLI / TOML schema
- [Cubed-sphere topology](/topologies/cubed-sphere) — body-fitted near-body mesh
- [Config files](/reference/config) — TOML schema and precedence
