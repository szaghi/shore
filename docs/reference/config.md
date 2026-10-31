# Config files

For complex meshing scenarios — many spacing laws active at once,
non-default $\theta_{\text{cap}}$, a non-default marching scheme — the
flat CLI surface becomes hard to read. SHORE supports an alternative
**TOML input file** that you pass via `-c / --config` to `shore mesh`.

Three config kinds share the same loader rules (closed schema, `version = 1`,
paths relative to the file): the **mesh** config of this page (`kind =
"mesh"`, the default), the **primitive** config of `shore primitive
box|annulus|flat-caps -c` (`kind = "primitive"`, see
[primitive TOML schema](primitive#toml-schema-primitive-configs)) and the
**split** config of `shore split -c` (`kind = "split"`, see
[split TOML schema](split#toml-schema-kind-split)). `shore config` handles
the mesh kind only.

The schema is intentionally small: every key maps 1:1 to an existing
CLI flag. Nothing in the file is interpreted, expanded, or templated.
What you write is what runs.

## Quick start

Write a config to a file and run:

```bash
shore config template > shore.toml     # bootstrap with sane defaults
$EDITOR shore.toml                     # edit grid / spacing / etc.
shore config validate shore.toml       # parse-check before running
shore mesh -c shore.toml               # full pipeline, no other flags
```

A minimal valid config is just two keys:

```toml
version = 1
input   = "fuselage.stl"
```

Everything else (grid resolution, spacing laws, smoothing) falls back
to the same built-in defaults that `shore mesh` uses.

## Precedence

When both a config file and CLI flags are passed, the rule is:

> **explicit CLI flag → config-file value → built-in default**

A flag counts as "explicit" only if you actually typed it on the
command line. CLI flags left at their default are *not* treated as
overrides — that is what lets `shore mesh -c shore.toml` work without
every default leaking through and clobbering the file.

```bash
# config sets ni=80; CLI default ni=40 is NOT used; effective ni = 80
shore mesh -c shore.toml

# config sets ni=80; explicit --ni 100 overrides it; effective ni = 100
shore mesh -c shore.toml --ni 100
```

The positional `INPUT` argument follows the same rule: if the config
has `input = "..."` and you don't pass a positional, the file's value
is used. Passing a positional STL on the command line overrides the
config's `input`.

## Schema (version 1)

```toml
version = 1                          # required, integer

kind    = "mesh"                     # optional; the default
input   = "body.stl"                 # required, path to the input STL or CAD file
output  = "body.geo"                 # optional, output stem for the 6 .geo files

[grid]
ni            = 40                   # latitude (equator meridional)
nj            = 60                   # longitude (multiple of 4)
nk            = 30                   # wall-normal layers
theta_cap_deg = 30.0                 # cap angle (degrees); default 30 (cubed sphere), 9 (O-grid)
# theta_cap   = 0.5236               # alternative: in radians (mutually exclusive)
topology      = "cubed_sphere"       # 'cubed_sphere' (default) or 'ogrid' only

[march]
ds              = 1.0e-3             # first-layer thickness (geometric)
growth          = 1.15               # geometric growth ratio
smooth          = 0.2                # Laplacian blend strength [0–1]
smooth_iters   = 2                   # Laplacian sweeps per layer
smooth_mode    = "limited"           # limited | legacy
march_scheme   = "algebraic"         # algebraic | hyperbolic

[surface]
i_spacing = "uniform"                # uniform | tanh | tanh2
i_beta    = 3.0                      # clustering strength
j_spacing = "uniform"                # uniform | tanh2
j_beta    = 3.0

[volume]
k_spacing   = "geometric"            # geometric | tanh | tanh2
k_beta      = 3.0
k_thickness = 5.0                    # required when k_spacing is tanh / tanh2
```

::: tip Closed schema
Unknown keys at the top level *or* in any section raise
`ParameterError` immediately. This is deliberate: a typo like
`grwoth = 1.15` should fail loudly, not silently fall back to the
default. Run `shore config validate <file>` to surface typos before a
real run.
:::

## Path resolution

Relative paths in `input` and `output` are resolved against the
**directory containing the config file**, not the current working
directory:

```toml
# /home/me/cases/fuselage/shore.toml
input  = "../meshes/raw.stl"          # → /home/me/cases/meshes/raw.stl
output = "out/fuselage"               # → /home/me/cases/fuselage/out/fuselage
```

This makes a config file portable: copy the directory, the paths still
resolve.

## `shore config` subcommand

`shore config <action>` is a small inspection utility:

| Action | Purpose |
|---|---|
```
shore config [OPTIONS] ACTION [PATH]
```

| Argument | Description |
|---|---|
| `ACTION` | `template` or `validate` (required). |
| `PATH` | The config to validate; required for `validate`, ignored for `template`. |

| Action | Purpose |
|---|---|
| `template` | Write an annotated default mesh config to stdout. Pipe to a file: `shore config template > shore.toml`. |
| `validate PATH` | Parse and validate a mesh config; print the resolved values. Useful in CI. |

`validate` exits 1 on any schema violation, with a single human-readable
error message; an unknown action also exits 1. For example, on the
template with `input` pointing at an existing STL:

```
✓ shore.toml parses cleanly
  input         = /home/me/cases/sphere.stl
  grid          ni=40 nj=60 nk=30 topology='cubed_sphere'
  march         ds=0.001 growth=1.15 smooth=0.2/2
  surface       i=uniform(3.0) j=uniform(3.0)
  volume        k=geometric(3.0) thickness=None
```

## Field reference

Every key has the same semantics as the corresponding `shore mesh` CLI
flag. The CLI is the source of truth for definitions; this table is
just the alias map.

| TOML key | CLI flag |
|---|---|
| `input` | `INPUT` (positional) |
| `output` | `-o` / `--output` |
| `grid.ni` | `--ni` |
| `grid.nj` | `--nj` |
| `grid.nk` | `--nk` |
| `grid.theta_cap` | `--theta-cap` |
| `grid.theta_cap_deg` | `--theta-cap-deg` |
| `grid.topology` | `--topology` (`cubed_sphere` or `ogrid`; the C-grid and hybrid topologies need CLI flags that have no config keys) |
| `march.ds` | `--ds` |
| `march.growth` | `--growth` |
| `march.smooth` | `--smooth` |
| `march.smooth_iters` | `--smooth-iters` |
| `march.smooth_mode` | `--smooth-mode` |
| `march.march_scheme` | `--march-scheme` |
| `march.blend_normals_k` | `--blend-normals-k` (removed, [ADR 0034](/adr/0034-remove-normal-blending): only `0` is accepted, with a notice) |
| `surface.i_spacing` | `--surface-i-spacing` |
| `surface.i_beta` | `--surface-i-beta` |
| `surface.j_spacing` | `--surface-j-spacing` |
| `surface.j_beta` | `--surface-j-beta` |
| `volume.k_spacing` | `--volume-k-spacing` |
| `volume.k_beta` | `--volume-k-beta` |
| `volume.k_thickness` | `--volume-k-thickness` |

## Worked example

A config representative of a real boundary-layer mesh — fixed
far-field thickness, tanh wall-normal clustering, equator clustering
at the cap seams:

```toml
# fuselage.toml
version = 1

input  = "fuselage.stl"
output = "out/fuselage"

[grid]
ni            = 80
nj            = 120
nk            = 50
theta_cap_deg = 30.0

[march]
ds              = 1.0e-5      # ignored when k_spacing != geometric
growth          = 1.15        # ignored when k_spacing != geometric
smooth          = 0.3
smooth_iters    = 3

[surface]
i_spacing = "tanh"
i_beta    = 4.0
j_spacing = "uniform"
j_beta    = 3.0

[volume]
k_spacing   = "tanh"
k_beta      = 4.0
k_thickness = 5.0
```

Run with:

```bash
shore mesh -c fuselage.toml -v
```

A ready-to-run version of this file ships at
[`examples/01-body-fitted/config_fuselage.toml`](https://github.com/szaghi/shore/blob/main/examples/01-body-fitted/config_fuselage.toml).

## Limitations

- **Single body, two topologies**. The schema holds one input file and
  the `cubed_sphere` and `ogrid` topologies; `cgrid`, `ch` and `oh` are
  CLI-only. Multi-component runs are built as projects
  (`--project DIR --component C`, [`shore xall`](xall)).
- **No project options.** `--project`, `--component`, `--journal` and
  `--adjacency-json` have no config keys; pass them on the command line.
- **No includes / inheritance / templating**. By design — share
  defaults via a Python script that emits TOML if you need it.
- **No environment-variable interpolation.** Paths resolve against
  the config file's directory; that's the only context-dependent
  behaviour.

## See also

- [`shore mesh`](mesh) — the command this config drives
- [Quick start](/guide/quickstart) — the same workflow with bare CLI
- [Spacing laws](/concepts/spacing-laws) — what the spacing knobs do
