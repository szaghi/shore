---
title: shore.config
outline: [2, 3]
---

# `shore.config`

TOML-driven configuration for the `shore mesh` pipeline.

A config file is a *single-body* meshing recipe.  Every key has a default
that mirrors the corresponding `shore mesh` CLI flag.  CLI flags
override config-file values; config-file values override built-in
defaults.

Schema (single body, `version = 1`)
-------------------------------------

```
version = 1                              # required, integer

input  = "fuselage.stl"                  # required, path to STL
output = "fuselage.geo"                  # optional output stem

[grid]
ni = 60                                  # latitude  resolution
nj = 80                                  # longitude resolution (multiple of 4)
nk = 40                                  # wall-normal layers
theta_cap_deg = 30.0                     # polar exclusion (degrees)

[march]
ds            = 1.0e-4                   # first-layer thickness
growth        = 1.15                     # geometric growth ratio
smooth        = 0.2                      # Laplacian blend strength
smooth_iters  = 2                        # Laplacian sweeps per layer
smooth_mode   = "limited"                # limited | legacy (unlimited, folds thin layers)
march_scheme  = "algebraic"              # algebraic | hyperbolic (Steger-Chaussee)

[surface]
i_spacing = "uniform"                    # uniform | tanh | tanh2
i_beta    = 3.0
j_spacing = "uniform"                    # uniform | tanh2
j_beta    = 3.0

[volume]
k_spacing   = "geometric"                # geometric | tanh | tanh2
k_beta      = 3.0
k_thickness = 5.0                        # required for tanh / tanh2

```
Unknown keys raise `ParameterError`: the schema is closed so that
typos surface immediately rather than being silently ignored.

## `GridConfig` <Badge type="tip" text="dataclass" /> {#GridConfig}

```python
GridConfig(
    ni: int = 40,
    nj: int = 60,
    nk: int = 30,
    theta_cap: float | None = None,
    topology: str = 'cubed_sphere',
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L97">source</a>

## `MarchConfig` <Badge type="tip" text="dataclass" /> {#MarchConfig}

```python
MarchConfig(
    ds: float = 0.001,
    growth: float = 1.15,
    smooth: float = 0.2,
    smooth_iters: int = 2,
    smooth_mode: str = 'limited',
    march_scheme: str = 'algebraic',
    blend_normals_k: int | None = None,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L106">source</a>

## `SurfaceConfig` <Badge type="tip" text="dataclass" /> {#SurfaceConfig}

```python
SurfaceConfig(
    i_spacing: str = 'uniform',
    i_beta: float = 3.0,
    j_spacing: str = 'uniform',
    j_beta: float = 3.0,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L117">source</a>

## `VolumeConfig` <Badge type="tip" text="dataclass" /> {#VolumeConfig}

```python
VolumeConfig(
    k_spacing: str = 'geometric',
    k_beta: float = 3.0,
    k_thickness: float | None = None,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L125">source</a>

## `GeometryConfig` <Badge type="tip" text="dataclass" /> {#GeometryConfig}

```python
GeometryConfig(
    type: str = 'box',
    vmin: tuple[float, float, float] | None = None,
    vmax: tuple[float, float, float] | None = None,
    axis: str | tuple[float, float, float] = 'z',
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    z0: float | None = None,
    z1: float | None = None,
    r_in: float | None = None,
    r_out: float | None = None,
    theta_start: float = 0.0,
    theta_end: float | None = None,
    inner_frac: float | None = None,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L132">source</a>

Geometric inputs for a primitive.

The `type` field selects which primitive: `"box"`, `"annulus"`,
or `"flat_caps"`.  Only the keys relevant to that primitive are
consulted; the rest are ignored *but still validated* against the
closed schema (typos surface).

## `PrimSpacingConfig` <Badge type="tip" text="dataclass" /> {#PrimSpacingConfig}

```python
PrimSpacingConfig(
    i_law: str = 'uniform',
    j_law: str = 'uniform',
    k_law: str = 'uniform',
    i_beta: float = 3.0,
    j_beta: float = 3.0,
    k_beta: float = 3.0,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L162">source</a>

Per-axis spacing laws for primitives.

Box / annulus / flat_caps share the same three axes (i, j, k for
box and annulus; i, c, k for flat_caps where `c` maps to the
central-square / arc direction — written as `j` in the config
for visual consistency with the other primitives).

## `PrimGridConfig` <Badge type="tip" text="dataclass" /> {#PrimGridConfig}

```python
PrimGridConfig(ni: int = 40, nj: int = 40, nk: int = 30)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L180">source</a>

Node counts for primitives.

## `PrimitiveConfig` <Badge type="tip" text="dataclass" /> {#PrimitiveConfig}

```python
PrimitiveConfig(
    output: Path | None = None,
    geometry: GeometryConfig = field(default_factory=GeometryConfig),
    grid: PrimGridConfig = field(default_factory=PrimGridConfig),
    spacing: PrimSpacingConfig = field(default_factory=PrimSpacingConfig),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L189">source</a>

Resolved configuration for a single `shore primitive` invocation.

## `SplitDirective` <Badge type="tip" text="dataclass" /> {#SplitDirective}

```python
SplitDirective(label: str, axis: str, at: tuple[int, ...])
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L199">source</a>

One split directive parsed from a TOML split config.

Mirrors :class:`shore.split.BlockSplit`; kept separate so the
config layer doesn't import the split module (avoids a circular
import for users who only want config validation).

## `SplitConfig` <Badge type="tip" text="dataclass" /> {#SplitConfig}

```python
SplitConfig(
    output: Path | None = None,
    splits: list[SplitDirective] = field(default_factory=list),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L213">source</a>

Resolved configuration for a single `shore split` invocation.

## `MeshConfig` <Badge type="tip" text="dataclass" /> {#MeshConfig}

```python
MeshConfig(
    input: Path | None = None,
    output: Path | None = None,
    grid: GridConfig = field(default_factory=GridConfig),
    march: MarchConfig = field(default_factory=MarchConfig),
    surface: SurfaceConfig = field(default_factory=SurfaceConfig),
    volume: VolumeConfig = field(default_factory=VolumeConfig),
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L221">source</a>

Resolved configuration for a single `shore mesh` invocation.

All sections are present with defaults; `input` is required.

## `load_config` <Badge type="info" text="function" /> {#load_config}

```python
load_config(path: str | Path) -> MeshConfig
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L235">source</a>

Load and validate a SHORE TOML config.

**Parameters**

- **`path`** — Path to the `.toml` file.

**Returns**

- `MeshConfig` — Fully populated config; missing keys take their dataclass defaults.

**Raises**

- `ParameterError` — On file-not-found, malformed TOML, missing required keys, unknown keys, or version mismatch.

## `dump_template` <Badge type="info" text="function" /> {#dump_template}

```python
dump_template(stream=None) -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L372">source</a>

Write an annotated single-body TOML template to *stream* (default: stdout).

The template uses the dataclass defaults so it stays in sync with
the code automatically.

## `load_primitive_config` <Badge type="info" text="function" /> {#load_primitive_config}

```python
load_primitive_config(path: str | Path) -> PrimitiveConfig
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L455">source</a>

Load and validate a SHORE TOML primitive config.

A primitive config is identified by `kind = "primitive"` at the
top level.  Mesh and primitive configs share the same versioning
but otherwise have disjoint schemas.

**Parameters**

- **`path`** — Path to the `.toml` file.

**Returns**

- `PrimitiveConfig` — Fully populated config; missing keys take their dataclass defaults.

**Raises**

- `ParameterError` — On file-not-found, malformed TOML, missing required keys, unknown keys, version mismatch, or unknown geometry type.

## `dump_primitive_template` <Badge type="info" text="function" /> {#dump_primitive_template}

```python
dump_primitive_template(stream=None, type_: str = 'box') -> None
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L555">source</a>

Write an annotated primitive TOML template to *stream*.

*type_* is one of `"box"`, `"annulus"`, `"flat_caps"` and
determines which geometry block is emitted.

## `load_split_config` <Badge type="info" text="function" /> {#load_split_config}

```python
load_split_config(path: str | Path) -> SplitConfig
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/config.py#L618">source</a>

Load and validate a SHORE TOML split config.

A split config is identified by `kind = "split"` at the top
level.  Schema (single example):

```
version = 1
kind    = "split"

# Optional output stem; the CLI supplies the rest.
# output = "split_out/mesh"

[[splits]]
label = "sub0"
axis  = "k"
at    = [4]

[[splits]]
label = "background_center"
axis  = "i"
at    = [8, 16, 24]
```

**Returns**

- `SplitConfig` — Fully populated config.

**Raises**

- `ParameterError` — On file-not-found, malformed TOML, missing required keys, unknown keys, or version / kind mismatch.
