---
title: shore.pipeline
outline: [2, 3]
---

# `shore.pipeline`

Meshing pipeline steps shared by the CLI and the command layer (ADR 0004).

`build_topology` turns a topology name plus parameters into wired blocks.
It is the single code path behind `shore mesh` and
`Project.build_topology`, so a journal replay produces the same blocks as
the CLI run that recorded it.

## `k_spacing` <Badge type="info" text="function" /> {#k_spacing}

```python
k_spacing(
    law: str = 'geometric',
    *,
    ds: float = 0.001,
    growth: float = 1.15,
    beta: float = 3.0,
    thickness: float | None = None,
) -> Spacing1D
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L33">source</a>

Wall-normal spacing law, as `shore mesh --volume-k-*` builds it.

## `make_topology` <Badge type="info" text="function" /> {#make_topology}

```python
make_topology(name: str, options: dict[str, Any] | None = None) -> Topology
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L47">source</a>

Topology object for *name*; *options* are its constructor arguments.

`ch` takes `chord_min/max`, `normal_min/max`, `nk_channel`,
`freestream_world`, `span_axis`, `k_law_channel`,
`k_beta_channel`, `i_split_body_upper/lower`, `ni_outflow`,
`c1_seam`; `oh` takes `box_vmin`, `box_vmax`, `nk_channel`,
`k_law_channel`, `k_beta_channel`, `c1_seam`; the others take none.

## `build_topology` <Badge type="info" text="function" /> {#build_topology}

```python
build_topology(
    body: Body,
    name: str,
    *,
    nk: int,
    spacing: Spacing1D,
    options: dict[str, Any] | None = None,
    **params: Any,
) -> tuple[Topology, list[HexBlock]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L75">source</a>

Build the wired blocks of topology *name* on *body*.

*params* go to the topology's `build`: `ni`, `nj`, `theta_cap`,
`i_spacing`, `i_beta`, `j_spacing`, `j_beta` for the ray-cast
topologies (`cubed_sphere`, `ogrid`, `oh`); `ni_body`,
`ni_wake`, `n_stations`, `wake_length`, `wake_growth`,
`span_axis`, `te_dihedral_deg`, `blocks` and the `cap_*`
options for `cgrid`; the C-grid subset for `ch`.

## `build_primitive` <Badge type="info" text="function" /> {#build_primitive}

```python
build_primitive(kind: str, **params: Any) -> Any
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L128">source</a>

Build parametric primitive *kind* (`shore primitive`); returns the
primitive object (`.blocks`, `.write_geo`).

## `common_base_stem` <Badge type="info" text="function" /> {#common_base_stem}

```python
common_base_stem(stems: list[str]) -> str
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L152">source</a>

Longest common `<base>_` prefix of *stems*, ending at an underscore.

A single stem gives everything up to and including its last underscore.

## `load_blocks` <Badge type="info" text="function" /> {#load_blocks}

```python
load_blocks(
    geo_paths: list[str] | tuple[str, ...],
    adjacency_paths: list[str] | tuple[str, ...] = (),
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L169">source</a>

Blocks from `.geo` files, wired from adjacency sidecars (`shore split`).

Labels come from the file stems with their longest common `<base>_`
prefix removed (`sphere_sub0.geo` -&gt; `sub0`); a single file drops
everything up to its last underscore.

## `apply_splits` <Badge type="info" text="function" /> {#apply_splits}

```python
apply_splits(
    blocks: list[HexBlock],
    splits: list[tuple[str, str, tuple[int, ...]]],
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L195">source</a>

Split *blocks* by `(label, axis, at)` directives in one pass.

## `build_tip_cap` <Badge type="info" text="function" /> {#build_tip_cap}

```python
build_tip_cap(
    body_geo: str,
    *,
    tip: str,
    le_frac: float,
    te_frac: float,
    span_axis: str,
    tip_inset: float | None,
    nk: int | None,
    n_cap_j: int,
    cap_j_spacing: str,
    cap_j_beta: float,
    auto_round_te: bool,
    te_thickness: float | None,
    te_thickness_frac: float | None,
) -> tuple[Any, dict[str, Any]]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/pipeline.py#L210">source</a>

Wing-tip cap from a C-grid body block `.geo` (`shore cap-tip`).

Returns the cap `Mesh` and the values resolved automatically
(`tip_inset`: 1e-3 chord when *None*; `nk`: the body's when *None*).
