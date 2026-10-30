---
title: shore.collar
outline: [2, 3]
---

# `shore.collar`

The collar block of two intersecting bodies (ADR 0010).

`build_collar` wraps one Chimera volume block around the junction of two
bodies:

1. **Wall.** `intersection_loop` (ADR 0009) gives the junction curve and
   `march_surface` (ADR 0008) grows a surface grid on each body away from
   the other; the two halves meet at the junction row `i = len(steps_b)`.
   j runs around the loop (periodic, seam node not repeated), i from B's far
   edge to A's far edge, k leaves the walls.
2. **Near-wall layers.** The hyperbolic scheme (ADR 0005) while the layers
   are thinner than half the junction spacing *h*. Beyond that, orthogonal
   offsetting collapses the cells beside the concave crease (at a 90 deg
   junction the cell next to it has width `h - d` at thickness *d*).
3. **Outer boundary.** Per loop station, the near-wall i-line is offset by
   the remaining thickness in its best-fit plane. The swallowtail that the
   concave junction makes is trimmed at its self-intersection, and the
   nodes are placed by arc length with the trimmed crease mapped to the
   junction row.
4. **Interior.** Per loop station, TFI between the near-wall line and the
   outer boundary with the k schedule, then Winslow smoothing with
   Thomas-Middlecoff control functions (the clustering of the boundaries is
   kept, as in ADR 0007), solved by Picard iteration with SuperLU.

Steps 2-4 are `extrude_collar`. With `plane` (ADR 0013) the collar is
the half-collar of a half model: the open curve of `intersection_curve`
(ADR 0011), both halves marched with their ends on the plane (ADR 0012),
and the volume extruded from the wall mirrored across the plane into a
closed loop, of which the first `n_loop` stations are kept; the j-faces
are put on the plane and marked `FaceBC.SYMMETRY`. `loop_spacing` and
`loop_curvature` distribute the nodes along the junction curve (ADR 0020).

`split_collar` (ADR 0024) is the collar of a loop with sharp corners
(turning more than `CollarParams.sharp_deg`: a sharp trailing edge), where
the surface march around the corner folds. The loop is cut at each sharp
corner into pieces. Each piece is extended at both ends by a wake line on A
along the corner's bisector, and the wall over the wake is a virtual sheet
B': the corner's B-line translated along the wake. The two sides of every
sheet are joined by a SHARED seam (a C-grid wake cut). Each piece is
extruded as one open block and split into six (wake, body, wake by B, A) so
that every face has one boundary condition. Loops without sharp corners
give `build_collar`'s block.

## `CollarParams` <Badge type="tip" text="dataclass" /> {#CollarParams}

```python
CollarParams(
    hyperbolic_fraction: float = 0.5,
    max_iterations: int = 400,
    tolerance: float = 1e-12,
    control_limit: float = 1.5,
    oblique_deg: float = 1.0,
    sharp_deg: float = 105.0,
    wake_ratio: float = 1.2,
    volume: str = 'smoothed',
    offset_smoothing: float = 0.5,
    backstep_ratio: float = 0.5,
    backstep_iterations: int = 200,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/collar.py#L68">source</a>

Collar settings (not per-case tuning knobs).

## `build_collar` <Badge type="info" text="function" /> {#build_collar}

```python
build_collar(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    *,
    n_loop: int,
    steps_a: np.ndarray,
    steps_b: np.ndarray,
    steps_k: np.ndarray,
    label: str = 'collar',
    params: CollarParams = DEFAULT,
    plane: tuple[np.ndarray, np.ndarray] | None = None,
    loop_spacing: Spacing1D | None = None,
    loop_curvature: float = 0.0,
) -> HexBlock
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/collar.py#L89">source</a>

The collar block around the junction of *body_a* and *body_b*.

**Parameters**

- **`body_a, body_b`** — The two bodies (`Geometry.project` with normals).
- **`seed`** (np.ndarray, shape (3,)) — A point near the junction; selects the loop when there are several.
- **`n_loop`** — Nodes around the loop (j), &gt;= 3.
- **`steps_a, steps_b`** — Surface marching steps away from the junction on A and on B.
- **`steps_k`** — Wall-normal steps; the collar is `sum(steps_k)` thick.
- **`label`** — Block label.
- **`params`** — Collar settings.
- **`plane`** — `(point, normal)` of a symmetry plane: build the half-collar of a half model, on the seed's side. `None`: the closed collar.
- **`loop_spacing, loop_curvature`** — Node distribution along the junction curve, passed to `intersection_loop` / `intersection_curve` as *spacing* and *curvature* (ADR 0020; default: equally spaced in arc length).

**Returns**

- `HexBlock` — Nodes `(len(steps_b) + len(steps_a) + 1, n_loop, len(steps_k) + 1, 3)`; face `k_lo` is WALL. Closed collar: `periodic_j`, the other faces FREE. Half-collar: open in j, faces `j_lo` / `j_hi` on the plane and SYMMETRY, the others FREE.

**Raises**

- `ParameterError` — Malformed steps or plane, `n_loop < 3`, or no intersection near *seed*.

**Warns**

- `UserWarning` — The half-collar's curve meets the plane obliquely (its end tangent more than `params.oblique_deg` off the plane normal): the bodies are not mirror-symmetric there, and the block's quality is not guaranteed.

## `split_collar` <Badge type="info" text="function" /> {#split_collar}

```python
split_collar(
    body_a: Geometry,
    body_b: Geometry,
    seed: np.ndarray,
    *,
    n_loop: int,
    steps_a: np.ndarray,
    steps_b: np.ndarray,
    steps_k: np.ndarray,
    label: str = 'collar',
    params: CollarParams = DEFAULT,
    loop_spacing: Spacing1D | None = None,
    loop_curvature: float = 0.0,
    wake_length: float | None = None,
) -> list[HexBlock]
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/collar.py#L193">source</a>

The collar of a junction whose loop has sharp corners, cut there (ADR 0024).

**Parameters**

- **`body_a, body_b, seed, n_loop, steps_a, steps_b, steps_k, label, params, loop_spacing, loop_curvature`** — As for `build_collar` (closed loops only).
- **`wake_length`** — Length of the wake lines on A behind each sharp corner (default `sum(steps_a)`, the width of the collar on A).

**Returns**

- `list[HexBlock]` — Without sharp corners, `[build_collar(...)]`. Otherwise six blocks per piece *p* of the loop between consecutive sharp corners, `f"{label}_{p}_{part}_{side}"` with *part* `wake_lo` / `body` / `wake_hi` (j: the wake behind the piece's last corner, the piece, the wake behind its first corner) and *side* `b` / `a` (i below / above the junction row). `k_lo` is WALL, except on the `wake_*_b` blocks: the sheet behind a corner, SHARED with the other side's sheet block (along i same, along j reversed). The seams inside a piece are SHARED; the other faces FREE.

**Raises**

- `ParameterError` — As `build_collar`; a bad *wake_length*; a sharp corner that is not an edge of B (on A's edge the wake has no surface to run on).

## `extrude_collar` <Badge type="info" text="function" /> {#extrude_collar}

```python
extrude_collar(
    wall: np.ndarray,
    steps_k: np.ndarray,
    junction: int,
    *,
    params: CollarParams = DEFAULT,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/collar.py#L383">source</a>

The collar volume over a j-periodic *wall* (steps 2-4 of the method).

**Parameters**

- **`wall`** (np.ndarray, shape (ni, nj, 3)) — The collar wall: i from B's far edge through the junction row to A's far edge, j around the loop (periodic, seam node not repeated).
- **`steps_k`** — Wall-normal steps.
- **`junction`** — The junction row (`len(steps_b)`), 0 &lt; junction &lt; ni - 1.
- **`params`** — Collar settings.

**Returns**

- `np.ndarray, shape (ni, nj, len(steps_k) + 1, 3)` — `[:, :, 0]` is *wall*.
