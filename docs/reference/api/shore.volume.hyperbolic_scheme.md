---
title: shore.volume.hyperbolic_scheme
outline: [2, 3]
---

# `shore.volume.hyperbolic_scheme`

Hyperbolic volume marching: one layer of the Steger-Chaussee scheme (ADR 0005).

The new layer `r_{k+1} = r_k + dr` satisfies, at every node,

    r_xi . r_zeta = 0,   r_eta . r_zeta = 0,   det[r_xi, r_eta, r_zeta] = dV,

(the k-lines are orthogonal to the layer, and the cell volume is prescribed).
Linearised about the local solution `r_zeta0 = dV * n / |n|^2` with
`n = r_xi x r_eta` (Steger & Chaussee 1980; Chan & Steger 1992), the system
for `dr` becomes

    (I + Ci A d_xi + Ci B d_eta - eps_i Lap) dr = r_zeta0 + eps_e T(Lap r_k)

with `A = [r_zeta0; 0; r_eta x r_zeta0]`, `B = [0; r_zeta0; r_zeta0 x r_xi]`,
`Ci = C^-1`, `C = [r_xi; r_eta; n]`, central differences `d`, the
5-point Laplacian `Lap` and `T` the projection onto the layer's tangent
plane. The coupling `Ci A`, `Ci B` scales with step / cell size: thin
near-wall layers march almost exactly along the local normal, thick layers
feel their neighbours, which spreads converging fronts at concave corners.

Stabilisation, all scaled by `rho = min(1, step / h)` (`h` the local
cell size), so the near-wall layers stay orthogonal and of exact thickness:

- volume smoothing: `dV` is blended toward its neighbour average;
- a grid-convergence sensor raises dissipation where neighbouring nodes of
  the local-solution layer approach each other (concave regions);
- implicit (`eps_i`) and tangential explicit (`eps_e`) second-difference
  dissipation.

Boundary rows (non-periodic edges) are Dirichlet: they take the local
solution (or the caller's *boundary_target*, for Dirichlet-pinned cap
edges). Seam continuity is the topology marcher's job, as for the algebraic
scheme. The linear system is solved per block layer with SuperLU
(deterministic).

## `HyperbolicParams` <Badge type="tip" text="dataclass" /> {#HyperbolicParams}

```python
HyperbolicParams(
    eps_implicit: float = 0.5,
    eps_explicit: float = 0.25,
    volume_smoothing: float = 0.2,
    volume_iters: int = 2,
    convergence_gain: float = 16.0,
    distance_clip: float = 2.0,
    excess_bound: float = 0.04,
)
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/hyperbolic_scheme.py#L49">source</a>

Defaults of the hyperbolic scheme (not per-case tuning knobs).

## `hyperbolic_layer` <Badge type="info" text="function" /> {#hyperbolic_layer}

```python
hyperbolic_layer(
    layer: np.ndarray,
    step: float,
    periodic: tuple[bool, bool],
    ghost_rows: dict[str, np.ndarray] | None = None,
    boundary_target: np.ndarray | None = None,
    params: HyperbolicParams = DEFAULT,
) -> np.ndarray
```
<a class="api-source" href="https://github.com/szaghi/shore/blob/main/src/shore/volume/hyperbolic_scheme.py#L144">source</a>

Advance *layer* (ni, nj, 3) by one hyperbolic step of size *step*.

**Parameters**

- **`layer`** — Current layer.
- **`step`** — Marching distance of this layer.
- **`periodic`** — `(periodic_i, periodic_j)`.
- **`ghost_rows`** — Neighbour rows across block seams (see `_compute_normals`).
- **`boundary_target`** — Full (ni, nj, 3) array whose non-periodic edge rows are the required new boundary positions (Dirichlet-pinned cap edges); interior values are ignored.
