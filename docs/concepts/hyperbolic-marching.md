# Volume marching

Marching extrudes the k=0 wall grid outward one layer at a time until the
block has `nk` layers. SHORE has two schemes, selected by `march_scheme`
([ADR 0005](/adr/0005-marching-schemes)):

| Scheme | What it solves | Default |
|---|---|---|
| `algebraic` | normal extrusion followed by Laplacian smoothing | yes |
| `hyperbolic` | the Steger–Chaussee system: k-lines orthogonal to each layer, prescribed cell volume | no |

The algebraic scheme is what older SHORE documentation called
"hyperbolic marching". It does not solve a hyperbolic system; ADR 0005
renamed it and added the real one. Both schemes share the step schedule,
the seam machinery and the Jacobian guard.

Select the scheme with `--march-scheme` on `shore mesh`, `[march]
march_scheme` in a config, or `march_scheme=` on `Mesh.march` /
`Project.march`.

<Figure name="marching-layers" caption="A marched cubed-sphere volume with a quarter cut away: k-layers coloured by depth from the wall." />

## Step schedule

The marching distances $\Delta s_k$, $k = 0 \dots n_k - 2$, come from the
wall-normal [spacing law](/concepts/spacing-laws):

| Law | $\Delta s_k$ | Knobs |
|---|---|---|
| `geometric` (default) | $\Delta s_0\, r^k$, outer distance free | `--ds` (0.001), `--growth` (1.15) |
| `tanh` | one-sided Roberts, clustered at the wall, total thickness fixed | `--volume-k-thickness`, `--volume-k-beta` (3.0) |
| `tanh2` | two-sided Roberts, clustered at wall and outer layer | `--volume-k-thickness`, `--volume-k-beta` (3.0) |

All blocks of a mesh march with **one** schedule, because blocks coupled
through k-seams must stay coincident layer by layer. A k-edge override
(a `Spacing1D` on the k-edges, detected by identity) on block 0 wins; other
blocks that disagree get a `UserWarning`.

## Algebraic scheme

Each layer is

$$
\mathbf{P}^{k+1}_{i,j} = \mathbf{P}^{k}_{i,j} + \Delta s_k\, \hat{\mathbf{n}}^{k}_{i,j},
\qquad
\hat{\mathbf{n}} = \frac{\mathbf{t}_i \times \mathbf{t}_j}{\lVert \mathbf{t}_i \times \mathbf{t}_j \rVert},
$$

with the tangents $\mathbf{t}_i$, $\mathbf{t}_j$ centred differences of
the layer. On a block edge with a partner block the tangent is a centred
difference against a ghost row from the partner
([Seam-aware normals](/concepts/seam-aware-normals)); a periodic direction
wraps; only an unconnected edge falls back to a one-sided difference.

The new layer is then relaxed by `smooth_iters` sweeps (`--smooth-iters`,
default 2) of

$$
\mathbf{P}_{i,j} \leftarrow \mathbf{P}_{i,j} + \alpha\,(\bar{\mathbf{P}}_{i,j} - \mathbf{P}_{i,j}),
$$

$\bar{\mathbf{P}}$ the 4-neighbour average and $\alpha$ the strength
(`--smooth`, default 0.2). The wall and the non-periodic boundary rows are
clamped.

### Smoothing limiter (`smooth_mode`)

Smoothing node positions on a curved layer also pulls the layer toward
its centre of curvature, by about $\alpha\, h^2 \kappa$ ($h$ the cell
size, $\kappa$ the curvature), independently of the step. When that
exceeds the step, the new layer falls behind the old one and the cells
fold: on a coarse sphere a first step of 1e-3 folded at k = 1
([ADR 0005](/adr/0005-marching-schemes), "Finding").

`smooth_mode='limited'` (default; `--smooth-mode`, `[march] smooth_mode`)
splits each node's smoothing motion into its normal and tangential parts
and limits them:

- the pull-back along the marching normal is at most $0.25\,\Delta s_k$,
  so every layer advances at least $0.75\,\Delta s_k$;
- the sideways motion is at most $0.25\,\Delta s_k$.

Nodes inside both limits keep the plain Laplacian result bit for bit.
(0.25, 0.25) is the tightest pair that leaves the 60×60×10 cubed-sphere
regression grid unchanged. `smooth_mode='legacy'` switches the limiter
off to reproduce grids made before it; it folds thin first layers.

### Layer sequence on the cubed sphere

`Mesh._march_caps` advances the six blocks together:

1. Snapshot layer k of all six blocks; ghost rows are taken from the
   snapshot, so no block reads a neighbour that has already moved.
2. Advance the four sub-blocks (normals with ghost rows, step, smoothing).
3. Average each shared sub–sub column, to remove round-off differences
   left by smoothing.
4. Check the Jacobian of the sub-blocks.
5. Advance the two caps, then overwrite their four edges with the new
   pole rows of the sub-blocks (`DIRICHLET`), and check their Jacobian.
6. Commit layer k+1.

The O-grid, C-grid and hybrid wraps follow the same pattern with their
own ghost rows (`_march_ogrid`, `_march_cgrid_three`).

## Hyperbolic scheme

The new layer $\mathbf{r}_{k+1} = \mathbf{r}_k + \Delta\mathbf{r}$
satisfies, at every node,

$$
\mathbf{r}_\xi \cdot \mathbf{r}_\zeta = 0, \qquad
\mathbf{r}_\eta \cdot \mathbf{r}_\zeta = 0, \qquad
\det[\mathbf{r}_\xi, \mathbf{r}_\eta, \mathbf{r}_\zeta] = \Delta V ,
$$

k-lines orthogonal to the layer and a prescribed cell volume (Steger &
Chaussee 1980; Chan & Steger 1992). Linearised about the local solution
$\mathbf{r}_\zeta^0 = \Delta V\,\mathbf{n}/\lVert\mathbf{n}\rVert^2$,
$\mathbf{n} = \mathbf{r}_\xi \times \mathbf{r}_\eta$, it becomes one
sparse linear system per block layer,

$$
\left(I + C^{-1}A\,\delta_\xi + C^{-1}B\,\delta_\eta - \varepsilon_i \nabla^2\right)\Delta\mathbf{r}
= \mathbf{r}_\zeta^0 + \varepsilon_e\, T\!\left(\nabla^2 \mathbf{r}_k\right),
$$

with central differences $\delta$, the 5-point Laplacian, and $T$ the
projection onto the layer's tangent plane. It is solved with SuperLU
(deterministic). The coupling $C^{-1}A$, $C^{-1}B$ scales with step over
cell size: thin near-wall layers march almost exactly along the local
normal; thick layers feel their neighbours, which spreads converging
fronts at concave corners.

Stabilisation, all ramped by $\rho = \min(1, \Delta s_k / h)$:

| Term | Default (`HyperbolicParams`) |
|---|---|
| implicit dissipation $\varepsilon_i$ | 0.5 |
| explicit tangential dissipation $\varepsilon_e$ (not applied across seams) | 0.25 |
| cell-volume smoothing toward the neighbour mean | 0.2, 2 iterations |
| front-convergence sensor gain | 16 |
| marching distance kept in $[\Delta s/2,\ 2\Delta s]$ | clip 2 |
| outer-face volume correction, at most this fraction of $\Delta V$ | 0.04 |

The defaults were chosen once over all spec cases, not per case. The
`smooth_*` arguments are not used by this scheme. Edges with a partner
block take ghost rows as unknowns; unconnected edges and the pinned cap
edges are Dirichlet.

### What the acceptance spec measured

From the frozen spec of ADR 0005 (13 passed, 0 skipped, 0 xfailed):

| Case | Algebraic (limited) | Hyperbolic |
|---|---|---|
| 60×60×10 cubed sphere: wall-angle mean / max | 2.05 / 4.6 deg | 0.96 / 2.33 deg |
| same: min scaled Jacobian | | 0.740 |
| plus-shaped prism, four concave 90° edges, nk = 15 | folds at k = 13 | completes, no folded cell |

Unverified, as the spec states: blunt-TE C-grids, CH/OH channels,
non-star-shaped CAD bodies, finer grids and run time.

## Jacobian guard

After every layer each block's cell Jacobians between layers k and k+1
are checked,

$$
J = (\mathbf{e}_i \times \mathbf{e}_j)\cdot\mathbf{e}_k ,
$$

on the open-block field (`_jacobian_field_open`) or the j-periodic one.
On the cubed sphere the four corner cells of each block are excluded (the
3-block junctions there give near-zero volumes) and the tolerance is
$-10^{-6}$. A failing layer raises `MarchingError` at once: SHORE never
writes a folded grid.

## Removed: normal blending

Earlier releases could blend body normals into the first `K` marching
directions (`blend_normals_k`). It showed no benefit over the seam-aware
grid normals and folded C-grids at the wake cut
([ADR 0019](/adr/0019-deprecate-normal-blending)); it is removed
([ADR 0034](/adr/0034-remove-normal-blending)). An explicit 0 is accepted
with a `DeprecationWarning`; any other value raises `ParameterError`.

## Failure modes

| Symptom | Likely cause | Remedy |
|---|---|---|
| `MarchingError` at k = 1 | smoothing pulls a thin layer behind the wall | keep `smooth_mode='limited'`, or lower `--smooth` |
| `MarchingError` at larger k near concave edges | converging fronts | `--march-scheme hyperbolic`, or a smaller `nk` / `growth` |
| kink along a block seam | not expected: seam-aware normals are always on | report it |
| first cell thinner than `ds` | smoothing pull-back (algebraic) | `--march-scheme hyperbolic` keeps the first cell within a few per cent of `ds` |

## Choosing a scheme

The algebraic scheme is the default because its grids are bit-identical
to earlier releases and it is fast. Use the hyperbolic scheme when wall
orthogonality matters, for concave walls, or when the first cell must
match `ds` closely. Both are single sweeps outward: the far-field
smoothness of an elliptic generator is not their goal, and in an overset
assembly the background grid takes over there.

## References

- Steger, J. L. & Chaussee, D. S. (1980). *Generation of body-fitted
  coordinates using hyperbolic partial differential equations.* SIAM J.
  Sci. Stat. Comput., 1(4), 431–437.
- Chan, W. M. & Steger, J. L. (1992). *Enhancements of a three-dimensional
  hyperbolic grid generation scheme.* Appl. Math. Comput., 51(2–3),
  181–205.

## See also

- [Seam-aware normals](/concepts/seam-aware-normals): ghost rows at block seams
- [Spacing laws](/concepts/spacing-laws): the $\Delta s_k$ formulas
- [Grid quality](/concepts/quality): the metrics the spec is written in
- [`Mesh.march`](/reference/api/shore.mesh.mesh#Mesh-march), [`shore.volume.hyperbolic_scheme`](/reference/api/shore.volume.hyperbolic_scheme), [`shore.topology.sphere_cap`](/reference/api/shore.topology.sphere_cap)
- [`shore mesh`](/reference/mesh), [Set the wall spacing](/cookbook/wall-spacing)
