# Hyperbolic marching

Hyperbolic marching (also called *normal extrusion* or *advancing
layers*) extrudes the k=0 surface outward layer by layer, with each
layer advanced along the local outward normal direction. The
governing equations are hyperbolic in the wall-normal direction —
hence the name.

In SHORE the marcher operates on all 6 blocks of the caps topology
**simultaneously** with seam-aware normals at multi-block boundaries
so that block-to-block spacing continuity is preserved at every layer.

![Marching layers](/figures/marching-layers.png)

## Governing equation

Given the k=0 layer $\mathbf{P}^0$ and a sequence of step sizes
$\Delta s_k$, each subsequent layer is

$$
\mathbf{P}^{k+1}_{i,j} = \mathbf{P}^{k}_{i,j} + \Delta s_k \, \hat{\mathbf{n}}^{k}_{i,j},
$$

where $\hat{\mathbf{n}}^{k}_{i,j}$ is the outward unit normal at node
$(i, j)$ on layer $k$, computed from the in-plane neighbours.

## Per-layer recipe

For each layer index `k = 0, 1, ..., nk − 2`, SHORE does:

1. **Snapshot** the k-layer of all 6 blocks (a consistent state for
   ghost-row construction across blocks).
2. **Compute outward normals** for every node of every block. Boundary
   tangents use centred differences against
   [ghost rows](seam-aware-normals) supplied by the partner blocks; no
   one-sided fallback at multi-block seams.
3. **Step outward**: $\mathbf{P}^{k+1} = \mathbf{P}^k + \Delta s_k \,
   \hat{\mathbf{n}}^k$.
4. **Laplacian smoothing**: a few sweeps of in-plane averaging:

   $$
   \mathbf{P}_{i,j} \leftarrow \mathbf{P}_{i,j} + \alpha
     \left( \bar{\mathbf{P}}_{i,j} - \mathbf{P}_{i,j} \right),
   $$

   with $\bar{\mathbf{P}}_{i,j}$ the 4-neighbour average and $\alpha \in
   [0, 1]$ the blend factor (`--smooth`, default 0.2). The k=0 wall row
   and non-periodic boundary rows are clamped.
5. **Cap re-pinning**: each cap's 4 boundary edges are overwritten by
   the equatorial sub-block pole rows for layer `k+1`
   (`_overwrite_cap_boundaries`). This keeps cap-equator conformity
   bit-exact at every layer, regardless of what the cap's interior
   normal march would have produced.
6. **Jacobian guard**: compute the minimum cell Jacobian between layers
   `k` and `k+1`. If it falls below the tolerance the marcher raises
   `MarchingError` immediately — no silent bad meshes.
7. **Sub-sub stitching**: average the SHARED j-columns between adjacent
   sub-blocks to absorb FP-noise from the Laplacian smoother. With
   seam-aware normals this should be sub-FP-noise; the averaging is a
   belt-and-braces guard.
8. **Commit** all 6 blocks to layer `k+1`.

## Step sizes (`spacing_k`)

The per-layer step sizes $\Delta s_k$ are precomputed from a
`Spacing1D`:

| Law | Behaviour | Required knobs |
|---|---|---|
| `geometric` | $\Delta s_k = \Delta s_0 \cdot r^k$. Free far-field. | `ds`, `growth` |
| `tanh` | One-sided Roberts clustering toward the wall. Fixed total thickness. | `total_thickness`, `beta` |
| `tanh2` | Two-sided Roberts clustering at wall and far-field. | `total_thickness`, `beta` |

See [Spacing laws](spacing-laws) for the closed-form expressions.

## Normal estimation

For an interior node, both tangents are centred differences:

$$
\mathbf{t}_i = \frac{\mathbf{P}_{i+1,j} - \mathbf{P}_{i-1,j}}{2}, \qquad
\mathbf{t}_j = \frac{\mathbf{P}_{i,j+1} - \mathbf{P}_{i,j-1}}{2},
$$
$$
\hat{\mathbf{n}}_{i,j} = \frac{\mathbf{t}_i \times \mathbf{t}_j}{\|\mathbf{t}_i \times \mathbf{t}_j\|}.
$$

In the legacy single-block O-grid the j-direction is periodic, so the
boundary at `j = 0` wraps to `j = nj − 1`. In the caps topology both i
and j are non-periodic — the boundary tangent uses centred differences
against ghost rows from the partner block, see
[Seam-aware normals](seam-aware-normals).

## Laplacian smoothing

After each step, the new layer is relaxed by `smooth_iters` Jacobi
sweeps (default 2) of in-plane averaging. The k=0 wall row is always
clamped; non-periodic boundary rows/columns are clamped too. The blend
strength `α` (`--smooth`, default 0.2) trades smoothness against cell-
size uniformity.

Smoothing suppresses high-frequency wiggles that accumulate when the
surface has local curvature gradients. It is **not** a substitute for
elliptic volume smoothing — for production-quality grids, a post-
processing elliptic smoother (Sorenson, Thompson–Thames–Mastin) can
be applied after marching.

## Jacobian check

The cell Jacobian between layers `k` and `k+1` is the scalar triple
product of three edge vectors:

$$
J_{i,j} = (\mathbf{e}_i \times \mathbf{e}_j) \cdot \mathbf{e}_k,
$$

where $\mathbf{e}_i = \mathbf{P}_{i+1,j} - \mathbf{P}_{i,j}$,
$\mathbf{e}_j = \mathbf{P}_{i,j+1} - \mathbf{P}_{i,j}$, and
$\mathbf{e}_k = \mathbf{P}^{k+1}_{i,j} - \mathbf{P}^{k}_{i,j}$. A
positive value indicates a valid cell; a negative or zero value means
the cell has inverted.

In the caps topology all 6 blocks have open boundaries, so the marcher
uses `_jacobian_field_open` (no periodic wrap). The four **corner
cells** of each cap (where 4 sub-blocks meet at one node) are excluded
from the check: those cells have near-zero signed volume by
construction. The tolerance for the remaining cells is `-1e-6` to
absorb floating-point noise from the centred-difference normals.

## Optional: STL-normal blending

For sharp-feature bodies the grid-computed normals near the body lag
the true STL surface normals (the structured grid only "sees" body
geometry through its k=0 layer). The `--blend-normals-k` flag blends
grid normals with STL surface normals queried via
`trimesh.proximity.closest_point` for the first `K` layers.

The STL normal at each marching node is the **barycentric blend of the
closest triangle's three vertex normals at the closest point** — C0-
continuous across triangle edges, which suppresses the facet-aligned
banding that flat face normals would produce on a coarse STL. trimesh's
`vertex_normals` are area-weighted across incident faces and require a
clean mesh (consistent face orientation, no degenerate faces, merged
vertices) — all three are enforced unconditionally by `load_stl`'s
hygiene step. The legacy flat-face behaviour is still available via the
`interpolate=False` kwarg of `_stl_normals` for regression debugging.

$$
\hat{\mathbf{n}}_{\text{blend}} = \alpha \, \hat{\mathbf{n}}_{\text{STL}} + (1 - \alpha)\,\hat{\mathbf{n}}_{\text{grid}},
\quad \alpha = 1 - \frac{k}{K}.
$$

At `k = 0` the layer extrudes along pure STL normals; the blend
linearly fades to pure grid normals at `k = K`.

::: info Interaction with the seam-aware kernel
The grid normal $\hat{\mathbf{n}}_{\text{grid}}$ used in the blend is
already the **seam-aware ghost-row centred-difference** normal — the
blend is applied on top of that, not on top of the one-sided
fallback.  With **barycentric-interpolated** STL normals the blended
normal field is C0-continuous across triangle edges, and the
alpha-weighted contribution adds no measurable kick to the cap-equator
C1 seam: the blend run tracks the no-blend run within ~0.13 deg per
layer on the icosphere fixture (vs. ~0.55 deg with the old flat
face-normal implementation).

A previous implementation also re-extruded each sub-block's
cap-equator seam rows along **pure STL normals every layer**
("pole-row override").  The pole-row override pre-dated the
seam-aware kernel and was needed to compensate for the Laplacian
boundary clamp pulling the seam rows inward; with the seam-aware
kernel in place that clamp produces only sub-FP-noise drift and the
override itself amplified the C1 seam deviation by a factor of ~2 at
`k = NK − 1`.  The override has been removed and the regression is
locked down by
[`TestBlendNormalsK_C1SeamSurvives`](https://github.com/szaghi/shore/blob/main/tests/test_cap_gnomonic.py).
:::

## Failure modes and remedies

| Symptom | Cause | Remedy |
|---|---|---|
| `MarchingError` at small k | High curvature; normals converging | Reduce `ds`; increase `smooth_strength` |
| jmin → 0 at large k | Outer layers collapsing | Reduce `nk` or `growth` |
| Visible kink at multi-block seam | Seam-aware normals disabled | Should not happen — seam-aware normals are always on for caps |
| Jagged grid near sharp body features | Grid normals not tracking body | Use `--blend-normals-k 5` to blend STL normals near the wall |

## Comparison to elliptic methods

| Property | Hyperbolic marching | Elliptic (Thompson–Thames–Mastin) |
|---|---|---|
| Direction | Wall-normal (one sweep) | All directions (iterative) |
| Near-wall orthogonality | Good by construction | Controlled via source terms |
| Far-field smoothness | Bounded by step law | Excellent |
| Computational cost | $O(n_i \, n_j \, n_k)$ | $O(n_i \, n_j \, n_k \cdot \text{iters})$ |
| Failure mode | Inversion at concave surfaces | Slow convergence; fold-over |

For Chimera near-body grids, hyperbolic marching is the preferred
choice: the near-wall region is the critical one for the solver,
far-field smoothness is handled by the background grid, and the
single-sweep cost is orders of magnitude lower than iterative methods.

## References

- Steger, J. L. & Sorenson, R. L. (1979). *Automatic mesh-point
  clustering near a boundary in grid generation with elliptic partial
  differential equations.* J. Comput. Phys., 33(3), 405–410.
- Chan, W. M. & Steger, J. L. (1992). *Enhancements of a three-
  dimensional hyperbolic grid generation scheme.* Appl. Math. Comput.,
  51(2–3), 181–205.
- Blazek, J. (2015). *Computational Fluid Dynamics: Principles and
  Applications*, Chapter 9.

## See also

- [Seam-aware normals](seam-aware-normals) — centred differences across
  multi-block seams
- [Spacing laws](spacing-laws) — `Δs_k` formulas
- [`shore.volume.hyperbolic`](/reference/api-hyperbolic) and
  [`shore.topology.sphere_cap`](/reference/api-sphere-cap) (Python API)
- [Caps topology](/topologies/caps) — block layout + connectivity
