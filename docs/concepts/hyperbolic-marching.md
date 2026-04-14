# Hyperbolic marching

Hyperbolic marching (also called *normal extrusion* or *advancing layers*) is a class of structured grid generation methods where a surface mesh is extruded layer by layer along the surface normal direction. The governing equations are hyperbolic in the wall-normal direction — hence the name.

## Governing equation

Given the k=0 layer $\mathbf{P}^0$, each subsequent layer is:

$$\mathbf{P}_{i,j}^{k+1} = \mathbf{P}_{i,j}^{k} + \Delta s_k \cdot \hat{\mathbf{n}}_{i,j}^{k}$$

where:

- $\Delta s_k = \Delta s_0 \cdot r^k$ is the layer thickness at step $k$ (geometric growth with ratio $r$)
- $\hat{\mathbf{n}}_{i,j}^{k}$ is the outward unit normal at node $(i, j)$ on layer $k$

## Normal estimation

At each node $(i, j, k)$ the normal is estimated from the four in-plane neighbours using the cross product of the two tangent vectors:

$$\hat{\mathbf{n}}_{i,j} = \frac{\mathbf{t}_i \times \mathbf{t}_j}{\|\mathbf{t}_i \times \mathbf{t}_j\|}$$

where:

$$\mathbf{t}_i = \mathbf{P}_{i+1,j} - \mathbf{P}_{i-1,j} \qquad \mathbf{t}_j = \mathbf{P}_{i,j+1} - \mathbf{P}_{i,j-1}$$

Central differences are used in the interior; one-sided differences at the i-boundaries (i=0 and i=ni-1). The j-direction is treated as **periodic** (O-grid topology), so $j=-1$ wraps to $j=n_j-1$.

## Laplacian smoothing

After each marching step, the new layer is smoothed by a single-pass Laplacian:

$$\mathbf{P}^{k+1,\text{smooth}}_{i,j} = \mathbf{P}^{k+1}_{i,j} + \alpha \left( \bar{\mathbf{P}}_{i,j} - \mathbf{P}^{k+1}_{i,j} \right)$$

where $\bar{\mathbf{P}}_{i,j}$ is the average of the four in-plane neighbours and $\alpha \in [0,1]$ is the blend strength (`--smooth` option, default 0.2). The k=0 wall nodes are never moved.

This suppresses high-frequency wiggles that accumulate when the surface has local curvature gradients. It does not preserve cell volume and is not a substitute for elliptic volume smoothing — for production-quality grids, a post-processing elliptic smoother (Sorenson, Thompson–Thames–Mastin) is recommended.

## Jacobian check

After smoothing, the minimum cell Jacobian is computed for every cell between layer $k$ and layer $k+1$:

$$J_{i,j} = (\mathbf{e}_2 \times \mathbf{e}_1) \cdot \mathbf{e}_3$$

where:

- $\mathbf{e}_1 = \mathbf{P}_{i,j+1}^k - \mathbf{P}_{i,j}^k$ (j-tangent, periodic)
- $\mathbf{e}_2 = \mathbf{P}_{i+1,j}^k - \mathbf{P}_{i,j}^k$ (i-tangent)
- $\mathbf{e}_3 = \mathbf{P}_{i,j}^{k+1} - \mathbf{P}_{i,j}^k$ (wall-normal step)

$J > 0$ for a valid cell; $J \le 0$ indicates an inverted cell. If any cell in a layer inverts, marching stops immediately with a `RuntimeError`.

::: tip Sign convention
The sign of $J$ depends on the parametrisation orientation. For a lat-lon O-grid the convention above (cross product of i-tangent and j-tangent pointing outward) gives $J > 0$ for correctly extruded cells. This is verified by the test `test_normals_sphere_radially_outward`.
:::

## Failure modes and remedies

| Symptom | Cause | Remedy |
|---|---|---|
| `RuntimeError: Negative Jacobian at layer k=1` | Very high surface curvature, normals immediately converging | Reduce `--ds`, increase `--smooth` |
| Jacobian drops near zero at high k | Outer layers collapsing — sphere-shaped outer boundary contracts at poles | Reduce `--nk`, or reduce `--growth` |
| Jagged grid near sharp edges | One-sided normal estimation at i-boundaries | Increase `--smooth-iters` |
| Slow march at large ni×nj | NumPy normal computation scales as O(ni·nj) | Typical values (ni=60, nj=80) run in < 1 s for nk=40 |

## Comparison to elliptic methods

| Property | Hyperbolic marching | Elliptic (Thompson–Thames–Mastin) |
|---|---|---|
| Direction | Wall-normal (one sweep) | All directions (iterative) |
| Near-wall orthogonality | Good (by construction) | Controlled via source terms |
| Far-field smoothness | Poor without morphing | Good |
| Computational cost | O(ni · nj · nk) | O(ni · nj · nk · iterations) |
| Failure mode | Inversion at concave surfaces | Slow convergence or fold-over |

For Chimera near-body grids, hyperbolic marching is the preferred method: the near-wall region is the critical one for the CFD solver, far-field smoothness is handled by the background grid, and the single-sweep cost is orders of magnitude lower than iterative elliptic methods.

## References

- Steger, J. L. & Sorenson, R. L. (1979). *Automatic mesh-point clustering near a boundary in grid generation with elliptic partial differential equations.* Journal of Computational Physics, 33(3), 405–410.
- Chan, W. M. & Steger, J. L. (1992). *Enhancements of a three-dimensional hyperbolic grid generation scheme.* Applied Mathematics and Computation, 51(2–3), 181–205.
- Blazek, J. (2015). *Computational Fluid Dynamics: Principles and Applications*, Chapter 9.
