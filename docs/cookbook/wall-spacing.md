# Set the wall-normal spacing

The wall-normal (k) distribution sets the first cell at the wall, the growth
through the boundary layer and how far the marched grid reaches. These
recipes size it with the geometric law or a tanh law, pick the first cell
from a target y+, and explain how the schedule is shared between blocks.

The code is
[`examples/cookbook/wall_spacing.py`](https://github.com/szaghi/shore/blob/main/examples/cookbook/wall_spacing.py)
(about 3 s).

<Figure name="spacing-laws" ext="svg" caption="Node distributions of the uniform, tanh, tanh2 and geometric laws on one axis, with their cell sizes." />

## Reach a thickness with geometric layers

The default law: the first step is `ds`, each next one `growth` times the
previous, `ds * growth**k` for `k = 0 .. nk-2`. The thickness follows from
`nk`, so solve for `nk`:

<<< @/../examples/cookbook/wall_spacing.py#geometric{python}

```text
geometric: nk = 32, reach 0.501, last step 6.621e-02
```

::: code-group
```bash [CLI]
shore mesh sphere.stl --ni 33 --nj 64 --nk 41 --ds 4.7e-5 --growth 1.2 \
    --theta-cap-deg 35 -o cs
```
```python [Python]
p.build_topology("cubed_sphere", nk=41, ds=4.7e-5, growth=1.2, ni=33, nj=64,
                 theta_cap=np.deg2rad(35.0))
```
:::

- Defaults: `ds = 1e-3`, `growth = 1.15` (`--ds`, `--growth`); `ds` is in
  the units of the geometry.
- Keep `growth` at or below about 1.2 inside the boundary layer; the last
  step grows fast (`6.6e-2` above, 66 times `ds`), which is what you want
  for a near-body Chimera grid that ends in an overlap.

## Fix the thickness with a tanh law

`tanh` (one-sided Roberts) clusters at the wall and `tanh2` (two-sided) at
both ends, over a fixed `total_thickness`. They take no `ds`: the first step
follows from the clustering strength `beta` (default 3.0). To hit a target
first cell, solve for `beta`; the first step falls monotonically as `beta`
grows, so bisection works:

<<< @/../examples/cookbook/wall_spacing.py#tanh{python}

```text
tanh : beta = 2.582, first 1.00e-03, last 4.07e-02, max ratio 1.173
tanh2: beta = 5.272, first 1.00e-03, last 1.00e-03, max ratio 1.384
```

::: code-group
```bash [CLI]
shore mesh sphere.stl --topology ogrid --ni 33 --nj 49 --nk 33 --theta-cap-deg 10 \
    --volume-k-spacing tanh --volume-k-beta 2.582 --volume-k-thickness 0.5 -o og
```
```python [Python]
p.build_topology("ogrid", nk=33, k_law="tanh", k_beta=2.582, k_thickness=0.5,
                 ni=33, nj=49, theta_cap=np.deg2rad(10.0))
```
:::

- `tanh` and `tanh2` need the thickness (`--volume-k-thickness`,
  `k_thickness=`); without it the build raises `ParameterError`.
- For the same first cell `tanh2` needs a larger `beta` and its largest
  cell-to-cell ratio is higher (1.38 against 1.17 here): it spends nodes
  on the outer end. Use it when the outer boundary is a seam to another
  block that needs fine cells, not for a free outer boundary.
- With a fixed thickness and `nk`, a smaller first cell raises the growth
  ratio. Add nodes when the maximum ratio passes about 1.2.

## Choose the first cell from a target y+

SHORE has no flow model; the first cell is your input. A common starting
point is the turbulent flat-plate estimate below, evaluated at a reference
length `L` (body length or chord) and `Re = U L / ν`:

$$
C_f \approx 0.026\,Re^{-1/7},\qquad
u_\tau = U\sqrt{C_f/2},\qquad
\frac{\Delta s}{L} = \frac{y^+\,\nu}{u_\tau L} = \frac{y^+\sqrt{2/C_f}}{Re}.
$$

<<< @/../examples/cookbook/wall_spacing.py#yplus{python}

```text
Re = 1e+05: first cell for y+ = 1 is 2.00e-04 L
Re = 1e+06: first cell for y+ = 1 is 2.35e-05 L
Re = 1e+07: first cell for y+ = 1 is 2.77e-06 L
```

::: warning This is an estimate
The flat-plate correlation ignores pressure gradients, laminar regions,
separation and curvature. The local y+ of the converged solution varies over
the body, often by a factor of 2 or more. Check y+ in the flow solution and
adjust `ds`.
:::

Feed the estimate into the build; here for a sphere of diameter 2 at
`Re = 1e6` based on the diameter:

<<< @/../examples/cookbook/wall_spacing.py#march{python}

```text
ds = 4.71e-05: first cell 3.53e-05..4.71e-05
```

The measured first cell is between `0.75 ds` and `ds`: the default
`smooth_mode="limited"` lets the smoothing pull a node back by at most a
quarter of the step along the marching normal. Size `ds` for the smallest
first cell you accept, or march with `smooth_iters=0` where the exact value
matters more than the smoothing.

## Know which schedule a block uses

<<< @/../examples/cookbook/wall_spacing.py#components{python}

```text
main: 6 blocks, nk = 41, steps 5.65e-05..5.76e-02
ogrid: 1 blocks, nk = 33, steps 8.80e-04..3.72e-02
```

- All blocks of one component (one `Mesh`) march with **one** k-schedule:
  k-coupled seams must stay coincident layer by layer.
- Different components of a project are independent; each has its own law,
  `ds` and `nk`.
- If you override k-edge spacings per block, block 0's override wins and a
  `UserWarning` names any block that disagrees. The precedence is on
  [Per-edge spacing](/cookbook/per-edge-spacing#precedence-rules).
- The hybrid topologies (`oh`, `ch`) have a separate law for their channel
  blocks (`--nk-channel`, `--k-law-channel`, `--k-beta-channel`).

## See also

- [Spacing laws](/concepts/spacing-laws): the formulas.
- [Marching schemes](/concepts/hyperbolic-marching): smoothing and its
  limiter.
- [Per-edge spacing](/cookbook/per-edge-spacing) and
  [Redistribute an edge](/cookbook/redistribute): spacing along the
  other axes.
- [Cluster a wall layer](/cookbook/blocking#cluster-a-wall-layer-with-an-absolute-first-cell):
  the first cell on a premeshed blocking.
- CLI: [`shore mesh`](/reference/mesh).
- API: [`build_k_steps`](/reference/api/shore.volume.spacing#build_k_steps),
  [`Spacing1D`](/reference/api/shore.mesh.spacing#Spacing1D),
  [`Project.build_topology`](/reference/api/shore.project.project#Project-build_topology).
