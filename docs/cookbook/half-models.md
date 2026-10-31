# Mesh a half model on a symmetry plane

Use this when the configuration is mirror-symmetric and you grid only one half:
the bodies are cut at the plane, the junction grid ends on it, and its faces on
the plane carry a SYMMETRY boundary condition through to `cc.par`.

The script runs in about 2 s.

## Cut a body into a half body

<<< @/../examples/cookbook/half_models.py#half_bodies{python}

- A half body is an open STL: the facets on the plane are removed, so the
  surface ends on it. Here a box (the "wing") and a vertical post, the half
  model of [ADR 0011](/adr/0011-open-intersection-curves), cut at y = 0.
- An open STL is not watertight, so loading it warns. The script loads it with
  `skip_hygiene=True`, as the half-collar acceptance tests do: the facets and
  their orientation stay as written (no vertex merge, degenerate-face removal
  or `fix_normals`). The CLI runs the hygiene pass; on these bodies both give
  the same grid.
- The `cubed_sphere`, `ogrid` and `oh` topologies ray-cast from an anchor
  inside a closed, star-shaped body and do not take half bodies. Half models
  serve the junction grids: intersection curves, surface marching and collars.

## Build a half-collar at the junction

::: code-group
<<< @/../examples/cookbook/half_models.py#half_collar{python} [Python]

```bash [CLI]
shore collar half_box.stl half_post.stl --seed 0.45 0.2 0.02 --n-loop 25 \
    --ds-surface 0.01 --growth-surface 1.12 --m-a 12 --m-b 12 \
    --ds-wall 1e-4 --growth-wall 1.2 --nk 25 --plane 0 0 0 0 1 0 \
    -o half_collar --adjacency-json half_collar.adjacency.json
```
:::

```text
nodes: (25, 25, 26, 3)
j-faces off the plane by 0.0
{'i_lo': 'FREE', 'i_hi': 'FREE', 'j_lo': 'SYMMETRY', 'j_hi': 'SYMMETRY', 'k_lo': 'WALL', 'k_hi': 'FREE'}
```

- `plane=(point, normal)` (CLI `--plane PX PY PZ NX NY NZ`) turns the collar
  into a half-collar: the junction curve is traced from the seed both ways
  until it reaches the plane, the surface grids march with both ends held on
  the plane, and the block is open in j with `j_lo` and `j_hi` on the plane.
- Block shape: `m_b + m_a + 1` nodes across the junction, `n_loop` along it,
  `nk + 1` wall-normal.
- On mirror-symmetric bodies the half-collar is the y >= 0 half of the whole
  collar built on the mirrored wall, and its j-lines meet the plane at a right
  angle.
- When the curve meets the plane obliquely (more than 1 deg off the plane
  normal) the build warns: the faces are still on the plane, but the quality
  near it is not guaranteed.
- The CLI prints the two not-watertight warnings of the half bodies; they are
  expected. Add `--project DIR` to record the build in a project.
- One plane only: a model cut by two planes is out of scope.

## Write the SYMMETRY faces to `cc.par`

::: code-group
<<< @/../examples/cookbook/half_models.py#symmetry_bc{python} [Python]

```bash [CLI]
shore cc-par half_collar.adjacency.json cc.par --grd half_collar.grd
```
:::

The two j-faces come out as natural-BC code `2` (symmetry), the wall as `1`,
the free faces as the chimera code `40`:

```text
        1         3         2         0         0        24         0         0         0        25 !         3
        1         4         2         0         0        24        24        24         0        25 !         4
        1         5         1         0         0        24         0        24         0         0 !         5
```

`FaceBC.SYMMETRY` round-trips through `.adjacency.json` and is written as
`BCSymmetryPlane` by the CGNS export ([Export CGNS and Plot3D](/cookbook/exports)).

## Declare the symmetry plane of a blocking

A blocked half model declares its cut plane once, after `create_blocking`:

```python
p.create_blocking(53, 33, 25, vmin=(-6, -4, 0), vmax=(7, 4, 6))
p.symmetry_plane((0, 0, 0), (0, 0, 1))
```

`premesh` then makes the faces on the plane SYMMETRY faces, and
`orthogonalize` keeps the wall-normal lines that start on it in it
([ADR 0053](/adr/0053-blocking-symmetry-planes)); without the declaration
those faces are FREE and drift off the plane. [Blocking a wing–body](/tutorials/blocking)
uses it.

## See also

- [Intersection curves](/concepts/intersection-curves),
  [Surface marching](/concepts/surface-marching), [Collar blocks](/concepts/collars)
- [`shore collar`](/reference/collar), [`shore cc-par`](/reference/cc-par)
- [`shore.collar`](/reference/api/shore.collar)
- [ADR 0011](/adr/0011-open-intersection-curves), [ADR 0012](/adr/0012-march-ends-on-a-plane),
  [ADR 0013](/adr/0013-half-collar-block)
