# Redistribute the nodes along an edge

Use this when a built or marched grid has the right shape but the wrong node
spacing along one direction: you pick a block edge and a spacing law, and every
grid line of that block along the edge's axis is resampled, through the seams.

## Cluster the nodes along one block edge

<<< @/../examples/cookbook/redistribute.py#one_block{python}

```text
tanh  : [0.     0.011  0.0342 0.0824 0.1807 0.3734 0.7234 1.2797 2.    ]
pinned: [0.     0.05   0.3083 0.5667 0.825  1.0833 1.3417 1.6    2.    ] first 0.05 last 0.3999999999999999
```

- The edge names the axis: `j_lo_k_lo` fixes j and k, so it runs along i.
  The twelve edge names are listed on the
  [block model](/concepts/block-model) page.
- Laws: `uniform`, `tanh` (one-sided, `beta`), `tanh2` (two-sided, `beta`,
  `cluster_at` as a fraction of the edge), `pinned` (`first_step` and
  `last_step` as lengths along the picked edge, `interior` law between them).
  The law is evaluated over the picked edge's arc length and turned into
  arc-length fractions.
- Every other line along the axis takes the same fractions of its own arc
  length. Nodes stay on the line's old polyline and the two end nodes do not
  move, so the two faces across the axis are unchanged bit for bit.
- It is one command: `p.undo()` restores the nodes bit for bit, and the
  journal records `p.redistribute('box', 'j_lo_k_lo', 'tanh')` (defaults, here
  `beta=3.0`, are left out).

## Redistribute through seams on a marched grid

<<< @/../examples/cookbook/redistribute.py#through_seams{python}

```text
sub0       changed: True
sub1       changed: True
sub2       changed: True
sub3       changed: True
cap_north  changed: False
cap_south  changed: False
```

- The redistribution follows every SHARED or DIRICHLET face that contains the
  axis, transitively, and maps the axis across each seam (reversed across a
  flipped one, such as the C-grid wake cut). Seams stay closed. Here the
  i-direction of `sub0` reaches the other three equatorial blocks; the caps
  only touch the equator through faces across i, so they do not change.
- Associated nodes are not re-projected: a wall node resampled along a
  polyline lies on its chord, not on the body. Run `p.snap()` after (see
  [Associate with geometry](/cookbook/associate)).
- It works on a project loaded without replay too: there the seams are copies,
  not views, and the command writes them through
  ([ADR 0040](/adr/0040-seams-kept-closed)).
- `march` keeps its own k schedule; to change the wall-normal spacing before
  marching see [Wall-normal spacing](/cookbook/wall-spacing).

## Know what is refused

<<< @/../examples/cookbook/redistribute.py#refusals{python}

```text
refused: redistribute: 'k_lo' is not an edge of 'box'; edges: ['j_lo_k_lo', ...].
refused: redistribute: law must be one of ('uniform', 'tanh', 'tanh2', 'pinned'), got 'spline'.
refused: first_step + last_step (3) exceeds interval span (2)
```

Each refusal is a `ParameterError`, raised before anything moves: the nodes
and the journal are unchanged. Refused as well:

- a law that does not give increasing nodes over the edge;
- a seam face shared with several blocks on the way (the CH hybrid's
  `body.k_hi`, which faces several channel blocks);
- a PERIODIC face containing the axis;
- seams that loop back and would need two different distributions on one
  block.

## See also

- [Spacing laws](/concepts/spacing-laws) and [Per-edge spacing](/cookbook/per-edge-spacing)
- [Block model](/concepts/block-model): edges, faces, axis maps
- [`shore.mesh.redistribute`](/reference/api/shore.mesh.redistribute)
- [ADR 0041](/adr/0041-redistribute-edge), [ADR 0040](/adr/0040-seams-kept-closed)
