# shore cap-tip

Build a wing-tip cap from a C-grid body block's `.geo` file.

The cap ([ADR 0059](/adr/0059-tip-cap-on-the-tip-face)) is an overset
component of five blocks on the wing's flat tip face: `cap_core` inside
the section, and the ring `cap_upper`, `cap_le`, `cap_lower` and `cap_te`
joining it to the wall.
- Their `k_lo` faces are the wall and cover the section exactly.
- k grows outwards, along the tip face's normal.
- The faces between cap blocks are `SHARED`; every other face is `FREE`, a
  fringe for the overset assembly with the C-grid wrap
  ([`shore cc-par`](./cc-par)).

See [C-grid topology — Tip caps](/concepts/topologies/c-grid#tip-caps) for
the layout. `shore mesh --cap-tips` builds the same caps together with the
C-grid.

## Synopsis

```
shore cap-tip [OPTIONS] BODY_GEO
```

## Arguments

| Argument | Description |
|---|---|
| `BODY_GEO` | Path to the C-grid body block's `.geo` file (typically `<stem>_body.geo` from `shore mesh --topology cgrid`). |

## Options

| Option | Default | Description |
|---|---|---|
| `--tip {jlo,jhi}` | `jlo` | Which spanwise tip to cap. `jlo` = j=0 station, `jhi` = j=n_stations-1. |
| `--le-frac FLOAT` | `0.05` | Chord fraction of the LE cuts, which bound `cap_le`. |
| `--te-frac FLOAT` | — | Deprecated, ignored since ADR 0059: `cap_te` spans the TE base. |
| `--tip-inset FLOAT` | auto | Distance of the cap's wall from the body's tip station, outwards. Default: the tip face of a C-grid made by `shore mesh` without caps, whose stations sit 2 % of the span inside the tips (`0.02 / 0.96` × the distance between the first and last stations). |
| `--core-frac FLOAT` | `0.5` | The core's half-thickness as a fraction of the section's, in (0, 1). |
| `--flare-deg FLOAT` | `30` | Outward lean of the ring's k-lines at the section's edge, in [0, 90) degrees. |
| `--n-cap-j INT` | `12` | Ring blocks' node count from the core to the wall. |
| `--nk INT` | auto | Cap k-axis node count. Default: the body block's `nk`. |
| `--cap-j-spacing {uniform,tanh,tanh2}` | `tanh` | Spacing law from the core to the wall (`tanh` clusters at the wall). |
| `--cap-j-beta FLOAT` | `4.0` | Clustering strength for tanh / tanh2. |
| `--te-thickness-frac FLOAT` | `0.002` | Rounding of a sharp TE, as a fraction of the chord (0.2 %). Applied in the cap only; the body wrap is left untouched. |
| `--te-thickness FLOAT` | — | Absolute rounding thickness (mutually exclusive with `--te-thickness-frac`). |
| `--no-auto-round-te` | off | Reject sharp-TE inputs instead of rounding them. |
| `--span-axis {x,y,z}` | `z` | World axis the wing spans along. Used to find the LE cuts. |
| `-o` / `--output PATH` | `cap_<tip>` | Output stem. Writes `<stem>_<label>.geo` per cap block. |
| `--adjacency-json PATH` | — | Also write a `.adjacency.json` sidecar (for `shore cc-par` consumption). |
| `--journal PATH` | — | Also write the run as a replayable journal ([`shore run`](run)). |

## Output

| File | Covers |
|---|---|
| `<stem>_cap_core.geo` | The inside of the section, inset from both sides. |
| `<stem>_cap_upper.geo` | From the core to the body's i = 0 side, between the LE cut and the TE. |
| `<stem>_cap_le.geo` | From the core round the leading edge, between the LE cuts. |
| `<stem>_cap_lower.geo` | From the core to the other side. |
| `<stem>_cap_te.geo` | From the core across the TE base. |

The optional `--adjacency-json` sidecar holds the 16 SHARED faces between
the cap blocks.

## Behaviour

### The section and the wall

`BODY_GEO` is a `(ni, n_stations, nk, 3)` block from `shore mesh --topology
cgrid`. The cap takes the wall nodes of its tip station (k = 0):

- `--tip jlo` → `body[:, 0, 0, :]`
- `--tip jhi` → `body[:, n_stations-1, 0, :]`

The next station inward gives the outward side of the section's plane. The
wall is that section moved out by `--tip-inset` along the plane's normal;
the cap's edge nodes are the section's nodes, moved. The layers follow the
body's k-lines at the tip station: their mean arc lengths, resampled when
`--nk` differs from the body's.

### LE cuts

The LE cut is found by walking from the LE point along the i = 0 side
until the chord position reaches `--le-frac`. The other cut is **snapped
symmetric about the section's middle node** (`lower_idx = ni - 1 -
upper_idx`), so the `upper` and `lower` blocks, and the core between them,
have the same node count.

### Auto-round on sharp TE

`cap_te` spans the TE base, so a sharp TE (i = 0 and i = ni − 1 of the
section coincide) is **rounded by default**, in the cap only: its node is
split into two points `te_thickness` apart. A one-line note goes to
stderr.

```bash
# Round to 0.5% chord
shore cap-tip wing_body.geo --tip jlo --te-thickness-frac 0.005

# Or pass an absolute thickness in body units
shore cap-tip wing_body.geo --tip jlo --te-thickness 0.001

# Reject sharp inputs entirely
shore cap-tip wing_body.geo --tip jlo --no-auto-round-te
```

A cap whose cells fold, or a block without cells, is refused with an error
naming the block ([ADR 0058](/adr/0058-no-silent-tip-cap-folds)).

## Examples

```bash
# Build the C-grid first
shore mesh naca0012.stl \
    --topology cgrid --ni-body 120 --ni-wake 30 --n-stations 5 \
    --wake-length 15 --nk 20 --ds 5e-4 --growth 1.15 -o wing

# Cap both tips
shore cap-tip wing_body.geo --tip jlo -o cap_jlo --adjacency-json cap_jlo.adjacency.json
shore cap-tip wing_body.geo --tip jhi -o cap_jhi --adjacency-json cap_jhi.adjacency.json

# Bundle for ParaView (13 blocks: 3 C-grid + 5 cap_jlo + 5 cap_jhi)
shore export wing_*.geo cap_jlo_cap_*.geo cap_jhi_cap_*.geo -o wing_capped.vtm
```

The full pipeline lives in
[`examples/06-c-grid/visualize_naca0012_cap_tip.sh`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012_cap_tip.sh).

## Python API

`Project.build_tip_cap` runs the same command. For finer control, drive
the topology directly with [`CGridTipCap`](/reference/api/shore.mesh.topology)
and [`Mesh.from_array`](/reference/api/shore.mesh.mesh):

```python
from shore.io.geo import read_geo
from shore.mesh import CGridTipCap, Mesh
from shore.topology.cgrid_cap import find_symmetric_le_cuts

body = read_geo("wing_body.geo")
cross = body[:, 0, 0, :]                       # jlo tip, k=0 (the wall)
le_up, le_lo = find_symmetric_le_cuts(cross, le_frac=0.05, span_axis_world=2)

topo = CGridTipCap(
    le_cut_upper_idx=le_up,
    le_cut_lower_idx=le_lo,
    span_dir=(0.0, 0.0, -1.0),     # outward at the jlo tip
    tip_inset=0.01,                 # station -> tip face
    n_cap_j=12,
)
cap = Mesh.from_array(surface=cross, topology=topo, nk=body.shape[2])
cap.write_geo("cap_jlo")            # 5 .geo files
```

## See also

- [C-grid topology — Tip caps](/concepts/topologies/c-grid#tip-caps)
- [ADR 0059](/adr/0059-tip-cap-on-the-tip-face): the cap's design.
- [`CGridTipCap`](/reference/api/shore.mesh.topology): Python API.
- [`Mesh.from_array`](/reference/api/shore.mesh.mesh): non-STL `Mesh` constructor.
- [`shore cc-par`](./cc-par): chimera assembly that composes the cap
  against the C-grid body wrap.
- [`examples/06-c-grid/`](https://github.com/szaghi/shore/tree/main/examples/06-c-grid).
