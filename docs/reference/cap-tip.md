# shore cap-tip

Build a wing-tip cap from a C-grid body block's `.geo` file.

The cap is an independent 3-block (sharp TE) or 4-block (blunt TE)
H-style structured component that closes the wing's airfoil cross-
section at one spanwise extreme.  It composes against the C-grid body
wrap as a **chimera component** via [`shore cc-par`](./cc-par) — the
cap's `k_lo` face is `WALL` (the actual wing-tip closure), `k_hi` is
`FREE` (chimera fringe receiving from the body wrap).

See [C-grid topology — Tip caps](/topologies/c-grid#tip-caps-closed-finite-wing)
for the geometric layout, the rationale for chimera-overlap (vs
SHARED-face) cap-body coupling, and the H3 T-fill block decomposition.

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
| `--le-frac FLOAT` | `0.05` | Chord fraction for the LE-cut. Must be in `(0, te_frac)`. |
| `--te-frac FLOAT` | `0.95` | Chord fraction for the TE-cut.  Splits the chord band into a "mid" region (cap_upper / cap_lower) and a "TE" region (cap_te_upper / cap_te_lower). |
| `--tip-inset FLOAT` | auto | Spanwise wall thickness of the cap.  Default = `1e-3 * chord`, where `chord` is the airfoil's largest in-plane extent at the tip cross-section. |
| `--n-cap-j INT` | `12` | Cap j-axis (airfoil-thickness) node count. |
| `--nk INT` | auto | Cap k-axis (spanwise) node count.  Default: matches the body block's `nk`. |
| `--cap-j-spacing {uniform,tanh,tanh2}` | `tanh` | Cap j-axis spacing law. |
| `--cap-j-beta FLOAT` | `4.0` | Clustering strength for tanh / tanh2 (cluster near the airfoil wall). |
| `--te-thickness-frac FLOAT` | `0.002` | Auto-round synthetic TE thickness as a fraction of chord (default 0.2%).  Only applied when the input cross-section is sharp.  The body wrap is left untouched. |
| `--te-thickness FLOAT` | — | Absolute auto-round TE thickness (mutually exclusive with `--te-thickness-frac`). |
| `--no-auto-round-te` | off | Reject sharp-TE inputs instead of auto-rounding. |
| `--span-axis {x,y,z}` | `z` | World axis the wing spans along. Used by the LE-cut chord-axis detector. |
| `-o` / `--output PATH` | `cap_<tip>` | Output stem.  Writes `<stem>_<label>.geo` per cap sub-block. |
| `--adjacency-json PATH` | — | Also write a `.adjacency.json` sidecar (for `shore cc-par` consumption). |
| `--journal PATH` | — | Also write the run as a replayable journal ([`shore run`](run)). |

## Output

The 5-block H-fill produces 5 `.geo` files per tip:

| File | Geometry |
|---|---|
| `<stem>_cap_le.geo` | LE-region wrap from upper-LE-cut to lower-LE-cut. |
| `<stem>_cap_upper.geo` | Suction-side chord band, between LE-cut and TE-cut, above the chord cut. |
| `<stem>_cap_lower.geo` | Pressure-side chord band, mirror of cap_upper. |
| `<stem>_cap_te_upper.geo` | Suction-side TE region, between TE-cut and the (synthetically blunt) TE-base. |
| `<stem>_cap_te_lower.geo` | Pressure-side TE region, mirror of cap_te_upper. |

Plus optionally `<adjacency-json-path>` describing the cap component's
internal SHARED seams (LE-cuts, TE-cuts, and chord-cut) for downstream
chimera assembly.

## Behaviour

### Cross-section extraction

Reads `BODY_GEO` (a `(ni, n_stations, nk, 3)` block from `shore mesh
--topology cgrid`) and extracts the spanwise-extreme cross-section
at the wall (k=0):

- `--tip jlo` → `body[:, 0, 0, :]`
- `--tip jhi` → `body[:, n_stations-1, 0, :]`

The next station inward (`body[:, 1, 0, :]` for jlo, `body[:, n_stations-2, 0, :]`
for jhi) is used to infer the **outward span direction** as the
negated, normalised vector between station centroids.

### LE-cut and TE-cut detection

The LE-cut is found by walking from the LE point along the suction
side until chord position equals `--le-frac`.  The lower-side LE-cut
is **snapped to be symmetric around the body's i-axis midpoint**
(`lower_idx = ni - 1 - upper_idx`), guaranteeing equal i-station
counts for `cap_upper` and `cap_lower` regardless of `ni_body` parity.
TE-cut detection uses the same chord-fraction snap with `--te-frac`.
This symmetry is required by the cap-fill's chord-cut SHARED face.

### Auto-round on sharp TE

Sharp-TE airfoils (i=0 and i=ni-1 of the cross-section coincide) are
**auto-rounded by default** — only inside the cap, the body wrap
remains sharp.  The cap's cross-section gets its TE node split into
two points spaced ``te_thickness`` apart so the H-fill has no singular
axis.  When auto-rounding fires, a one-line note is printed to
stderr.  See [C-grid topology — Auto-round on sharp TE](/topologies/c-grid#auto-round-on-sharp-te)
for the geometric details.

To override the synthetic thickness:

```bash
# Round to 0.5% chord (more visible artifact, more robust mesh)
shore cap-tip wing_body.geo --tip jlo --te-thickness-frac 0.005

# Or pass an absolute thickness in body units
shore cap-tip wing_body.geo --tip jlo --te-thickness 0.001

# Reject sharp inputs entirely
shore cap-tip wing_body.geo --tip jlo --no-auto-round-te
```

## Examples

```bash
# Build the C-grid first
shore mesh naca0012.stl \
    --topology cgrid --ni-body 120 --ni-wake 30 --n-stations 5 \
    --wake-length 15 --nk 20 --ds 5e-4 --growth 1.15 -o wing

# Cap the jlo tip
shore cap-tip wing_body.geo \
    --tip jlo --le-frac 0.05 --te-frac 0.95 --tip-inset 0.02 \
    --n-cap-j 12 --span-axis z \
    -o cap_jlo --adjacency-json cap_jlo.adjacency.json

# Cap the jhi tip
shore cap-tip wing_body.geo \
    --tip jhi --le-frac 0.05 --te-frac 0.95 --tip-inset 0.02 \
    --n-cap-j 12 --span-axis z \
    -o cap_jhi --adjacency-json cap_jhi.adjacency.json

# Bundle for ParaView (13 blocks: 3 C-grid + 5 cap_jlo + 5 cap_jhi)
shore export \
    wing_upper.geo wing_body.geo wing_lower.geo \
    cap_jlo_cap_le.geo cap_jlo_cap_upper.geo cap_jlo_cap_lower.geo \
    cap_jlo_cap_te_upper.geo cap_jlo_cap_te_lower.geo \
    cap_jhi_cap_le.geo cap_jhi_cap_upper.geo cap_jhi_cap_lower.geo \
    cap_jhi_cap_te_upper.geo cap_jhi_cap_te_lower.geo \
    -o wing_capped.vtm
```

The full pipeline lives in
[`examples/06-c-grid/visualize_naca0012_cap_tip.sh`](https://github.com/szaghi/shore/blob/main/examples/06-c-grid/visualize_naca0012_cap_tip.sh).

## Python API

For blunt TE or for finer control over the cap construction, drive
the topology directly from Python via [`CGridTipCap`](./api-topology#cgridtipcap)
and [`Mesh.from_array`](./api-mesh#meshfrom_array):

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
    tip_inset=0.02,
    n_cap_j=12,
    is_blunt=False,                 # set True + wake_base_cross_section for blunt
)
cap = Mesh.from_array(surface=cross, topology=topo, nk=body.shape[2])
cap.write_geo("cap_jlo")            # 3 .geo files
```

## See also

- [C-grid topology — Tip caps](/topologies/c-grid#tip-caps-closed-finite-wing)
- [`CGridTipCap`](./api-topology#cgridtipcap) — Python API.
- [`Mesh.from_array`](./api-mesh#meshfrom_array) — non-STL `Mesh` constructor.
- [`shore cc-par`](./cc-par) — chimera assembly that composes the cap
  against the C-grid body wrap.
- [`examples/06-c-grid/`](https://github.com/szaghi/shore/tree/main/examples/06-c-grid).
