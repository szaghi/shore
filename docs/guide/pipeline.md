# Pipeline overview

SHORE turns a triangulated STL surface into a structured body-fitted
volume mesh through a fixed pipeline. This page is the high-level map;
the [Algorithm](/algorithm/) section drills into each stage, and the
[Topologies](/topologies/) section catalogues how the stages are
assembled into block layouts.

![Pipeline](/figures/pipeline.png)

## Stages

| Stage | Output | Algorithm reference |
|---|---|---|
| 1. Surface projection | `(N, 3)` k=0 surface nodes | [projection](/algorithm/projection) |
| 2. k=0 mesh construction | 6-block surface (cubed-sphere topology) | [k0-construction](/algorithm/k0-construction) |
| 3. Hyperbolic marching | 6-block volume `(ni, nj, nk, 3)` per block | [hyperbolic-marching](/algorithm/hyperbolic-marching) |
| 4. Seam-aware refinement | Continuous spacing across multi-block seams at every layer | [seam-aware-normals](/algorithm/seam-aware-normals) |

Stages 1 and 2 produce the **k=0 layer** — the body-conforming
surface mesh from which the volume is extruded. Stage 3 is the
hyperbolic extrusion. Stage 4 is integrated into stage 3: it changes
how each block's normals are computed at multi-block seams so that
spacing continuity is preserved layer by layer.

## Topology

The pipeline above is run **per topology**. The topology determines
the block layout and how stages 2 and 3 stitch the blocks together.
SHORE's production default is the
[cubed-sphere topology](/topologies/cubed-sphere) (6 blocks). Other
CLI-selectable topologies are the single-block
[O-grid](/topologies/o-grid) (`--topology ogrid`), the airfoil
[C-grid](/topologies/c-grid) (`--topology cgrid`, 3-/4-block sharp /
blunt TE), and the airfoil [CH hybrid](/topologies/ch-hybrid)
(`--topology ch`, C-grid wrap + far-field channel, 12-/14-block).
The [OH hybrid](/topologies/oh-hybrid) (`--topology oh`) extends the
cubed-sphere wrap with a 6-face H-channel. The full catalogue lives
in the [Topologies](/topologies/) section.

## Inputs

| | |
|---|---|
| **STL** | A single watertight `trimesh.Trimesh`. The body must be **star-shaped from the anchor** (default centroid). See [Surface projection](/algorithm/projection#failure-modes) for the failure modes. |
| **Resolution** | `ni` (equator meridional), `nj` (longitude), `nk` (wall-normal layers including k=0). |
| **Cap angle** | `theta_cap` (default $0.05\,\pi$ = 9°) — the colatitude of the cap-equator junction parallel. |
| **Spacing laws** | `i_spacing`, `j_spacing`, `volume_k_spacing` — see [Spacing laws](/algorithm/spacing-laws). |

## Outputs

The CLI writes 6 `.geo` files (one per block) plus optionally a
ParaView `.vtm` MultiBlock if `shore export` is used:

```
<stem>_sub0.geo .. <stem>_sub3.geo   # 4 equatorial sub-blocks
<stem>_cap_north.geo                 # north polar cap
<stem>_cap_south.geo                 # south polar cap
```

The `.geo` format is Xall compatible: ASCII header `ni nj nk`,
then x, y, z arrays in Fortran column-major order.

## Visualization

Two independent paths after the volume is generated:

### ParaView (VTK export)

```bash
# Per-block .vts
shore export sphere_sub0.geo

# Combined MultiBlock for the 6-block cubed-sphere topology
python -c "
from shore.io.geo import read_geo
from shore.io.vtk import write_caps_vtm
labels = ['sub0', 'sub1', 'sub2', 'sub3', 'cap_north', 'cap_south']
blocks = [read_geo(f'sphere_{lbl}.geo') for lbl in labels]
write_caps_vtm('sphere.vtm', blocks, with_jacobian=True)
"
paraview sphere.vtm
```

`shore export` attaches a `Jacobian` cell-data array; the
`write_caps_vtm` helper bundles all 6 blocks under one MultiBlock
with the same metadata.

### PyVista (interactive / Jupyter)

```bash
shore view sphere_sub0.geo   # desktop window
```

From Python or Jupyter:

```python
import pyvista as pv
from shore.io.geo import read_geo
from shore.io.vtk import open_to_pyvista

grid = read_geo("sphere_sub0.geo")
sg = open_to_pyvista(grid, with_jacobian=True)
pl = pv.Plotter()
pl.add_mesh(sg, scalars="Jacobian", cmap="RdYlGn", show_edges=True, opacity=0.7)
pl.show()
```

## Limitations

### Star-shaped constraint

The surface projection requires that every ray from the body anchor
hits the STL exactly once. Bodies with re-entrant cavities, holes,
or multiple disconnected components must be repaired or decomposed
first. See [Surface projection](/algorithm/projection#failure-modes).

### No smooth far field

The outer boundary of the volume (k=nk-1) currently inherits the
shape of the cap-equator seam topology — it is not morphed to a box
or matched to the Xall background grid. Far-field morphing is
the next planned feature.

### Cap-equator size compatibility

The C1 cap-equator seam pin (per-j first-step matching) requires
`(ni, nj, theta_cap)` to land in a compatibility window. Out of
window, the equator falls back to its unpinned `i_spacing` law and
emits one `UserWarning` per build. See
[k=0 mesh construction](/algorithm/k0-construction#compatibility-window-and-fallback)
for the rule of thumb.

## See also

- [Algorithm overview](/algorithm/)
- [Topology catalogue](/topologies/)
- [Quick start](/guide/quickstart) — worked CLI examples
- [Chimera context](/guide/chimera) — overset mesh generation
