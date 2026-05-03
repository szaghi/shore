# SHORE examples

Standalone scripts and a Jupyter notebook that demonstrate every
feature SHORE ships. Each script is invokable from the project
root after activating `.venv`:

```bash
python examples/<subdir>/<script>.py [--no-view] [-o my_run]
```

Most scripts accept `--no-view` to skip the interactive PyVista
window (CI-friendly), `--screenshot path.png` to save an
off-screen render, and `-o output_dir` / `-c config.toml` where
applicable. See each script's `--help` for the full option list.

The full prose tour is at the
[Running examples and tests](https://szaghi.github.io/shore/guide/running)
guide.

## Layout

The 5 subdirectories are ordered from "first SHORE script to read"
to "complete end-to-end Xall pipeline". Read top-to-bottom on a
fresh checkout.

### [`01-body-fitted/`](01-body-fitted/) — body-fitted meshing of an STL

| File | What it does |
|---|---|
| `visualize_cubed_sphere.py` | Build a 6-block cubed-sphere or single-block O-grid mesh from an STL; accepts every `shore mesh` flag on the CLI. Without `--generate` it acts as a generic `.geo` / `.npz` viewer. |
| `config_fuselage.toml` | Annotated TOML config for `shore mesh -c`. Demonstrates every section (tanh wall-normal spacing, equator clustering at the cap seams, blended STL normals). |
| `tour.ipynb` | Jupyter notebook covering generation, Jacobian analysis, spacing-law comparison, and the C1 cap-equator seam-matching demo. |

### [`02-primitives/`](02-primitives/) — parametric backgrounds

No STL needed — primitives build their nodes from analytic
geometry.

| File | Primitive | Default output |
|---|---|---|
| `visualize_box.py` | Cartesian box | `box.geo`, `box.vtm` |
| `visualize_annulus.py` | Annular cylinder | `annulus.geo`, `annulus.vtm` |
| `visualize_flat_caps.py` | 5-block flat-caps cylinder | `flat_caps_<label>.geo` × 5, `flat_caps.vtm` |
| `config_primitive_box.toml` | TOML config for `shore primitive box -c` | — |
| `config_primitive_annulus.toml` | TOML config for `shore primitive annulus -c` | — |
| `config_primitive_flat_caps.toml` | TOML config for `shore primitive flat-caps -c` | — |

### [`03-per-edge-spacing/`](03-per-edge-spacing/) — per-edge spacing data model

Three scripts that demonstrate the per-edge spacing override
machinery — `Spacing1D` instances replaced on individual edges,
with TFI blends or identity-check overrides depending on the
topology. See the
[per-edge spacing guide](https://szaghi.github.io/shore/guide/per-edge-spacing).

| File | What it does |
|---|---|
| `per_edge_ij_box.py` | Per-edge i / j on `BoxMesh` via TFI blend between the 4 parallel edges of an axis. Compare the pinned and uniform corners' first-cell spacing. |
| `per_edge_k_spacing_cubed_sphere.py` | Per-edge k on a cubed-sphere mesh. Replace all 4 k-edges of every block with a shared `Spacing1D` and `Mesh.march()` picks it up via identity check. |
| `per_edge_k_spacing_ogrid.py` | Same idea on a single-block O-grid. |

### [`04-split/`](04-split/) — parallel decomposition

`shore split` subdivides existing structured blocks along ijk
indices, preserving the adjacency graph (internal SHARED seams,
external partner rewiring). Used for parallel decomposition of
large blocks.

| File | What it does |
|---|---|
| `split_cubed_sphere.py` | Build a cubed-sphere mesh, split each of the 6 blocks along k, and write the per-chunk `.geo` files plus `.adjacency.json` sidecar. |
| `config_split_cubed_sphere.toml` | TOML config for `shore split -c`. |

### [`05-chimera-assembly/`](05-chimera-assembly/) — full Xall pipeline

The end-to-end overset assembly: a body-fitted near-body wall + a
flat-caps cylinder background, packed into the Xall solver bundle
(`.grd` + `cc.par` + marker boxes) and merged into a single
multi-component job.

| File | What it does |
|---|---|
| `assembly_sphere_in_cylinder.py` | Build the two meshes (cubed-sphere wall + flat-caps cylinder background) sharing one world domain and bundle them into a single `assembly.vtm` for ParaView. The reference for the [assembly walkthrough](https://szaghi.github.io/shore/guide/assembly). |
| `chimera_assembly.py` | Continuation of the above: emit per-component `.grd`, `.adjacency.json`, marker boxes from the body STL, and `cc.par`; merge into `assembly.grd` + `assembly.cc.par`. The reference for the [cc.par + boxes pipeline](https://szaghi.github.io/shore/guide/cc-par-pipeline). |
| `config_assembly_near_body.toml` | TOML config for the near-body cubed-sphere mesh. |
| `config_assembly_background.toml` | TOML config for the flat-caps cylinder background. |

### [`08-oh-hybrid/`](08-oh-hybrid/) — OH topology (cubed-sphere wrap + H-channel)

Style-2 conforming hybrid: 6 cubed-sphere wrap blocks (marched
outward from the STL) plus 6 channel-frustum blocks (ruled-surface
from the wrap envelope to a user-supplied bounding box).  Wrap and
channel are SHARED-by-view at the seam — one numpy buffer per
pair, bit-exact with no copy.

| File | What it does |
|---|---|
| `visualize_oh_hybrid.py` | Build a 12-block OH hybrid mesh on a unit icosphere; CLI knobs for box extents, channel depth, and channel k-spacing law.  Writes 12 `oh_<label>.geo` blocks plus an `oh.vtm` MultiBlock. |

`OHTopology` is Python-only at the moment; `shore mesh` does not
yet accept ``--topology oh-hybrid``.  Drive it from Python via
:class:`shore.mesh.OHTopology`.

## Adding a new example

Pick the subdirectory that matches the feature being demonstrated.
If the example doesn't fit in any existing bucket, propose a new
numbered subdir in the PR (don't drop scripts at the top level).
Every example script should:

- Accept `--no-view` for headless / CI use.
- Accept `-o <stem>` for the output directory.
- Carry a docstring at the top with the invocation patterns and a
  one-paragraph rationale (what the script demonstrates and why
  it matters).
