---
name: shore
description: >
  How to operate the SHORE mesh generator (`shore` CLI, PyPI package
  `shore-mesh`) end-to-end: turn a triangulated STL surface (or a parametric
  primitive) into a structured body-fitted volume mesh via hyperbolic normal
  extrusion, then package it for the Xall / Xnavis Chimera (overset) CFD
  solver. Use this skill whenever the user wants to: generate a volume mesh
  from an STL or an analytic body; pick a topology (cubed-sphere, O-grid,
  C-grid, OH/CH hybrid); build a NACA section or extrude a wing; produce Xall
  solver inputs (`.geo`, `.grd`, `cc.par`, `proc.input`, marker boxes); split
  blocks or plan an MPI load balance; merge multiple Chimera components; or
  inspect / Jacobian-check / export a mesh to VTK. Also trigger on mentions of
  "shore mesh", "shore-mesh", "hyperbolic marching mesh", "cubed-sphere mesh",
  "Xall .grd / cc.par", or "Chimera overset preprocessing". This skill is the
  guide to *running* shore; for editing shore's own source, the repository
  CLAUDE.md is the architecture reference.
---

# Operating SHORE

SHORE (Structured Hyperbolic Overlapping Remeshing Engine) converts triangulated
STL surfaces into structured body-fitted volumes by hyperbolic normal extrusion,
and packages the result for the **Xall / Xnavis** Chimera (overset) CFD solver.
PyPI/install name: **`shore-mesh`**. Importable package and CLI: **`shore`**.

This skill is about **driving the tool** — the workflows, the command order, the
gotchas. It is not the architecture reference.

## Scope — operation vs. source

This skill and the repository `CLAUDE.md` are **complementary**, not competing:

- **`CLAUDE.md`** (loaded automatically in this repo) owns the *architecture*:
  module layout, the OOP block model, topology internals, conventions, and how
  to **edit** SHORE's source. Use it when changing the code.
- **This skill** owns *operation*: how to drive the `shore` CLI end-to-end to
  produce meshes and Xall files. Use it when running the tool, whether inside
  this repo or from any other directory.

When a question is "how do I make X with shore", this skill leads. When it is
"how is X implemented / where do I change it", `CLAUDE.md` leads.

## Flags are authoritative from `--help`, not from memory

The CLI evolves. This skill freezes the **workflows and command order**, which are
stable. For the **exact current flags** of any command, run:

```bash
shore <command> --help          # e.g. shore mesh --help
shore <group> <sub> --help      # e.g. shore body capsule --help
```

Never invent a flag. If a workflow below names a flag and `--help` disagrees,
`--help` wins — report the drift to the user.

## Environment

- SHORE requires a Python virtualenv. WSL2 system Python is externally managed.
- If a `.venv` exists in the project, activate it first:
  `source .venv/bin/activate`. Otherwise install into one:
  `python -m venv .venv && source .venv/bin/activate && pip install shore-mesh`.
  (In this repo, `make dev` creates `.venv` and installs editable with extras.)
- Optional extras: `shore-mesh[vis]` (PyVista interactive `view`),
  `shore-mesh[dev]` (test/lint). VTK *export* needs `pyvista`; the `view` window
  needs a display.
- Arrays are always `(ni, nj, nk, 3)`, last axis `(x, y, z)`. **`ni`/`nj`/`nk`
  are NODE counts** everywhere except the `.grd` header (cell counts = node − 1,
  converted automatically).

## The command surface

`shore --help` top-level commands, grouped by what they're for:

**Generate a mesh**
- `mesh` — full pipeline: STL → structured volume → `.geo`. The main entry point.
- `project` — k=0 surface layer only (lat-lon sphere ray-cast onto STL), no marching.

**Make geometry to feed it** (sub-apps)
- `body` — analytic STL bodies: `sphere`, `cylinder`, `torus`, `box`, `capsule`,
  `extrude` (stack a 2D section into a wing).
- `profile naca` — 2D NACA 4-/5-digit airfoil section `.dat`.
- `primitive` — parametric structured grids, no STL: `box`, `annulus`,
  `flat-caps`, `h-grid-rect`, `h-grid-tube`.

**Inspect / validate / export**
- `info` — metadata of a `.geo` block or `.npz` surface.
- `check` — layer-by-layer Jacobian validation of a `.geo`.
- `export` — `.geo`/`.npz` → VTK XML for ParaView.
- `view` — interactive PyVista window (needs `[vis]` + display).

**Xall solver bundle** (the overset preprocessing pipeline)
- `grd` — pack `.geo` blocks into a Fortran-binary `.grd`.
- `cc-par` — `.adjacency.json` sidecar → Xnavis `cc.par` connectivity.
- `proc-input` — MPI rank/group/body work distribution from a `.grd`.
- `*-merge` (`grd-merge`, `cc-par-merge`, `proc-input-merge`) — concatenate
  independent Chimera components into one Xall run.
- `boxes-from-stl` / `boxes-export` — Chimera hole-cutting marker boxes.

**Parallel decomposition**
- `split` — subdivide blocks along ijk, preserving the adjacency graph.
- `balance` — *plan* topology-safe splits to hit an MPI imbalance tolerance
  (emits a `shore split` plan; never mutates files).
- `cap-tip` — build a wing-tip cap from a C-grid body block.

**Config**
- `config` — inspect / template SHORE TOML configs.
- Precedence everywhere: explicit CLI flag → `--config` TOML → built-in default.

## Topology — the first real decision

`shore mesh --topology` selects one of four:

- **`cubed_sphere`** (default) — 6-block (4 equatorial + 2 polar caps). For
  closed, **star-shaped** bodies (sphere, fuselage, blob). Gnomonic flat-square
  cap seed → AR≈1 cap quads, no polar singularity pattern.
- **`ogrid`** — single block, j-periodic, polar singular axis rows. Simplest;
  one `.geo`. Also star-shaped bodies.
- **`cgrid`** — C-grid wrap for **sharp- or blunt-TE** bodies (airfoils,
  hydrofoils). 1/3/4-block. Needs `--ni-body`/`--ni-wake`/`--n-stations`/
  `--wake-length`. Optional `--cap-tips` for finite wings.
- **`ch`** — C-grid wrap + far-field H-channel hybrid (12-block sharp / 14
  blunt). Adds the frame flags `--ch-chord-min`/`--ch-chord-max`/
  `--ch-normal-min`/`--ch-normal-max` + `--nk-channel`.

Two more (`OHTopology`, `CGridTipCap`) are **Python-API-only** — not exposed on
`shore mesh`. Drive them from `shore.mesh`.

**Star-shaped requirement**: cubed-sphere and o-grid ray-cast from the body
centroid. If a ray from the centroid can't reach a surface point in one shot
(concave/occluded geometry), projection fails. C-grid/CH slice cross-sections
instead and tolerate non-star bodies.

## Spacing — three independent directions

- `--surface-i-spacing`, `--surface-j-spacing` — in-surface distributions.
- `--volume-k-spacing` — wall-normal extrusion. Laws: `geometric` (default),
  `tanh`, `tanh2`. **`tanh`/`tanh2` k-spacing require `--volume-k-thickness`.**
- First-cell control via `--ds` (first wall-normal step) and `--growth`.
- All blocks in one mesh share ONE k-step schedule (k-coupled seams must stay
  coincident layer-by-layer); block 0 wins, `UserWarning` if others disagree.

## Theta-cap (cubed-sphere / o-grid)

The polar-cap half-angle has mutually-exclusive aliases: `--theta-cap`
(radians) **xor** `--theta-cap-deg` (degrees). Pass exactly one.

## Workflows

### A. STL → body-fitted volume → `.geo`

```bash
shore mesh body.stl -o run/wall \
  --topology cubed_sphere \
  --ni 64 --nj 64 --nk 48 \
  --ds 1e-3 --growth 1.08 \
  --volume-k-spacing geometric
# check flags: shore mesh --help
```

Then validate before trusting it:

```bash
shore check run/wall.geo          # Jacobian per layer; negative = tangled
shore info  run/wall.geo
shore export run/wall.geo -o run/wall.vts   # ParaView
```

Marching aborts immediately on a negative Jacobian. If it fails: reduce `--ds`,
lower `--growth`, increase smoothing, or reduce `--nk`. A non-star body is the
other common cause — switch to `cgrid` or fix the geometry.

### B. No CAD? Build the body first

```bash
shore body sphere   -o sphere.stl --radius 1
shore body capsule  -o capsule.stl --help     # see params
# Wing: section then extrude. The section is passed via --section (a flag),
# NOT positionally.
shore profile naca 0012 -o n0012.dat
shore body extrude --section n0012.dat -o wing.stl --span 4 --n-stations 20
shore mesh wing.stl -o run/wing --topology cgrid \
  --ni-body 128 --ni-wake 32 --n-stations 20 --wake-length 5
```

### C. Background grid with no STL (parametric primitive)

```bash
shore primitive flat-caps   -o bg/cyl --help    # 5-block Chimera-free background
shore primitive box         -o bg/box --help
shore primitive h-grid-tube -o bg/tube --help   # internal flow, j-periodic
```

### D. Full Xall Chimera bundle (the payoff)

Order matters. From per-block `.geo` + the `.adjacency.json` sidecar SHORE
writes alongside a multi-block mesh. **`cc-par` and `proc-input` take the input
and output as POSITIONAL args** (`ADJACENCY OUTPUT`, `GRD OUTPUT`) — not `-o`.

```bash
# 1. pack geometry
shore grd run/wall*.geo -o run/wall.grd
# 2. hole-cutting marker boxes from the body STL (feeds cc-par next).
#    --block is a SHORE block LABEL (e.g. sub0), not a .geo path;
#    --max-outward-gap is required.
shore boxes-from-stl body.stl --block sub0 --max-outward-gap 0.10 \
  -o run/boxes.json
# 3. connectivity  (.adjacency.json -> cc.par); the .grd name is --grd (NOT
#    --grd-filename), boxes via --boxes.
shore cc-par run/wall.adjacency.json run/cc.par \
  --grd wall.grd --boxes run/boxes.json
# 4. MPI work distribution (weights = .grd cell counts); summary on stderr
shore proc-input run/wall.grd run/proc.input --np 16
```

Multi-component overset (near-body wall + background, each its own mesh):

```bash
shore grd-merge        wall.grd bg.grd       -o assembly.grd
shore cc-par-merge     wall.cc.par bg.cc.par -o assembly.cc.par
shore proc-input-merge wall.proc bg.proc     -o assembly.proc.input
```

`merge` shifts block/patch indices across components; it does **not** re-balance.

### E. MPI load balancing

```bash
shore balance run/wall.grd --np 32 --tolerance 0.1 \
  --adjacency run/wall.adjacency.json -o split.toml   # PLAN only
shore split run/wall*.geo -c split.toml -o run/balanced/  # apply the plan
```

`balance` prefers k-axis cuts (k-faces are WALL/FREE, never SHARED), falling
back to paired j/i co-splits. It never mutates files — it emits a `shore split`
plan you then run.

### F. Config-driven runs

```bash
shore config template > mesh.toml          # emit the default mesh config
shore mesh body.stl -c mesh.toml -o run/wall   # CLI flags still override TOML
```

`config` takes an `ACTION` (`template` | `validate`) and, for `validate`, a
path: `shore config validate mesh.toml`. `template` always emits the `shore
mesh` config to stdout.

## Common failure modes

- **`MarchingError` / negative Jacobian**: step too aggressive or non-star body.
  Smaller `--ds`, smaller `--growth`, fewer layers, more smoothing; or switch
  topology (cgrid for foils/non-star).
- **`MeshInputError` on a non-watertight / multi-mesh STL**: SHORE wants a single
  watertight mesh. It warns on non-watertight; repair upstream or accept the warn.
- **`tanh`/`tanh2` k-spacing without `--volume-k-thickness`**: required pair.
- **Both `--theta-cap` and `--theta-cap-deg`**: mutually exclusive; pass one.
- **`view` does nothing / errors**: needs `shore-mesh[vis]` and a display; in
  headless/CI use `export` to `.vts`/`.vtm` and open in ParaView instead.
- **Wrong importable name**: install is `pip install shore-mesh`, import/CLI is
  `shore` — not `import shore_mesh`.

## When you need depth

The canonical end-to-end scripts live under `examples/` (numbered
`01-body-fitted` … `09-ch-hybrid`), each a runnable, documented workflow. The
prose guide is at `https://szaghi.github.io/shore/`. Point users at
`examples/<n>-.../` for a copy-pasteable reference rather than reconstructing a
workflow from scratch.
