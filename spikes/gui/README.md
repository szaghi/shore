# GUI spike (ADR 0006)

Throwaway code that checks the PySide6 + pygfx/wgpu stack against the criteria
in [ADR 0006](../../docs/adr/0006-gui-stack.md) before the ADR is accepted.
It is not part of the `shore` package and is not imported by it.

| File | Checks |
|---|---|
| `synthetic.py` | 20-block spherical-shell assembly, N nodes, geometric wall spacing `ds_wall`, world offset |
| `viewport.py` | `Viewport` protocol + pygfx backend (the only pygfx import); three origin modes |
| `bench_headless.py` | items 2–3 offscreen: upload/first frame, orbit frame time, float32 pick trap → `results/*.json` |
| `latency.py` | item 4: UI heartbeat lateness during numpy, GIL-holding, and subprocess workloads |
| `app.py` | items 1–4 on screen: docked layout, OS dark/light theme, orbit bench, double-click pick, fake march, animated zoom to a wall cell with rebase log |
| `interaction.py` | simulated user zoom (camera + controller target only) for testing the rebase trigger |

## Setup (each OS)

```bash
python -m venv .venv                 # Windows: py -m venv .venv
.venv/bin/pip install -r requirements.txt   # Windows: .venv\Scripts\pip
```

## Run

```bash
python bench_headless.py --nodes 10000000          # frames, trap, rebase trigger, drift, sweep -> results/*.json
python bench_headless.py --skip-frames --skip-sweep # trap only (quick)
python latency.py --workload numpy --seconds 30
python latency.py --workload gil --seconds 10
python latency.py --workload gil-process --seconds 10
python app.py --nodes 10000000 --no-vsync          # Load, Orbit bench, Fake march; double-click to pick
QT_QPA_PLATFORM=offscreen python app.py --smoke    # scripted end-to-end check
```

Check the adapter line first. A `(CPU)` adapter (llvmpipe/lavapipe) means
wgpu did not find the GPU, and the frame times are meaningless for the budget.

## Results

### adam: WSL2 (Ubuntu 24.04, WSLg), 2x RTX 4070, 2026-10-25

| Item | Result |
|---|---|
| wgpu adapter | **llvmpipe (CPU) only.** Ubuntu 24.04 Mesa 25.2 ships no `dzn` Vulkan ICD, so the RTX 4070 is invisible to wgpu |
| 1. Window, docks, theme | Works (offscreen QPA smoke run); on-screen look still to be checked by eye under WSLg |
| 2. 10⁷ nodes (749k surface nodes), upload + first frame | 0.89 s (budget 2 s): **PASS** |
| 2. Orbit frame time, 10⁷ nodes, offscreen incl. readback | mean 158 ms, p95 171 ms (budget 16 ms): **FAIL**, CPU rasteriser, expected |
| 2. Orbit frame time, 10⁵ nodes, Qt widget | mean 49 ms, p95 55 ms: **FAIL**, CPU rasteriser |
| 3. Pick latency | 0.8 ms (budget 50 ms): **PASS** |
| 3. Float32 trap, `ds_wall` 1e-6 at offset 1e3 | `absolute` FAIL, `block` FAIL, `floating` **PASS** |
| 3. Extent x `ds_wall` sweep, AR 1e4, `ds/extent` 1e-7 to 1e-10 | `floating` fails below ~1e-8; `tiled` partial and costly; **`rebase` 18/18 everywhere** (see below) |
| 3. Rebase cost, 10⁷ nodes (749k surface vertices) | 3.8 ms CPU, no measurable upload penalty on the next draw |
| 3. Rebase trigger, user-style zoom onto 3 wall cells (`ds/extent` 1e-9) | factor 100: 3/3 picks, 1 rebase per zoom; trigger off: **0/3** |
| 3. Origin drift tolerance (trigger off), 20 px wall cell | exact to 1e4 view widths (predicted 0.86 px); 2/3 at 1e5 (8.6 px); 0/3 at 1e6 (86 px) |
| 4. numpy march in worker, UI lateness | max 1.8 ms (budget 100 ms): **PASS** |
| 4. one long GIL-holding C call in worker | max 705 ms: **FAIL**. The same call in a subprocess: max 1.9 ms, **PASS** |

### adam: Windows 11 native Python 3.14.3, RTX 4070 (wgpu via Vulkan), 2026-10-29

Run from WSL through `py.exe` on a copy of this directory (`results/adam-NVIDIA-GeForce-RTX-4070-DiscreteGPU-via-Vulkan.json`).

| Item | Result |
|---|---|
| wgpu adapter | **NVIDIA GeForce RTX 4070 (DiscreteGPU) via Vulkan** (high-performance and low-power alike) |
| 1. Window, docks, theme | offscreen QPA smoke run passes; **on-screen look still to be checked by eye** |
| 2. 10⁷ nodes (749k surface nodes), upload + first frame | 1.15 s (budget 2 s): **PASS** |
| 2. Orbit frame time, 10⁷ nodes, offscreen incl. readback | mean 6.1 ms, p95 8.9 ms, max 11.3 ms (budget 16 ms): **PASS** |
| 2. Orbit frame time, app smoke (98k nodes, Qt widget offscreen) | mean 5.2 ms, p95 7.1 ms: **PASS** |
| 3. Float32 trap | `absolute` FAIL, `block` FAIL, `floating` PASS (as on WSL) |
| 3. Sweep `ds/extent` 1e-7 to 1e-10, extents 2.85 and 285 | **`rebase` 18/18 everywhere**, draw 8-14 ms, rebase 1.0-1.2 ms; `floating` 2-12/18 below 1e-8; `tiled` up to 7012 meshes, 410 ms draws |
| 3. Rebase trigger (factor 100) | 3/3 picks, 3 rebases; trigger off 0/3 |
| 3. App zoom to a wall cell (2.85e-7) | 1 rebase, pick error 0: **PASS** |
| 4. numpy march in worker, UI lateness | max 14.9 ms (budget 100 ms): **PASS** |
| 4. GIL-holding call in a subprocess | **PASS** |

### Still to run

- [x] Windows 11 native Python, RTX 4070 (Vulkan), 2026-10-29: above
- [x] Windows: the on-screen window (docks, theme, HiDPI) checked by eye: `python app.py --nodes 10000000` (user, 2026-10-29: looks fine)
- [ ] macOS, Apple Silicon (Metal)
- [ ] Linux native, discrete GPU (Vulkan)
- [ ] WSL2 with a `dzn` ICD (Mesa ≥ 23 built with `-Dvulkan-drivers=microsoft-experimental`), to see if WSL2 can be a GUI dev box at all
- [ ] Rebase trigger checked by eye on a native GPU: pan and zoom on the default `app.py` assembly, watching the console's `rebase(...)` lines

### Extent x `ds_wall` sweep (adam, offset 1e3, wall-cell aspect ratio 1e4)

18 picks per row: 3 blocks x the 2 wall cells adjacent to a tangential
cluster centre x k = 0..2, on `i_lo` faces of a non-axis-aligned (`n_phi=5`)
assembly. Polar angle is clustered so the picked cells have a realistic
aspect ratio (tangential / normal = 1e4) at every `ds_wall`.
`res32 = eps32 * max|p - block centroid|`. `draw` is offscreen on llvmpipe
(CPU), so only its relative growth with mesh count means anything.

| extent | `ds/extent` | `ds/res32` | `floating` | `tiled` (meshes, setup, draw) | `rebase` (rebase cost) |
|---:|---:|---:|---:|---|---|
| 2.85 | 1e-7 | 4 | 18/18 | 18/18 (200, 0.1 s, 45 ms) | 18/18 (0.5 ms) |
| 2.85 | 1e-8 | 0.4 | 18/18 | 18/18 (1050, 0.6 s, 75 ms) | 18/18 (0.5 ms) |
| 2.85 | 1e-9 | 0.04 | 10/18 | 14/18 (4852, 2.9 s, 203 ms) | 18/18 (0.5 ms) |
| 2.85 | 1e-10 | 0.004 | 1/18 | 12/18 (7012, 4.1 s, 282 ms) | 18/18 (0.5 ms) |
| 285 | 1e-7 | 4 | 18/18 | 18/18 (200, 0.1 s, 46 ms) | 18/18 (0.5 ms) |
| 285 | 1e-8 | 0.4 | 18/18 | 18/18 (1050, 1.4 s, 81 ms) | 18/18 (0.6 ms) |
| 285 | 1e-9 | 0.04 | 12/18 | 17/18 (4852, 2.9 s, 210 ms) | 18/18 (0.6 ms) |
| 285 | 1e-10 | 0.004 | 2/18 | 17/18 (7012, 4.5 s, 271 ms) | 18/18 (0.5 ms) |

`floating` (20 meshes) holds while `ds/res32 >~ 0.4`, i.e.
`ds/extent >~ 1e-8`. With realistic aspect ratios the threshold is lower than
the 1.4 of the first sweep, whose coarse assembly gave aspect ratios up to
~1e9. `tiled` needs thousands of meshes at `ds/extent <= 1e-9` and still
misses picks (tiles near the cluster bottom out at a few cells, where the
cell's own extent dominates). `rebase` is exact at every ratio with no extra
meshes.

Why `rebase` wins: after `float32(p - origin)`, a vertex at distance `d`
from the camera target is in error by ~`eps32 * d`. Anything visible lies
within about one view width `V` of the target, and a pixel is ~`V / 1000`.
So visible geometry is exact to ~1e-4 pixel at any scale, and a rebase is
only needed when the target drifts more than ~`1 / (1000 eps32)` ~ 8e3 view
widths from the origin, after a large pan at high zoom.

### Rebase trigger (adam)

`Viewport.maybe_rebase()` runs at the start of every draw. It reads the orbit
controller's live target, so user pan and zoom are covered, and moves the
origin to the target when the drift exceeds `rebase_factor` view widths at
the target (default 100). The factor comes from the bound below: the drift
error is `eps32 * drift`, against a 0.1 px goal at ~1000 px, which gives
`drift <~ 840 V`. The factor of 100 keeps a ~8x margin.

| Check | Result |
|---|---|
| `drift_tolerance_check`: origin D view widths from a 20 px wall cell, trigger off | D = 1e2 / 1e3 / 1e4: 3/3 (predicted error 0.009 / 0.086 / 0.86 px); 1e5: 2/3 (8.6 px); 1e6: 0/3 (86 px). Failures start where the predicted error reaches the ~10 px pick tolerance, so the `eps32 * D` model holds |
| `rebase_policy_check`: 60-frame zoom from 3 radii onto the k=0 wall cell of 3 blocks, `ds/extent = 1e-9` | factor 100: 3/3, one rebase per zoom, fired at ~100 view widths; factor inf (off): 0/3 |
| `app.py --smoke` | zoom to wall: 1 rebase at 102 view widths, 0.6 ms, pick exact |

A zoom cannot probe the factor bound. It keeps zooming after the trigger
fires, so it always ends with the origin at the target (factors up to 1e6
all pass a zoom). The bound only bites after a long pan at high zoom, which
is what the drift check measures directly. In `floating` mode the trigger
moves the origin without re-uploading, which is why its logged cost is 0 ms
and why one origin per block still limits its precision. The single-cell app
demo does not discriminate the modes (that cell happens to pass in
`floating`); the three-cell bench check does.

## Findings that change ADR 0006

1. **A per-block centroid offset does not avoid the float32 trap.** pygfx
   applies `world_transform` and then `cam_transform` on the GPU in float32,
   so a 1e3 translation cancels in float32 anyway. A **floating origin**
   (object translations relative to a camera-following origin, computed in
   float64 on the CPU) is required.
2. **Long OCP calls must run in a subprocess, not a thread.** pybind11 does
   not release the GIL by default, and a single long call such as
   `BRepMesh_IncrementalMesh` on a large STEP freezes the UI for its whole
   duration (simulated: 705 ms).
3. **WSL2 is not a GUI development box** in its default Ubuntu 24.04 setup.
   GUI work happens on native Windows Python against the same repo.
4. `OrbitController.target`'s setter calls `camera.look_at` in float32, and
   `look_at` is degenerate when the view direction is parallel to `up`. Both
   bit the spike. The `Viewport` must own camera state and treat the pygfx
   controller as an input device only.
5. **Camera-relative rebase replaces per-tile local origins.** One origin
   per block holds only while `ds_wall / extent >~ 1e-8`. For y+ ~ 1 the
   first cell is `ds/L ~ 1 / (Re * sqrt(Cf/2))`: ~4e-8 at Re = 1e9 and ~5e-9
   at Re = 1e10. (An earlier version of this README said 1e-9 to 1e-10, which
   was wrong.) Large full-scale cases are therefore at or past the threshold.
   - Per-tile origins (prototyped as `tiled`) need thousands of meshes and
     multi-second setup, and still miss picks.
   - **Rebase** (`float32(p - origin)` re-uploaded from float64 when the
     origin moves) is exact at every ratio tested, keeps one mesh per block,
     and costs 3.8 ms for 10⁷ nodes.

   Relative-to-eye rendering with emulated double precision is no longer
   needed.
