#!/usr/bin/env bash
# chimera_assembly.sh — bash equivalent of chimera_assembly.py.
#
# Builds the same sphere-in-cylinder Chimera assembly using only the
# `shore` CLI (no Python imports).  Demonstrates the pure-CLI path
# from STL through assembled Xall job:
#
#   1. shore body sphere     -> body STL
#   2. shore mesh ... --topology cubed_sphere --adjacency-json ...
#                            -> 6-block near-body wall + sidecar
#   3. shore primitive flat-caps ... --adjacency-json ...
#                            -> 5-block background + sidecar
#   4. shore boxes-from-stl  -> marker boxes JSON
#   5. shore export          -> ParaView .vtm files
#   6. shore boxes-export    -> standalone boxes .vtm
#   7. shore grd             -> per-component binary .grd
#   8. shore cc-par          -> per-component cc.par (with --meta + --boxes)
#   9. shore proc-input      -> per-component proc.input
#  10. shore grd-merge / cc-par-merge / proc-input-merge
#                            -> assembly.grd + assembly.cc.par + assembly.proc.input
#
# Usage:
#   ./examples/05-chimera-assembly/chimera_assembly.sh             # default examples_out/
#   ./examples/05-chimera-assembly/chimera_assembly.sh -o out_dir
#   ./examples/05-chimera-assembly/chimera_assembly.sh --np 16
#
# This script is the canonical reference for "what does an end-to-end
# chimera build look like in shell?".  See chimera_assembly.py for the
# Python API equivalent (with optional --balance pass + interactive
# viewer).

set -euo pipefail

# ── Mesh parameters (match chimera_assembly.py defaults exactly) ──────────────
BODY_RADIUS=0.5
ICOSPHERE_SUBDIVISIONS=4

# Near-body cubed-sphere
NB_NI=30
NB_NJ=40
NB_NK=13
NB_THETA_CAP_DEG=30.0
NB_DS=0.01
NB_GROWTH=1.1

# Background flat-caps cylinder
BG_R_OUT=10.0
BG_Z0=-10.0
BG_Z1=50.0
BG_INNER_FRAC=0.4
BG_NI=16
BG_NC=24
BG_NK=60

# ── Argument parsing ──────────────────────────────────────────────────────────
OUT_DIR="examples_out"
NP=4
SHORE="${SHORE:-shore}"

usage() {
  sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o | --output)   OUT_DIR="$2"; shift 2 ;;
    --np)            NP="$2"; shift 2 ;;
    -h | --help)     usage 0 ;;
    *) echo "Unknown option: $1" >&2; usage 1 ;;
  esac
done

cyan()  { printf '\033[1;36m%s\033[0m\n' "$*"; }
red()   { printf '\033[1;31m%s\033[0m\n' "$*" >&2; }
info()  { cyan "==> $*"; }
err()   { red "ERROR: $*"; exit 1; }

command -v "$SHORE" >/dev/null 2>&1 \
  || err "'$SHORE' not on PATH.  Install with: pipx install shore-mesh"

mkdir -p "$OUT_DIR"
info "output directory: $OUT_DIR"

# ── 1. Body STL (sphere) ──────────────────────────────────────────────────────
info "1. body STL"
"$SHORE" body sphere \
  --radius "$BODY_RADIUS" \
  --subdivisions "$ICOSPHERE_SUBDIVISIONS" \
  -o "$OUT_DIR/sphere.stl"

# ── 2. Near-body cubed-sphere mesh + adjacency sidecar ────────────────────────
info "2. near-body cubed-sphere (${NB_NI}x${NB_NJ}, nk=${NB_NK})"
"$SHORE" mesh "$OUT_DIR/sphere.stl" \
  --topology cubed_sphere \
  --ni "$NB_NI" --nj "$NB_NJ" --nk "$NB_NK" \
  --ds "$NB_DS" --growth "$NB_GROWTH" \
  --theta-cap-deg "$NB_THETA_CAP_DEG" \
  --smooth 0.2 --smooth-iters 2 \
  -o "$OUT_DIR/wall" \
  --adjacency-json "$OUT_DIR/wall.adjacency.json"

# ── 3. Background flat-caps cylinder + adjacency sidecar ──────────────────────
info "3. background flat-caps cylinder (ni=${BG_NI}, nc=${BG_NC}, nk=${BG_NK})"
"$SHORE" primitive flat-caps \
  --r-out "$BG_R_OUT" \
  --z0 "$BG_Z0" --z1 "$BG_Z1" \
  --inner-frac "$BG_INNER_FRAC" \
  --ni "$BG_NI" --nc "$BG_NC" --nk "$BG_NK" \
  -o "$OUT_DIR/background" \
  --adjacency-json "$OUT_DIR/background.adjacency.json"

# ── 4. Marker boxes from the body STL ─────────────────────────────────────────
info "4. marker boxes from STL"
"$SHORE" boxes-from-stl "$OUT_DIR/sphere.stl" \
  --block sub0 \
  --max-outward-gap 0.10 \
  --method voxel-fill \
  -o "$OUT_DIR/sphere_boxes.json"

# ── 5. ParaView .vtm exports for the per-component meshes ─────────────────────
info "5. ParaView .vtm exports"
"$SHORE" export \
  "$OUT_DIR/wall_sub0.geo" "$OUT_DIR/wall_sub1.geo" \
  "$OUT_DIR/wall_sub2.geo" "$OUT_DIR/wall_sub3.geo" \
  "$OUT_DIR/wall_cap_north.geo" "$OUT_DIR/wall_cap_south.geo" \
  -o "$OUT_DIR/wall.vtm" || info "  (pyvista not installed; skipped wall.vtm)"

"$SHORE" export \
  "$OUT_DIR/background_center.geo" \
  "$OUT_DIR/background_sub_e.geo" "$OUT_DIR/background_sub_n.geo" \
  "$OUT_DIR/background_sub_w.geo" "$OUT_DIR/background_sub_s.geo" \
  -o "$OUT_DIR/background.vtm" || info "  (pyvista not installed; skipped background.vtm)"

# ── 6. Standalone marker-boxes .vtm ───────────────────────────────────────────
info "6. boxes-export"
"$SHORE" boxes-export "$OUT_DIR/sphere_boxes.json" \
  -o "$OUT_DIR/sphere_boxes.vtm" || info "  (pyvista not installed; skipped sphere_boxes.vtm)"

# ── 7. Per-component .grd ─────────────────────────────────────────────────────
info "7. binary .grd files"
"$SHORE" grd \
  "$OUT_DIR/wall_sub0.geo" "$OUT_DIR/wall_sub1.geo" \
  "$OUT_DIR/wall_sub2.geo" "$OUT_DIR/wall_sub3.geo" \
  "$OUT_DIR/wall_cap_north.geo" "$OUT_DIR/wall_cap_south.geo" \
  -o "$OUT_DIR/wall.grd"

"$SHORE" grd \
  "$OUT_DIR/background_center.geo" \
  "$OUT_DIR/background_sub_e.geo" "$OUT_DIR/background_sub_n.geo" \
  "$OUT_DIR/background_sub_w.geo" "$OUT_DIR/background_sub_s.geo" \
  -o "$OUT_DIR/background.grd"

# ── 8. Per-component cc.par ───────────────────────────────────────────────────
info "8. cc.par files"

# Wall cc.par: chimera handoff (free_face_bc=40) on every block, with
# the cap blocks promoted to overset level=1 to break ties cleanly.
# --boxes hole-cuts the background mesh's cells inside the body STL.
WALL_META='{
  "sub0":      {"free_face_bc": 40, "comment": "sub0"},
  "sub1":      {"free_face_bc": 40, "comment": "sub1"},
  "sub2":      {"free_face_bc": 40, "comment": "sub2"},
  "sub3":      {"free_face_bc": 40, "comment": "sub3"},
  "cap_north": {"free_face_bc": 40, "level": 1, "comment": "cap_north"},
  "cap_south": {"free_face_bc": 40, "level": 1, "comment": "cap_south"}
}'
"$SHORE" cc-par "$OUT_DIR/wall.adjacency.json" "$OUT_DIR/wall.cc.par" \
  --grd wall.grd \
  --meta "$WALL_META" \
  --boxes "$OUT_DIR/sphere_boxes.json"

# Background cc.par: far-field extrapolation everywhere (free_face_bc=9).
"$SHORE" cc-par "$OUT_DIR/background.adjacency.json" "$OUT_DIR/background.cc.par" \
  --grd background.grd \
  --free-bc 9

# ── 9. Per-component proc.input ───────────────────────────────────────────────
info "9. proc.input files (np=${NP})"

# Wall proc.input: every wall block belongs to physical body=1 (rigid
# body for force / moment integration).  Background blocks stay at
# body=0 by default.
WALL_PROC_META='{
  "1": {"body": 1}, "2": {"body": 1}, "3": {"body": 1},
  "4": {"body": 1}, "5": {"body": 1}, "6": {"body": 1}
}'
"$SHORE" proc-input "$OUT_DIR/wall.grd" "$OUT_DIR/wall.proc.input" \
  --np "$NP" --meta "$WALL_PROC_META"

"$SHORE" proc-input "$OUT_DIR/background.grd" "$OUT_DIR/background.proc.input" \
  --np "$NP"

# ── 10. Merge into a single Xall multi-mesh job ───────────────────────────────
# Order: background first (5 blocks → indices 1..5), wall second
# (6 blocks → indices 6..11).  All three merge calls must use the
# same ordering.
info "10. merge into assembly"
"$SHORE" grd-merge \
  "$OUT_DIR/background.grd" "$OUT_DIR/wall.grd" \
  -o "$OUT_DIR/assembly.grd"

"$SHORE" cc-par-merge \
  "$OUT_DIR/background.cc.par" "$OUT_DIR/wall.cc.par" \
  -o "$OUT_DIR/assembly.cc.par" \
  --grd assembly.grd

"$SHORE" proc-input-merge \
  "$OUT_DIR/background.proc.input" "$OUT_DIR/wall.proc.input" \
  -o "$OUT_DIR/assembly.proc.input" \
  --grd "$OUT_DIR/assembly.grd"

info "Done — outputs in $OUT_DIR"
echo
echo "Key files:"
ls -1 "$OUT_DIR"/{assembly.grd,assembly.cc.par,assembly.proc.input} 2>/dev/null || true
echo
echo "Inspect with:"
echo "  $SHORE info $OUT_DIR/assembly.grd"
echo "  paraview $OUT_DIR/wall.vtm $OUT_DIR/background.vtm $OUT_DIR/sphere_boxes.vtm"
