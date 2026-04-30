#!/usr/bin/env bash
# visualize_naca0012.sh — bash equivalent of visualize_naca0012.py.
#
# Builds a C-grid mesh around a NACA 0012 wing using only `shore`
# CLI calls (no Python imports).  Produces the same outputs as the
# Python script:
#
#     examples_out/
#         naca0012.dat               # 2D Selig section (intermediate)
#         naca0012.stl               # 3D wing STL
#         naca0012.geo               # marched C-grid block
#         naca0012.vtm               # ParaView MultiBlock
#
# Pipeline:
#   1. shore profile naca 0012   ->  2D Selig section .dat
#   2. shore body extrude        ->  3D wing STL with end caps
#   3. shore mesh --topology cgrid -> single-block C-grid .geo
#   4. shore export              ->  ParaView .vtm
#
# Usage:
#   ./examples/06-c-grid/visualize_naca0012.sh                 # default examples_out/
#   ./examples/06-c-grid/visualize_naca0012.sh -o /tmp/out
#   ./examples/06-c-grid/visualize_naca0012.sh --ni-body 160 --nk 30

set -euo pipefail

# ── Defaults (match visualize_naca0012.py) ────────────────────────────────────
NI_BODY=120
NI_WAKE=30
N_STATIONS=5
WAKE_LENGTH=15.0
WAKE_GROWTH=1.15
SPAN_AXIS="z"
NK=20
DS=5e-4
GROWTH=1.15
SMOOTH=0.1
SMOOTH_ITERS=1

# Wing-STL parameters: span and section density.  The STL needs more
# span stations than the C-grid slicer (so the slice planes land
# between STL stations, not on the centroid-fan caps).
STL_SPAN=0.5
STL_N_STATIONS=8
STL_N_PTS=40                 # chord samples per side (LE..TE)

OUT_DIR="examples_out"
SHORE="${SHORE:-shore}"

usage() {
  sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o | --output)   OUT_DIR="$2"; shift 2 ;;
    --ni-body)       NI_BODY="$2"; shift 2 ;;
    --ni-wake)       NI_WAKE="$2"; shift 2 ;;
    --n-stations)    N_STATIONS="$2"; shift 2 ;;
    --wake-length)   WAKE_LENGTH="$2"; shift 2 ;;
    --wake-growth)   WAKE_GROWTH="$2"; shift 2 ;;
    --span-axis)     SPAN_AXIS="$2"; shift 2 ;;
    --nk)            NK="$2"; shift 2 ;;
    --ds)            DS="$2"; shift 2 ;;
    --growth)        GROWTH="$2"; shift 2 ;;
    --smooth)        SMOOTH="$2"; shift 2 ;;
    --smooth-iters)  SMOOTH_ITERS="$2"; shift 2 ;;
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

# Make the STL span-station count generous enough that the C-grid
# slicer (which places n_stations planes between the body's span
# extremes with a 2% margin) doesn't graze the centroid-fan caps.
if [[ "$STL_N_STATIONS" -lt $((N_STATIONS + 2)) ]]; then
  STL_N_STATIONS=$((N_STATIONS + 2))
fi

# ── 1. 2D Selig section ──────────────────────────────────────────────────────
info "1. NACA 0012 section -> $OUT_DIR/naca0012.dat"
"$SHORE" profile naca 0012 \
  --n-pts "$STL_N_PTS" \
  -o "$OUT_DIR/naca0012.dat"

# ── 2. 3D wing STL (extruded along --span-axis with end caps) ────────────────
info "2. wing STL (span=$STL_SPAN, n_stations=$STL_N_STATIONS, --end-caps)"
"$SHORE" body extrude \
  --section "$OUT_DIR/naca0012.dat" \
  --span "$STL_SPAN" \
  --n-stations "$STL_N_STATIONS" \
  --span-axis "$SPAN_AXIS" \
  --end-caps \
  -o "$OUT_DIR/naca0012.stl"

# ── 3. C-grid mesh ───────────────────────────────────────────────────────────
info "3. C-grid (ni_body=$NI_BODY, ni_wake=$NI_WAKE, n_stations=$N_STATIONS, nk=$NK)"
"$SHORE" mesh "$OUT_DIR/naca0012.stl" \
  --topology cgrid \
  --ni-body "$NI_BODY" \
  --ni-wake "$NI_WAKE" \
  --n-stations "$N_STATIONS" \
  --wake-length "$WAKE_LENGTH" \
  --wake-growth "$WAKE_GROWTH" \
  --span-axis "$SPAN_AXIS" \
  --nk "$NK" \
  --ds "$DS" \
  --growth "$GROWTH" \
  --smooth "$SMOOTH" \
  --smooth-iters "$SMOOTH_ITERS" \
  -o "$OUT_DIR/naca0012"

# ── 4. ParaView .vtm ─────────────────────────────────────────────────────────
info "4. ParaView .vtm"
"$SHORE" export "$OUT_DIR/naca0012.geo" -o "$OUT_DIR/naca0012.vtm" \
  || info "  (pyvista not installed; skipped naca0012.vtm)"

info "Done — outputs in $OUT_DIR"
echo
echo "Inspect with:"
echo "  $SHORE info $OUT_DIR/naca0012.geo"
echo "  $SHORE check $OUT_DIR/naca0012.geo"
echo "  paraview $OUT_DIR/naca0012.vtm"
