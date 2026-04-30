#!/usr/bin/env bash
# visualize_h_grid.sh — bash equivalent of visualize_h_grid.py.
#
# Builds two H-grid variants (rectangular duct + circular tube)
# using only `shore` CLI calls.  Outputs:
#
#     examples_out/
#         h_grid_rect.geo      # 1m-long rectangular duct
#         h_grid_rect.vtm      # ParaView MultiBlock (1 block)
#         h_grid_tube.geo      # 0.5m-long circular tube
#         h_grid_tube.vtm      # ParaView MultiBlock (1 block)
#
# Usage:
#   ./examples/07-h-grid/visualize_h_grid.sh             # default examples_out/
#   ./examples/07-h-grid/visualize_h_grid.sh -o /tmp/out

set -euo pipefail

OUT_DIR="examples_out"
SHORE="${SHORE:-shore}"

usage() {
  sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o | --output)  OUT_DIR="$2"; shift 2 ;;
    -h | --help)    usage 0 ;;
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

# ── 1. Rectangular duct ──────────────────────────────────────────────────────
info "1. Rectangular duct (1m × 0.1m × 0.2m, streamwise=x)"
"$SHORE" primitive h-grid-rect \
    --length 1.0 --width 0.1 --height 0.2 \
    --ni 80 --nj 20 --nk 24 \
    --streamwise-axis x \
    --j-spacing tanh2 --j-beta 4.0 \
    --k-spacing tanh2 --k-beta 4.0 \
    -o "$OUT_DIR/h_grid_rect"

"$SHORE" export "$OUT_DIR/h_grid_rect.geo" -o "$OUT_DIR/h_grid_rect.vtm" \
  || info "  (pyvista not installed; skipped h_grid_rect.vtm)"

# ── 2. Circular tube ─────────────────────────────────────────────────────────
info "2. Circular tube (0.5m long, radius 0.05m, streamwise=x)"
"$SHORE" primitive h-grid-tube \
    --length 0.5 --radius 0.05 \
    --n-axial 60 --n-circ 32 --n-radial 10 \
    --streamwise-axis x \
    --k-spacing tanh --k-beta 3.0 \
    -o "$OUT_DIR/h_grid_tube"

"$SHORE" export "$OUT_DIR/h_grid_tube.geo" -o "$OUT_DIR/h_grid_tube.vtm" \
  || info "  (pyvista not installed; skipped h_grid_tube.vtm)"

info "Done — outputs in $OUT_DIR"
echo
echo "Inspect with:"
echo "  $SHORE info $OUT_DIR/h_grid_rect.geo"
echo "  $SHORE info $OUT_DIR/h_grid_tube.geo"
echo "  paraview $OUT_DIR/h_grid_rect.vtm $OUT_DIR/h_grid_tube.vtm"
