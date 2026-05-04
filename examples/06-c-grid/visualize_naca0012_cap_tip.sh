#!/usr/bin/env bash
# visualize_naca0012_cap_tip.sh — conformal C-grid + tip caps via shore CLI.
#
# Builds a C-grid mesh around a NACA 0012 wing, with **conformal**
# wing-tip caps at both spanwise extremes.  The C-grid is auto-promoted
# to its 4-block layout and the body's TE is rounded over its full
# span by ``--cap-te-thickness-frac`` so the cap-body interface is
# bit-coincident — no Chimera fringe needed at the cap-body seam.
#
# Outputs (default ``examples_out/``):
#
#     naca0012.dat                # 2D Selig section (intermediate)
#     naca0012.stl                # 3D wing STL
#     wing_upper.geo              # C-grid upper-wake block
#     wing_body.geo               # C-grid body wrap (rounded TE)
#     wing_lower.geo              # C-grid lower-wake block
#     wing_wake_base.geo          # C-grid wake-base block (TE bridge)
#     wing_cap_jlo_le.geo         # jlo cap blocks (6 of them)
#     wing_cap_jlo_upper.geo
#     wing_cap_jlo_lower.geo
#     wing_cap_jlo_te_upper.geo
#     wing_cap_jlo_te_lower.geo
#     wing_cap_jlo_wake_base.geo
#     wing_cap_jhi_*.geo          # jhi cap blocks (6, mirror)
#     wing_capped.vtm             # ParaView MultiBlock (16 blocks)
#
# Pipeline:
#   1. shore profile naca 0012   ->  2D Selig section .dat
#   2. shore body extrude        ->  3D wing STL with end caps
#   3. shore mesh --topology cgrid --cap-tips both
#                                ->  4 C-grid + 12 cap blocks (conformal)
#   4. shore export              ->  ParaView .vtm (16 blocks)
#
# Usage:
#   ./examples/06-c-grid/visualize_naca0012_cap_tip.sh
#   ./examples/06-c-grid/visualize_naca0012_cap_tip.sh -o /tmp/out
#   ./examples/06-c-grid/visualize_naca0012_cap_tip.sh --cap-te-thickness-frac 0.005

set -euo pipefail

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

# Cap parameters
CAP_TIPS="both"
CAP_LE_FRAC=0.05
CAP_TE_FRAC=0.95
CAP_TE_THICKNESS_FRAC=0.002
CAP_TIP_INSET=0.02
CAP_N_CAP_J=12

STL_SPAN=0.5
STL_N_STATIONS=8
STL_N_PTS=40

OUT_DIR="examples_out"
SHORE="${SHORE:-shore}"

usage() {
  sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o | --output)             OUT_DIR="$2"; shift 2 ;;
    --ni-body)                 NI_BODY="$2"; shift 2 ;;
    --ni-wake)                 NI_WAKE="$2"; shift 2 ;;
    --n-stations)              N_STATIONS="$2"; shift 2 ;;
    --wake-length)             WAKE_LENGTH="$2"; shift 2 ;;
    --nk)                      NK="$2"; shift 2 ;;
    --ds)                      DS="$2"; shift 2 ;;
    --growth)                  GROWTH="$2"; shift 2 ;;
    --cap-tips)                CAP_TIPS="$2"; shift 2 ;;
    --cap-le-frac)             CAP_LE_FRAC="$2"; shift 2 ;;
    --cap-te-frac)             CAP_TE_FRAC="$2"; shift 2 ;;
    --cap-te-thickness-frac)   CAP_TE_THICKNESS_FRAC="$2"; shift 2 ;;
    --cap-tip-inset)           CAP_TIP_INSET="$2"; shift 2 ;;
    --cap-n-cap-j)             CAP_N_CAP_J="$2"; shift 2 ;;
    -h | --help)               usage 0 ;;
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

if [[ "$STL_N_STATIONS" -lt $((N_STATIONS + 2)) ]]; then
  STL_N_STATIONS=$((N_STATIONS + 2))
fi

info "1. NACA 0012 section -> $OUT_DIR/naca0012.dat"
"$SHORE" profile naca 0012 --n-pts "$STL_N_PTS" -o "$OUT_DIR/naca0012.dat"

info "2. wing STL (span=$STL_SPAN, n_stations=$STL_N_STATIONS, --end-caps)"
"$SHORE" body extrude \
  --section "$OUT_DIR/naca0012.dat" \
  --span "$STL_SPAN" \
  --n-stations "$STL_N_STATIONS" \
  --span-axis "$SPAN_AXIS" \
  --end-caps \
  -o "$OUT_DIR/naca0012.stl"

info "3. C-grid + conformal caps (auto-rounded TE = $CAP_TE_THICKNESS_FRAC * chord)"
"$SHORE" mesh "$OUT_DIR/naca0012.stl" \
  --topology cgrid \
  --ni-body "$NI_BODY" --ni-wake "$NI_WAKE" \
  --n-stations "$N_STATIONS" --wake-length "$WAKE_LENGTH" \
  --wake-growth "$WAKE_GROWTH" --span-axis "$SPAN_AXIS" \
  --nk "$NK" --ds "$DS" --growth "$GROWTH" \
  --smooth "$SMOOTH" --smooth-iters "$SMOOTH_ITERS" \
  --cap-tips "$CAP_TIPS" \
  --cap-le-frac "$CAP_LE_FRAC" \
  --cap-te-frac "$CAP_TE_FRAC" \
  --cap-te-thickness-frac "$CAP_TE_THICKNESS_FRAC" \
  --cap-tip-inset "$CAP_TIP_INSET" \
  --cap-n-cap-j "$CAP_N_CAP_J" \
  -o "$OUT_DIR/wing"

info "4. ParaView .vtm (16 blocks: 4 C-grid + 6 cap_jlo + 6 cap_jhi)"
"$SHORE" export \
  "$OUT_DIR/wing_upper.geo" \
  "$OUT_DIR/wing_body.geo" \
  "$OUT_DIR/wing_lower.geo" \
  "$OUT_DIR/wing_wake_base.geo" \
  "$OUT_DIR/wing_cap_jlo_le.geo" \
  "$OUT_DIR/wing_cap_jlo_upper.geo" \
  "$OUT_DIR/wing_cap_jlo_lower.geo" \
  "$OUT_DIR/wing_cap_jlo_te_upper.geo" \
  "$OUT_DIR/wing_cap_jlo_te_lower.geo" \
  "$OUT_DIR/wing_cap_jlo_wake_base.geo" \
  "$OUT_DIR/wing_cap_jhi_le.geo" \
  "$OUT_DIR/wing_cap_jhi_upper.geo" \
  "$OUT_DIR/wing_cap_jhi_lower.geo" \
  "$OUT_DIR/wing_cap_jhi_te_upper.geo" \
  "$OUT_DIR/wing_cap_jhi_te_lower.geo" \
  "$OUT_DIR/wing_cap_jhi_wake_base.geo" \
  -o "$OUT_DIR/wing_capped.vtm" \
  || info "  (pyvista not installed; skipped wing_capped.vtm)"

info "Done — outputs in $OUT_DIR"
echo
echo "Inspect with:"
echo "  $SHORE info $OUT_DIR/wing_body.geo            # body block (rounded TE)"
echo "  $SHORE check $OUT_DIR/wing_cap_jlo_te_upper.geo"
echo "  paraview $OUT_DIR/wing_capped.vtm"
