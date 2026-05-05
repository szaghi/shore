#!/usr/bin/env bash
# visualize_ch_hybrid.sh — bash equivalent of visualize_ch_hybrid.py.
#
# Builds a CH hybrid mesh (C-grid wrap + H-channel far-field) around
# a NACA 0012 wing using only `shore` CLI calls (no Python imports).
# CH-3 (sharp TE, default) writes 12 .geo files; CH-4 (blunt TE, when
# the body's TE is detected as blunt) writes 14.  In either case the
# script exports a single ParaView .vtm bundling every block.
#
# Outputs (sharp-TE default, 12 blocks):
#
#     examples_out/
#         naca0012.dat                       # 2D Selig section (intermediate)
#         naca0012.stl                       # 3D wing STL
#         ch_upper.geo                       # wrap blocks (3)
#         ch_body.geo
#         ch_lower.geo
#         ch_chan_upper.geo                  # k_hi-side channels (5)
#         ch_chan_body_top.geo
#         ch_chan_body_front.geo
#         ch_chan_body_bot.geo
#         ch_chan_lower.geo
#         ch_chan_outflow_top_upper.geo      # outflow sub-blocks (4)
#         ch_chan_outflow_wake_upper.geo
#         ch_chan_outflow_wake_lower.geo
#         ch_chan_outflow_top_lower.geo
#         ch.vtm                             # ParaView MultiBlock
#
# Pipeline:
#   1. shore profile naca 0012     -> 2D Selig section .dat
#   2. shore body extrude          -> 3D wing STL with end caps
#   3. shore mesh --topology ch    -> 12 (sharp) or 14 (blunt) .geo files
#   4. shore export                -> ParaView .vtm
#
# Usage:
#   ./examples/09-ch-hybrid/visualize_ch_hybrid.sh                  # default examples_out/
#   ./examples/09-ch-hybrid/visualize_ch_hybrid.sh -o /tmp/out
#   ./examples/09-ch-hybrid/visualize_ch_hybrid.sh --ni-body 100 --nk 16
#   ./examples/09-ch-hybrid/visualize_ch_hybrid.sh --chord-min -3 --chord-max 18 \
#                                                  --normal-min -3.5 --normal-max 3.5 \
#                                                  --nk-channel 12 --k-law-channel tanh

set -euo pipefail

# ── Defaults (match visualize_ch_hybrid.py) ──────────────────────────────────
NI_BODY=80
NI_WAKE=20
N_STATIONS=5
WAKE_LENGTH=15.0
WAKE_GROWTH=1.15
SPAN_AXIS="z"
NK=12
DS=5e-4
GROWTH=1.08
SMOOTH=0.0
SMOOTH_ITERS=0

# CH far-field frame extents (chord/normal plane, span comes from the body).
CHORD_MIN=-2.0
CHORD_MAX=20.0
NORMAL_MIN=-2.5
NORMAL_MAX=2.5

# Channel knobs.
NK_CHANNEL=8
K_LAW_CHANNEL="uniform"
K_BETA_CHANNEL=3.0

# Wing-STL parameters (must straddle the C-grid slicer's planes).
STL_SPAN=0.5
STL_N_STATIONS=8
STL_N_PTS=40                 # chord samples per side (LE..TE)

OUT_DIR="examples_out"
SHORE="${SHORE:-shore}"

usage() {
  sed -n '2,33p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o | --output)        OUT_DIR="$2"; shift 2 ;;
    --ni-body)            NI_BODY="$2"; shift 2 ;;
    --ni-wake)            NI_WAKE="$2"; shift 2 ;;
    --n-stations)         N_STATIONS="$2"; shift 2 ;;
    --wake-length)        WAKE_LENGTH="$2"; shift 2 ;;
    --wake-growth)        WAKE_GROWTH="$2"; shift 2 ;;
    --span-axis)          SPAN_AXIS="$2"; shift 2 ;;
    --nk)                 NK="$2"; shift 2 ;;
    --ds)                 DS="$2"; shift 2 ;;
    --growth)             GROWTH="$2"; shift 2 ;;
    --smooth)             SMOOTH="$2"; shift 2 ;;
    --smooth-iters)       SMOOTH_ITERS="$2"; shift 2 ;;
    --chord-min)          CHORD_MIN="$2"; shift 2 ;;
    --chord-max)          CHORD_MAX="$2"; shift 2 ;;
    --normal-min)         NORMAL_MIN="$2"; shift 2 ;;
    --normal-max)         NORMAL_MAX="$2"; shift 2 ;;
    --nk-channel)         NK_CHANNEL="$2"; shift 2 ;;
    --k-law-channel)      K_LAW_CHANNEL="$2"; shift 2 ;;
    --k-beta-channel)     K_BETA_CHANNEL="$2"; shift 2 ;;
    -h | --help)          usage 0 ;;
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

# Generous STL station count so the C-grid slicer doesn't hit the
# centroid-fan end caps.
if [[ "$STL_N_STATIONS" -lt $((N_STATIONS + 2)) ]]; then
  STL_N_STATIONS=$((N_STATIONS + 2))
fi

# ── 1. 2D Selig section ──────────────────────────────────────────────────────
info "1. NACA 0012 section -> $OUT_DIR/naca0012.dat"
"$SHORE" profile naca 0012 \
  --n-pts "$STL_N_PTS" \
  -o "$OUT_DIR/naca0012.dat"

# ── 2. 3D wing STL ───────────────────────────────────────────────────────────
info "2. wing STL (span=$STL_SPAN, n_stations=$STL_N_STATIONS, --end-caps)"
"$SHORE" body extrude \
  --section "$OUT_DIR/naca0012.dat" \
  --span "$STL_SPAN" \
  --n-stations "$STL_N_STATIONS" \
  --span-axis "$SPAN_AXIS" \
  --end-caps \
  -o "$OUT_DIR/naca0012.stl"

# ── 3. CH hybrid mesh ────────────────────────────────────────────────────────
info "3. CH hybrid (ni_body=$NI_BODY, ni_wake=$NI_WAKE, n_stations=$N_STATIONS, "
info "             nk_wrap=$NK, nk_channel=$NK_CHANNEL, "
info "             frame chord=[$CHORD_MIN,$CHORD_MAX] normal=[$NORMAL_MIN,$NORMAL_MAX])"
"$SHORE" mesh "$OUT_DIR/naca0012.stl" \
  --topology ch \
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
  --ch-chord-min "$CHORD_MIN" \
  --ch-chord-max "$CHORD_MAX" \
  --ch-normal-min "$NORMAL_MIN" \
  --ch-normal-max "$NORMAL_MAX" \
  --nk-channel "$NK_CHANNEL" \
  --k-law-channel "$K_LAW_CHANNEL" \
  --k-beta-channel "$K_BETA_CHANNEL" \
  -o "$OUT_DIR/ch"

# ── 4. ParaView .vtm ─────────────────────────────────────────────────────────
# Bundle every ch_*.geo into a MultiBlock for ParaView.  Globbing picks
# up CH-3 (12 blocks) and CH-4 (14 blocks) without further branching.
info "4. ParaView .vtm bundling all ch_*.geo"
shopt -s nullglob
GEO_FILES=("$OUT_DIR"/ch_*.geo)
shopt -u nullglob
if [[ ${#GEO_FILES[@]} -eq 0 ]]; then
  err "no ch_*.geo files found in $OUT_DIR — mesh step failed silently."
fi
"$SHORE" export "${GEO_FILES[@]}" -o "$OUT_DIR/ch.vtm" \
  || info "  (pyvista not installed; skipped ch.vtm)"

info "Done — outputs in $OUT_DIR (${#GEO_FILES[@]} blocks)"
echo
echo "Inspect with:"
echo "  $SHORE info $OUT_DIR/ch_body.geo"
echo "  $SHORE check $OUT_DIR/ch_body.geo"
echo "  paraview $OUT_DIR/ch.vtm"
