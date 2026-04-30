#!/usr/bin/env bash
# install.sh — install shore-mesh on Linux / macOS / WSL.
#
# Usage:
#   curl -fsSL https://raw.githubusercontent.com/szaghi/shore/main/install.sh | bash
#   curl -fsSL .../install.sh | bash -s -- --vis           # with pyvista extras
#   curl -fsSL .../install.sh | bash -s -- --version 0.1.0 # pin a version
#   curl -fsSL .../install.sh | bash -s -- --help
#
# What it does:
#   1. Verifies Python >= 3.11 is available.
#   2. Installs pipx if missing (apt on Debian/Ubuntu, pip elsewhere —
#      handles PEP 668 "externally managed" environments).
#   3. Runs `pipx install shore-mesh[<extras>]`.
#   4. Verifies the `shore` CLI is on PATH and prints next steps.
#
# This script never touches your system Python.  pipx puts shore-mesh
# in its own venv under ~/.local/share/pipx/ and symlinks the `shore`
# entry point to ~/.local/bin/shore.

set -euo pipefail

PKG="shore-mesh"
EXTRAS=""
PIN=""

usage() {
  sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vis)        EXTRAS="[vis]" ;;
    --version)    [[ $# -ge 2 ]] || { echo "--version requires an argument" >&2; exit 1; }
                  PIN="==$2"; shift ;;
    -h | --help)  usage 0 ;;
    *)            echo "Unknown option: $1" >&2; usage 1 ;;
  esac
  shift
done

cyan()  { printf '\033[1;36m%s\033[0m\n' "$*"; }
red()   { printf '\033[1;31m%s\033[0m\n' "$*" >&2; }
info()  { cyan "==> $*"; }
err()   { red "ERROR: $*"; exit 1; }

# ── 1. Python detection (>= 3.11 required) ────────────────────────────────────
detect_python() {
  for p in python3.13 python3.12 python3.11 python3; do
    if command -v "$p" >/dev/null 2>&1; then
      ver=$("$p" -c 'import sys; print("%d.%d" % sys.version_info[:2])' 2>/dev/null) || continue
      maj=${ver%.*}; min=${ver#*.}
      if [[ "$maj" -ge 3 ]] && [[ "$min" -ge 11 ]]; then
        echo "$p"; return 0
      fi
    fi
  done
  return 1
}

PY=$(detect_python) || err "Python >= 3.11 not found.
Install Python first (https://www.python.org/downloads/) and re-run."
info "Found Python: $($PY --version) at $(command -v "$PY")"

# ── 2. pipx — install if missing ──────────────────────────────────────────────
ensure_pipx() {
  if command -v pipx >/dev/null 2>&1; then
    info "pipx is already installed: $(pipx --version)"
    return 0
  fi

  info "pipx not found — installing it"

  # Detect Debian/Ubuntu: their system Python is PEP 668 externally
  # managed, so `pip install --user pipx` fails.  Use apt when we can.
  if [[ -r /etc/os-release ]]; then
    # shellcheck disable=SC1091
    . /etc/os-release
    case "${ID:-}${ID_LIKE:-}" in
      *debian*|*ubuntu*)
        if command -v apt-get >/dev/null 2>&1; then
          info "Detected Debian/Ubuntu — installing pipx via apt"
          if [[ "$(id -u)" -eq 0 ]]; then
            apt-get update -qq && apt-get install -y -qq pipx
          else
            sudo apt-get update -qq && sudo apt-get install -y -qq pipx
          fi
          return 0
        fi
        ;;
    esac
  fi

  # Fallback: pip --user.  Works on macOS, WSL, Fedora, Arch, etc.
  # On a PEP 668 system without apt available, the user gets a clear
  # error message from pip with instructions.
  "$PY" -m pip install --user --quiet pipx
  "$PY" -m pipx ensurepath >/dev/null 2>&1 || true
}

ensure_pipx

# pipx ensurepath modifies ~/.bashrc / ~/.zshrc but won't affect the
# *current* shell.  Add the user-bin dir to PATH for this session so
# the install + verify steps see it.
export PATH="${HOME}/.local/bin:${PATH}"

command -v pipx >/dev/null 2>&1 \
  || err "pipx is still not on PATH after install.
Open a new shell and re-run, or add this line to your shell profile:
    export PATH=\"\${HOME}/.local/bin:\${PATH}\""

# ── 3. Install shore-mesh ─────────────────────────────────────────────────────
SPEC="${PKG}${EXTRAS}${PIN}"
info "Installing ${SPEC} via pipx"
pipx install --force "${SPEC}"

# ── 4. Verify ─────────────────────────────────────────────────────────────────
if ! command -v shore >/dev/null 2>&1; then
  err "Installation completed but 'shore' is not on PATH.
Add this to your shell profile and reload:
    export PATH=\"\${HOME}/.local/bin:\${PATH}\""
fi

INSTALLED_VER=$(shore --version 2>/dev/null || echo "(installed)")
info "shore ${INSTALLED_VER} is ready"

cat <<'EOF'

Next steps:
  shore --help                # CLI overview
  shore mesh --help           # body-fitted meshing
  shore primitive --help      # parametric backgrounds

Docs:    https://szaghi.github.io/shore/
Issues:  https://github.com/szaghi/shore/issues

To upgrade later:  pipx upgrade shore-mesh
To remove:         pipx uninstall shore-mesh
EOF
