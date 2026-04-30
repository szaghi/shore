#!/usr/bin/env bash
# release.sh — create a new SHORE release
#
# Usage:
#   ./release.sh (--major | --minor | --patch | <X.Y.Z>)
#
#   --major, -M     X.Y.Z → X+1.0.0
#   --minor, -m     X.Y.Z → X.Y+1.0
#   --patch, -p     X.Y.Z → X.Y.Z+1
#   <X.Y.Z>         set an explicit version
#
# What it does:
#   1. Pre-flight: clean tree, on main, up-to-date with remote, tag free
#   2. Lint gate (ruff check + format --check)
#   3. Confirms with user
#   4. Bumps version in pyproject.toml + src/shore/__init__.py
#   5. Updates CHANGELOG.md and docs/guide/changelog.md
#       - First release (no prior v* tag): rewrite the manual
#         "## [Unreleased]" header in docs/guide/changelog.md to
#         "## [X.Y.Z] — YYYY-MM-DD" and seed CHANGELOG.md from it,
#         preserving the curated prose.
#       - Subsequent releases: regenerate both via git-cliff from
#         conventional commits.
#   6. Runs test suite
#   7. Commits + annotated tag + push (--follow-tags → triggers CI → PyPI)

set -euo pipefail

TRUNK="main"
INIT_FILE="src/shore/__init__.py"
DOCS_CHANGELOG="docs/guide/changelog.md"
ROOT_CHANGELOG="CHANGELOG.md"

die()  { echo "ERROR: $*" >&2; exit 1; }
info() { echo "==> $*"; }

usage() {
  sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

current_version() {
  grep -E '^version = ' pyproject.toml | head -1 | sed 's/version = "\(.*\)"/\1/'
}

bump() {
  local ver="$1" part="$2"
  local major minor patch
  IFS='.' read -r major minor patch <<< "$ver"
  case "$part" in
    major) echo "$((major + 1)).0.0" ;;
    minor) echo "${major}.$((minor + 1)).0" ;;
    patch) echo "${major}.${minor}.$((patch + 1))" ;;
  esac
}

# ── stage tracking + recovery trap ────────────────────────────────────────────
STAGE="preflight"
NEW_VER=""

on_error() {
  echo ""
  echo "================================================================"
  echo "  release.sh FAILED at stage: $STAGE"
  echo "================================================================"
  case "$STAGE" in
    preflight | lint | confirm)
      echo "  Nothing was changed. Fix the issue above and re-run."
      ;;
    bumped)
      echo "  Version bump and changelog were modified locally but not committed."
      echo "  To discard and start over:"
      echo "    git checkout -- pyproject.toml ${INIT_FILE} ${ROOT_CHANGELOG} ${DOCS_CHANGELOG}"
      ;;
    committed)
      echo "  Version bump was committed but not tagged/pushed. To resume:"
      echo "    git tag -a v${NEW_VER} -m \"Release v${NEW_VER}\""
      echo "    git push origin ${TRUNK} --follow-tags"
      ;;
    tagged)
      echo "  Tag v${NEW_VER} was created locally but not pushed. To resume:"
      echo "    git push origin ${TRUNK} --follow-tags"
      ;;
  esac
  echo "================================================================"
}

trap 'on_error' ERR

# ── argument parsing ───────────────────────────────────────────────────────────
[[ $# -ge 1 ]] || usage

BUMP_ARG=""

for arg in "$@"; do
  case "$arg" in
    --major | -M)          BUMP_ARG=major ;;
    --minor | -m)          BUMP_ARG=minor ;;
    --patch | -p)          BUMP_ARG=patch ;;
    [0-9]*.[0-9]*.[0-9]*)  BUMP_ARG="$arg" ;;
    *) usage ;;
  esac
done

[[ -n "$BUMP_ARG" ]] || usage

CUR_VER=$(current_version)
[[ -n "$CUR_VER" ]] || die "could not read current version from pyproject.toml"

case "$BUMP_ARG" in
  major | minor | patch) NEW_VER=$(bump "$CUR_VER" "$BUMP_ARG") ;;
  *)                     NEW_VER="$BUMP_ARG" ;;
esac

# ── pre-flight checks ──────────────────────────────────────────────────────────
STAGE="preflight"
[[ -f "$INIT_FILE" ]]           || die "$INIT_FILE not found — run from the repo root"
command -v git-cliff >/dev/null || die "'git-cliff' not found (install: pipx install git-cliff)"

CURRENT_BRANCH=$(git symbolic-ref --short HEAD 2>/dev/null || true)
[[ "$CURRENT_BRANCH" == "$TRUNK" ]] \
  || die "must be on '${TRUNK}' branch (currently on '${CURRENT_BRANCH}')"

[[ -z "$(git status --porcelain)" ]] \
  || die "working tree is not clean — commit or stash changes first"

git fetch --tags --quiet
[[ -z "$(git tag -l "v${NEW_VER}")" ]] || die "tag v${NEW_VER} already exists"

BEHIND=$(git rev-list --count HEAD..origin/${TRUNK} 2>/dev/null || echo 0)
[[ "$BEHIND" -eq 0 ]] \
  || die "${TRUNK} is ${BEHIND} commit(s) behind origin/${TRUNK} — run: git pull origin ${TRUNK}"

# Detect first release: no v* tag exists yet.
FIRST_RELEASE=0
if [[ -z "$(git tag -l 'v[0-9]*')" ]]; then
  FIRST_RELEASE=1
fi

# ── lint gate ──────────────────────────────────────────────────────────────────
STAGE="lint"
info "Running ruff lint + format check"
.venv/bin/ruff check src/ tests/         || die "lint failed — run 'make fmt' to fix, then retry"
.venv/bin/ruff format --check src/ tests/ || die "format check failed — run 'make fmt' to fix, then retry"

# ── confirm ────────────────────────────────────────────────────────────────────
STAGE="confirm"
echo "  Current version : $CUR_VER"
echo "  New version     : $NEW_VER"
if [[ "$FIRST_RELEASE" -eq 1 ]]; then
  echo "  Mode            : first release — preserving curated [Unreleased] prose"
else
  echo "  Mode            : regenerating changelog from conventional commits"
fi
echo
read -r -p "Proceed? [y/N] " confirm
[[ "$confirm" =~ ^[Yy]$ ]] || { echo "Aborted."; exit 0; }

# ── bump version ───────────────────────────────────────────────────────────────
STAGE="bumped"
if [[ "$CUR_VER" == "$NEW_VER" ]]; then
  if [[ "$FIRST_RELEASE" -eq 1 ]]; then
    info "Version unchanged ($CUR_VER) — first release, tagging existing version"
  else
    die "version is already ${NEW_VER} — nothing to bump"
  fi
else
  info "Bumping version ($CUR_VER → $NEW_VER)"

  sed -i "s/^version = \"${CUR_VER}\"/version = \"${NEW_VER}\"/" pyproject.toml
  grep -q "version = \"${NEW_VER}\"" pyproject.toml || die "pyproject.toml was not updated"

  sed -i "s/__version__ = \"${CUR_VER}\"/__version__ = \"${NEW_VER}\"/" "$INIT_FILE"
  grep -q "__version__ = \"${NEW_VER}\"" "$INIT_FILE" || die "$INIT_FILE was not updated"
fi

# ── update changelogs ──────────────────────────────────────────────────────────
TODAY=$(date +%Y-%m-%d)

if [[ "$FIRST_RELEASE" -eq 1 ]]; then
  info "First release — converting [Unreleased] header to [${NEW_VER}] in ${DOCS_CHANGELOG}"
  grep -q '^## \[Unreleased\]' "$DOCS_CHANGELOG" \
    || die "${DOCS_CHANGELOG} has no '## [Unreleased]' header to rewrite"
  sed -i "s/^## \[Unreleased\]\$/## [${NEW_VER}] — ${TODAY}/" "$DOCS_CHANGELOG"

  info "Seeding ${ROOT_CHANGELOG} from ${DOCS_CHANGELOG}"
  {
    printf '# Changelog\n\n'
    printf 'All notable changes to SHORE are documented here.\n'
    printf 'Versions follow [Semantic Versioning](https://semver.org/).\n'
    printf 'Format follows [Keep a Changelog](https://keepachangelog.com/).\n\n'
    awk '/^## \[/{found=1} found' "$DOCS_CHANGELOG"
  } > "$ROOT_CHANGELOG"
else
  info "Regenerating ${ROOT_CHANGELOG} via git-cliff for v${NEW_VER}"
  git cliff --tag "v${NEW_VER}" -o "$ROOT_CHANGELOG"

  info "Mirroring to ${DOCS_CHANGELOG} (with VitePress front-matter)"
  {
    printf -- "---\ntitle: Changelog\n---\n\n"
    awk '/^## \[/{found=1} found' "$ROOT_CHANGELOG"
  } > "$DOCS_CHANGELOG"
fi

# ── run tests ──────────────────────────────────────────────────────────────────
info "Running test suite"
.venv/bin/pytest

# ── commit, tag, push ──────────────────────────────────────────────────────────
STAGE="committed"
info "Committing version bump"
git add pyproject.toml "$INIT_FILE" "$ROOT_CHANGELOG" "$DOCS_CHANGELOG"
git commit -m "chore(release): bump version to ${NEW_VER}"

STAGE="tagged"
git tag -a "v${NEW_VER}" -m "Release v${NEW_VER}"

info "Pushing ${TRUNK} + tag to origin (CI will run tests and publish to PyPI)"
git push origin "${TRUNK}" --follow-tags

info "Done — v${NEW_VER} released"
