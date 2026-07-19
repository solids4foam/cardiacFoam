#!/usr/bin/env bash
# Capture reproducibility provenance for a case directory as JSON on stdout.
# Pure bash + git + coreutils. Reads WM_PROJECT_VERSION (OpenFOAM env) if sourced.
set -euo pipefail

CASE_DIR="${1:?usage: capture_provenance.sh <caseDir>}"
REPO_ROOT="$(git -C "$CASE_DIR" rev-parse --show-toplevel 2>/dev/null || echo "")"

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
  elif command -v shasum   >/dev/null 2>&1; then shasum -a 256 "$1" | awk '{print $1}'
  else echo "no-sha256-tool"; fi
}
hash_or_absent() { [ -f "$1" ] && printf 'sha256:%s' "$(sha256 "$1")" || printf 'absent'; }
git_sha() { git -C "$1" rev-parse HEAD 2>/dev/null || echo "absent"; }
esc() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }

OF_VERSION="${WM_PROJECT_VERSION:-absent}"
CF_SHA="$( [ -n "$REPO_ROOT" ] && git_sha "$REPO_ROOT" || echo absent )"
S4F_SHA="absent"; [ -n "$REPO_ROOT" ] && [ -d "$REPO_ROOT/modules/solids4foam" ] && S4F_SHA="$(git_sha "$REPO_ROOT/modules/solids4foam")"
CC_SHA="absent"; [ -n "${CARDIACCORE_ROOT:-}" ] && CC_SHA="$(git_sha "$CARDIACCORE_ROOT")"
PY_VERSION="$(python3 --version 2>&1 || echo absent)"
HOSTID="$(hostname 2>/dev/null || echo absent)"
OSID="$(uname -a 2>/dev/null || echo absent)"

EP="$(hash_or_absent "$CASE_DIR/constant/electroProperties")"
CD="$(hash_or_absent "$CASE_DIR/system/controlDict")"
FS="$(hash_or_absent "$CASE_DIR/system/fvSchemes")"
FV="$(hash_or_absent "$CASE_DIR/system/fvSolution")"
MESH_FILE="$CASE_DIR/system/blockMeshDict"; [ -f "$MESH_FILE" ] || MESH_FILE="$CASE_DIR/constant/polyMesh/points"
MESH="$(hash_or_absent "$MESH_FILE")"

cat <<EOF
{
  "environment": {
    "openfoam_version": "$(esc "$OF_VERSION")",
    "cardiacfoam_sha": "$(esc "$CF_SHA")",
    "solids4foam_sha": "$(esc "$S4F_SHA")",
    "cardiaccore_sha": "$(esc "$CC_SHA")",
    "python_version": "$(esc "$PY_VERSION")",
    "host": "$(esc "$HOSTID")",
    "os": "$(esc "$OSID")"
  },
  "inputs": {
    "electroProperties": "$EP",
    "controlDict": "$CD",
    "fvSchemes": "$FS",
    "fvSolution": "$FV",
    "mesh": "$MESH"
  },
  "provenance": {
    "command": "$(esc "${PROVENANCE_CMD:-capture_provenance.sh $CASE_DIR}")",
    "timestamp_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "exit_status": "$(esc "${PROVENANCE_EXIT:-0}")"
  }
}
EOF
