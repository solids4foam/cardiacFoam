#!/bin/bash
# Reproduce the Frontal-family diagonal/rotated monodomain MMS matrix.
set -euo pipefail

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
  if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
    source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
  else
    echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
    exit 2
  fi
fi

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"
AGGREGATE="$REPO_ROOT/applications/scripts/paperI_results/aggregate.py"
STUDY="$CASE_DIR/setup/studies/tetConvergence"

rm -rf "$STUDY/results/sweepCasesFrontal" \
       "$STUDY/results/sweepRunFrontal"
mkdir -p "$STUDY/results"
"$DRIVER" sweep-run \
  --spec "$STUDY/sweep_tet_frontal.json" \
  --output-dir "$STUDY/results/sweepRunFrontal"

python3 "$AGGREGATE" mono_tet_frontal --repo-root "$REPO_ROOT"
