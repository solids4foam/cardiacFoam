#!/bin/bash
# Reproduce the Frontal-family rotated eikonal MMS matrix.
if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
  if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
    source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
  else
    echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
    exit 2
  fi
fi

set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"
AGGREGATE="$REPO_ROOT/applications/scripts/paperI_results/aggregate.py"
STUDY="$CASE_DIR/setup/studies/tetConvergence"

mkdir -p "$CASE_DIR/setup/mesh/tet"
cp "$REPO_ROOT/tutorials/manufacturedSolutions/monodomainPseudoECG/setup/mesh/tet/box.geo.template.optimised" "$CASE_DIR/setup/mesh/tet/box.geo.template"

rm -rf "$STUDY/results/sweepCasesFrontal" \
       "$STUDY/results/sweepRunFrontal"
mkdir -p "$STUDY/results"

"$DRIVER" sweep-run \
  --spec "$STUDY/sweep_tet_frontal.json" \
  --output-dir "$STUDY/results/sweepRunFrontal"

python3 "$AGGREGATE" eikonal_tet_frontal --repo-root "$REPO_ROOT"
