#!/bin/bash
# Canonical driverFOAM path for the axis-aligned FDA tet pseudo-ECG study,
# broken out per electrode rather than collapsed to max/mean/min.
#
# Supersedes the hand-rolled run_per_electrode_tet.sh in this same
# directory: that script called the tet ladder directly with Allclean
# suppressed so postProcessing/manufacturedPseudoECGSummary.dat survived to
# the next N, and pointed ELECTRO_PROPERTIES at
# configs/electroProperties.fdaDiagonalECG to get the axis-aligned FDA
# verifier and conductivity. Both workarounds are unnecessary here:
#   - driverFOAM's generic sweepCases/<case_id>/ archive (openfoam_driver's
#     snapshot/diff collector) already preserves each case's own
#     postProcessing/ output before the next case's Allclean runs, so
#     nothing needs to be rescued manually.
#   - the manufacturedFDA entry's own defaults already carry the same
#     verifier (manufacturedFDAMonodomainVerifier), quadrature settings, and
#     five electrode positions the deleted configs/ file duplicated; this
#     sweep's own "conductivity"/"verification_model_type" base values (see
#     sweep_tet_per_electrode.json) are the only override actually needed.
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

rm -rf "$STUDY/results/sweepCasesPerElectrode" \
       "$STUDY/results/sweepRunPerElectrode"
mkdir -p "$STUDY/results"
"$DRIVER" sweep-run \
  --spec "$STUDY/sweep_tet_per_electrode.json" \
  --output-dir "$STUDY/results/sweepRunPerElectrode"

python3 "$AGGREGATE" mono_tet_per_electrode --repo-root "$REPO_ROOT"
