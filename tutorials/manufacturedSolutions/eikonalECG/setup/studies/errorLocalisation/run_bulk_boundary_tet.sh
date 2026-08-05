#!/bin/bash
# Canonical driverFOAM path for eikonalECG's solved-field bulk/boundary
# error decomposition across the tet N-ladder (both gradient schemes).
#
# Complementary to ../gradientVerification/run_error_localisation.sh, which
# is a single-case (default N=40, leastSquares) spatial-correlation deep
# dive kept for direct/manual use -- see that script's header. This one
# instead sweeps N=10,20,40,80 x {leastSquares,gaussLinear} through the
# manufacturedEikonalECG driverFOAM entry with writeErrorField enabled (see
# sweep_tet_error_localisation.json), to get the L2_bulk/L2_boundary/
# L2_total trend the paper reports, using the same bulk/boundary split
# convention (manufacturedEikonalVerifier.C's computeBoundaryBulkNorms) as
# the standalone gradientReconstructionOrder utility's own decomposition.
if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
  if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
    source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
  else
    echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
    exit 2
  fi
fi

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"

rm -rf "$SCRIPT_DIR/results/sweepCases" "$SCRIPT_DIR/results/sweepRun"
mkdir -p "$SCRIPT_DIR/results"
"$DRIVER" sweep-run \
  --spec "$SCRIPT_DIR/sweep_tet_error_localisation.json" \
  --output-dir "$SCRIPT_DIR/results/sweepRun"

python3 "$SCRIPT_DIR/aggregate_bulk_boundary.py"
