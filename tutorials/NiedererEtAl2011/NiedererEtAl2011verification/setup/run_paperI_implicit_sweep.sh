#!/bin/bash
#
# Re-run the Paper I Niederer spatial ladder on the IMPLICIT algorithm.
#
# The published dx = 0.2 mm and dx = 0.1 mm results were produced with
# solutionAlgorithm = explicit, while dx = 0.5 mm used implicit. This sweep
# regenerates the two finer levels implicitly so that the reported spatial
# ladder uses a single time-integration path throughout, matching the
# implicit MMS verification.
#
# Each config also carries dt = 0.01, 0.005, 0.001 ms, so the temporal
# sensitivity is measured at all three mesh spacings rather than at 0.5 mm
# only.
#
# Usage (OpenFOAM environment must be sourced first):
#   source /Volumes/OpenFOAM-v2412/etc/bashrc
#   ./run_paperI_implicit_sweep.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"

if [ -z "${WM_PROJECT_VERSION:-}" ]; then
    echo "ERROR: OpenFOAM environment not sourced (WM_PROJECT_VERSION unset)." >&2
    echo "       source /Volumes/OpenFOAM-v2412/etc/bashrc" >&2
    exit 3
fi

# Sourcing RunFunctions keeps the OpenFOAM library paths visible to child
# processes; without it cardiacFoam aborts on macOS with a missing
# libOpenFOAM.dylib under SIP.
#
# RunFunctions dereferences DYLD_LIBRARY_PATH, which is unset on a clean macOS
# shell because System Integrity Protection strips it. Under "set -u" that is a
# fatal unbound-variable error, so seed it and relax nounset across the source.
export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}"
set +u
# shellcheck disable=SC1091
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
set -u

echo "== Paper I Niederer implicit re-run =="
echo "   OpenFOAM : $WM_PROJECT_VERSION"
echo "   solver   : $(command -v cardiacFoam || echo 'NOT FOUND')"
echo "   started  : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo

for config in \
    paperI_dx02_implicit.json \
    paperI_dx01_implicit.json
do
    echo "== $config =="
    "$DRIVER" all \
        --entry niederer2012 \
        --config "$SCRIPT_DIR/$config"
    echo "== $config done at $(date -u '+%Y-%m-%dT%H:%M:%SZ') =="
    echo
done

echo "Paper I Niederer implicit sweep complete."
