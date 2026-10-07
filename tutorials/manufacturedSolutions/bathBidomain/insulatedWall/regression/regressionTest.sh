#!/usr/bin/env bash
set -euo pipefail

# Shared regression library (tutorials/regression/lib.sh), found by walking up
# from this case; CARDIAC_REGRESSION_LIB overrides the lookup.
regressionLib="${CARDIAC_REGRESSION_LIB:-}"
if [[ -z "${regressionLib}" ]]; then
    dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    while [[ "${dir}" != / && ! -f "${dir}/regression/lib.sh" ]]; do
        dir="$(dirname "${dir}")"
    done
    regressionLib="${dir}/regression/lib.sh"
fi
. "${regressionLib}"

# The 2D mesh (system/blockMeshDict, N = 20) with sealedHeartBoundary true
# and bathHeartPhiETrace global, matching constant/electroProperties. The
# verifier writes postProcessing/2D_<N>_cells.dat; its summary values and
# error norms are compared against the reference.
regression_init "Bath-bidomain insulated-wall manufactured-solution regression test" \
    regression/insulatedWall.reference "$@"

regression_run_or_fail

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Bath-bidomain manufactured solution error summary' 'postProcessing/*_cells.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "bath-bidomain manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

regression_compare || true
regression_finish
