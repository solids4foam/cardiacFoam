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

# The G1 level-set solution on the sheared slab (N = 20, configuration ATBI),
# run in parallel on four ranks with the manual decomposition so that
# interface faces are both internal to a rank and on processor patches. The
# verifier writes postProcessing/2D_<N>_cells.dat; its summary values and
# error norms are compared against the reference.
regression_init "Bath-bidomain oblique-slab transmission manufactured-solution regression test" \
    regression/obliqueSlab.reference "$@"

# regression/skips-on-openfoam: the regression is suppressed on OpenFOAM v2312,
# where the cellZone myocardium sub-mesh fails its fvSchemes lookups.
if [[ "${WM_PROJECT_VERSION:-unknown}" == *2312* ]]; then
    echo "SKIP: oblique-slab regression is suppressed on OpenFOAM v2312 (regression/skips-on-openfoam)."
    exit "${REGRESSION_SKIP_CODE}"
fi

regression_run_or_fail parallel -N 20 -endTime 0.02

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Bath-bidomain manufactured solution error summary' 'postProcessing/*_cells.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "bath-bidomain transmission manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

regression_compare || true
regression_finish
