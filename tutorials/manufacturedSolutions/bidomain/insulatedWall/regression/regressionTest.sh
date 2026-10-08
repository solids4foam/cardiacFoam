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

# Compares the summary values and error norms of the manufactured-solution
# error summary written by the verifier against the reference.

regression_init "Bidomain insulated-wall manufactured-solution regression test" \
    regression/insulatedWall.reference "$@"

regression_run_or_fail parallel

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'anufactured-solution error summary' 'postProcessing/*.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" 'Field     L1-error'

regression_compare || true
regression_finish
