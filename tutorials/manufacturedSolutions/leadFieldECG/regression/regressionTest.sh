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

# Runs the case defaults (slab, hex, N 10, monodomain) in parallel on two
# ranks and compares the rows of the lead-field ECG verifier summary
# postProcessing/manufacturedLeadFieldECGSummary_ECG.dat against the
# reference: the reciprocity gap leadField-direct and the static errors
# against the exact potential.

regression_init "Lead-field ECG manufactured-solution regression test" \
    regression/leadFieldECG.reference "$@"

regression_run_or_fail parallel

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Manufactured lead-field ECG summary' \
    'postProcessing/manufacturedLeadFieldECGSummary_*.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "lead-field ECG verifier summary" \
    "${REGRESSION_SUMMARY_FILE}" '# row L1 L2 Linf'

regression_compare || true
regression_finish
