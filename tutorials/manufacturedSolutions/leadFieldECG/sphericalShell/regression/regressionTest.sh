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

# The curved lead-field solution with the case defaults (hex cubed sphere,
# N = 12), run serially; the verifier writes
# postProcessing/manufacturedLeadFieldECGSummary_ECG.dat, whose reciprocity
# gap and static errors against the exact potential are compared against the
# reference.
regression_init "Lead-field ECG curved-shell manufactured-solution regression test" \
    regression/sphericalShell.reference "$@"

regression_run_or_fail

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Manufactured lead-field ECG summary' \
    'postProcessing/manufacturedLeadFieldECGSummary_*.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "lead-field ECG verifier summary" \
    "${REGRESSION_SUMMARY_FILE}" '# row L1 L2 Linf'

regression_compare || true
regression_finish
