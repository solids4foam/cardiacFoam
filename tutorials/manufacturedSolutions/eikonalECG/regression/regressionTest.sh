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

# Runs the case defaults (3D hex, rotated conductivity, N 20) in parallel and
# compares the activation-time error summary, the final ECG values of E1-E5
# and their manufactured ECG errors against the references.

regression_init "Eikonal ECG manufactured-solution regression test" \
    regression/eikonalECG.reference "$@"

regression_run_or_fail parallel

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Eikonal manufactured activation-time summary' \
    'postProcessing/manufacturedEikonalActivationTime.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "manufactured activation-time summary" \
    "${REGRESSION_SUMMARY_FILE}" activationTime

ecgSummary="$(regression_find_output 'Manufactured eikonal ECG summary' \
    'postProcessing/manufacturedEikonalECGSummary_*.dat')" || ecgSummary=""
regression_require_output "manufactured eikonal ECG summary" "${ecgSummary}"
ecgFile="$(regression_resolve_file postProcessing/eikonalECG.dat)" || ecgFile=""
regression_require_output "eikonalECG trace" "${ecgFile}" '^#.*time'

regression_compare || true

# The manufactured ECG errors against the reference integral
REGRESSION_SUMMARY_FILE="${ecgSummary}"
regression_compare regression/eikonalECG.ecg.reference || true
regression_finish
