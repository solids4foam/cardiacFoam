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
# error summary and the final pseudo-ECG electrode values against the reference.

EXPECTED_ELECTRODES=(E1 E2 E3 E4 E5)

# constant/electroProperties must name every electrode the reference samples.
checkElectrodeConfiguration()
{
    local dict="constant/electroProperties" electrode

    if ! grep -q 'ecgDomains' "${dict}" || ! grep -q 'electrodePositions' "${dict}"; then
        regression_fail "pseudo-ECG electrode configuration not found in ${dict}"
        return 1
    fi
    for electrode in "${EXPECTED_ELECTRODES[@]}"; do
        if ! grep -Eq "^[[:space:]]*${electrode}[[:space:]]*\\(" "${dict}"; then
            regression_fail "electrode ${electrode} not found in ${dict}"
            return 1
        fi
    done
    echo "PASS: pseudo-ECG electrode configuration contains ${EXPECTED_ELECTRODES[*]}"
}

regression_init "Monodomain pseudo-ECG manufactured-solution regression test" \
    regression/monodomainPseudoECG.reference "$@"
checkElectrodeConfiguration || regression_finish

regression_run_or_fail parallel

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'anufactured-solution error summary' 'postProcessing/*.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" 'Field     L1-error'

ecgSummary="$(regression_find_output 'Manufactured pseudo-ECG summary' \
    'postProcessing/manufacturedPseudoECGSummary_*.dat')" || ecgSummary=""
regression_require_output "manufactured pseudo-ECG summary" "${ecgSummary}"
ecgFile="$(regression_resolve_file postProcessing/pseudoECG.dat)" || ecgFile=""
regression_require_output "pseudoECG trace" "${ecgFile}" '^#.*time'

regression_compare || true
regression_finish
