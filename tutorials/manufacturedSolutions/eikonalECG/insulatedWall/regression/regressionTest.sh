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

# Compares the activation-time error norms, the final ECG electrode values
# and sampled ECG trace points against the reference.

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

# The verifier writes one activation-time summary per convergence-study
# resolution. The file matching the mesh in system/blockMeshDict.3D is the
# one compared; the middle entry of the sorted set is the fallback.
findErrorSummary()
{
    local meshN candidate
    local -a fallbacks=()

    meshN="$(grep -E 'hex \(' system/blockMeshDict.3D \
        | grep -oE '\) \([0-9]+' | grep -oE '[0-9]+' | head -1 || true)"

    for candidate in postProcessing/*.dat processor*/postProcessing/*.dat; do
        [[ -s "${candidate}" ]] || continue
        grep -q 'Eikonal manufactured activation-time summary' "${candidate}" || continue
        grep -q 'dimension 3D' "${candidate}" || continue
        if [[ -n "${meshN}" && "${candidate}" == *"3D_${meshN}_cells_"* ]]; then
            echo "${candidate}"
            return 0
        fi
        fallbacks+=("${candidate}")
    done

    if (( ${#fallbacks[@]} > 0 )); then
        printf '%s\n' "${fallbacks[@]}" | sort -V | head -$(( ${#fallbacks[@]} / 2 + 1 )) | tail -1
        return 0
    fi
    return 1
}

regression_init "Eikonal ECG insulated-wall manufactured-solution regression test" \
    regression/insulatedWall.reference "$@"
REGRESSION_TIME_WINDOW=1e-6
checkElectrodeConfiguration || regression_finish

# The tutorial samples 161 electrodes, and the ECG verifier integrates a
# reference ECG for each by Gauss quadrature at every sample, which is most of
# the case's cost. The regression keeps E1-E5, the electrodes it checks; each
# electrode's trace is computed on its own, so theirs do not change.
regression_edit constant/electroProperties '/^[[:space:]]*R[0-9]+[[:space:]]*\(/d'

regression_run_or_fail parallel

REGRESSION_SUMMARY_FILE="$(findErrorSummary)" || REGRESSION_SUMMARY_FILE=""
regression_require_output "manufactured activation-time summary" \
    "${REGRESSION_SUMMARY_FILE}" activationTime
ecgSummary="$(regression_find_output 'Manufactured eikonal ECG summary' \
    'postProcessing/manufacturedEikonalECGSummary_*.dat')" || ecgSummary=""
regression_require_output "manufactured eikonal ECG summary" "${ecgSummary}"
ecgFile="$(regression_resolve_file postProcessing/eikonalECG.dat)" || ecgFile=""
regression_require_output "eikonalECG trace" "${ecgFile}" '^#.*time'

regression_compare || true
regression_finish
