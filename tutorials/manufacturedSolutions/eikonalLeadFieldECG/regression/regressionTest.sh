#!/usr/bin/env bash
set -euo pipefail
regressionLib="${CARDIAC_REGRESSION_LIB:-}"
if [[ -z "$regressionLib" ]]; then
    dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    while [[ "$dir" != / && ! -f "$dir/regression/lib.sh" ]]; do dir="$(dirname "$dir")"; done
    regressionLib="$dir/regression/lib.sh"
fi
. "$regressionLib"
set +u
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
set -u
regression_init "Eikonal dual ECG manufactured verification" regression/leadFieldECG.reference "$@"
regression_run_or_fail parallel
regression_require_output "reciprocal recording" postProcessing/eikonalLeadFieldECG.dat '^#.*time'
regression_require_output "independent boundary integral" postProcessing/manufacturedEikonalLeadFieldECG_LeadField.dat '^#.*time'
regression_require_output "concurrent pseudo recording" postProcessing/eikonalECG.dat '^#.*time'
regression_compare || true
if ! python3 "$(dirname "$regressionLib")/checkEikonalDualECG.py" localCase; then
    regression_fail "eikonal observer or independent reference control failed"
fi
for guard in regionalTemplate unknownTemplate invalidSampling; do
    case "$guard" in
        regionalTemplate) message='requires a spatially uniform template';;
        unknownTemplate) message='templateType must be';;
        invalidSampling) message='Invalid eikonalLeadFieldECG sampling interval';;
    esac
    if setManufacturedLeadFields -case localCase -eikonalGuard "$guard" \
        > "localCase/log.guard.$guard" 2>&1; then
        regression_fail "eikonal guard $guard accepted an invalid configuration"
    elif ! grep -q "$message" "localCase/log.guard.$guard"; then
        regression_fail "eikonal guard $guard failed for an unexpected reason"
    else
        echo "PASS: eikonal guard $guard"
    fi
done
regression_finish
