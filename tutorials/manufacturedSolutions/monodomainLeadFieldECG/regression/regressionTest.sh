#!/usr/bin/env bash
set -euo pipefail
regressionLib="${CARDIAC_REGRESSION_LIB:-}"
if [[ -z "$regressionLib" ]]; then
    dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    while [[ "$dir" != / && ! -f "$dir/regression/lib.sh" ]]; do dir="$(dirname "$dir")"; done
    regressionLib="$dir/regression/lib.sh"
fi
. "$regressionLib"
regression_init "Manufactured lead-field ECG" regression/leadFieldECG.reference "$@"
regression_run_or_fail parallel -pseudoECG yes
regression_require_output "lead-field trace" postProcessing/leadFieldECG.dat '^#.*time'
regression_require_output "independent functional" postProcessing/manufacturedLeadFieldECG_ECG.dat '^#.*time'
regression_compare || true
REGRESSION_SUMMARY_FILE="$(regression_find_output 'anufactured-solution error summary' 'postProcessing/*.dat')" || REGRESSION_SUMMARY_FILE=""
regression_require_output "source field oracle" "$REGRESSION_SUMMARY_FILE" 'Field     L1-error'
regression_compare regression/source.reference || true
regression_require_output "concurrent pseudo trace" postProcessing/pseudoECG.dat '^#.*time'
regression_require_output "independent pseudo references" postProcessing/manufacturedPseudoECG_Pseudo.dat '^#.*time'
if ! python3 "$(dirname "$regressionLib")/checkDualECG.py" .; then
    regression_fail "concurrent ECG functional or quadrature check failed"
fi
regression_finish
