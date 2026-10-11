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

# Coupled monodomain + total-Lagrangian solid manufactured solution. Needs the
# full solids4foam build: exits 77 under CARDIAC_REGRESSION_BUILD_MODE=lightweight.
#
# Runs the case at N 20 with -dt 0.00224215, the step of the reference values,
# and compares the final-time error norms that manufacturedElectromechanicsVerifier
# writes for Vm, D, lambda and Ta against the reference.

regression_init "Electromechanics manufactured-solution regression test" \
    regression/monodomainTotalLagrangianEM.reference "$@"
regression_require_solids4foam

regression_run_or_fail -dt 0.00224215

REGRESSION_SUMMARY_FILE="postProcessing/manufacturedElectromechanicsSummary.dat"
regression_require_output "manufactured electromechanics error summary" \
    "${REGRESSION_SUMMARY_FILE}" 'Field     L1-error'

regression_compare || true
regression_finish
