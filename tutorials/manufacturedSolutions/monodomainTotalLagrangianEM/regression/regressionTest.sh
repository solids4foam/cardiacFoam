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
# Compares the final-time error norms that manufacturedElectromechanicsVerifier
# writes for Vm, D, lambda and Ta against the reference.
#
# The tracked case runs 40^3 cells at the N=80 time step. The regression runs
# 20^3 cells at the matching N=20 time step from setup/driver_config.json.

regression_init "Electromechanics manufactured-solution regression test" \
    regression/monodomainTotalLagrangianEM.reference "$@"
regression_require_solids4foam

regression_edit system/blockMeshDict \
    's/^([[:space:]]*hex[[:space:]]*\([^)]*\)[[:space:]]*)\([^)]*\)/\1(20 20 20)/'
regression_set system/controlDict deltaT 0.00224215

regression_run_or_fail

REGRESSION_SUMMARY_FILE="postProcessing/manufacturedElectromechanicsSummary.dat"
regression_require_output "manufactured electromechanics error summary" \
    "${REGRESSION_SUMMARY_FILE}" 'Field     L1-error'

regression_compare || true
regression_finish
