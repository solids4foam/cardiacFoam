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
# The tracked case runs 40^3 cells at the N=80 time step, too slow for a
# regression. This script runs 20^3 cells at the matching N=20 time step from
# setup/driver_config.json: blockMeshDict and controlDict are rewritten for
# this invocation only and restored on exit.

REGRESSION_CELLS="20 20 20"
REGRESSION_DELTA_T="0.00224215"
BLOCKMESH_DICT="system/blockMeshDict"
CONTROL_DICT="system/controlDict"
BACKUP_SUFFIX=".regressionTest.bak"

regression_init "Electromechanics manufactured-solution regression test" \
    regression/monodomainTotalLagrangianEM.reference "$@"
regression_require_solids4foam

if (( ! REGRESSION_CHECK_ONLY )); then
    cp "${BLOCKMESH_DICT}" "${BLOCKMESH_DICT}${BACKUP_SUFFIX}"
    cp "${CONTROL_DICT}" "${CONTROL_DICT}${BACKUP_SUFFIX}"
    regression_add_exit_hook \
        'mv -f "${BLOCKMESH_DICT}${BACKUP_SUFFIX}" "${BLOCKMESH_DICT}"; mv -f "${CONTROL_DICT}${BACKUP_SUFFIX}" "${CONTROL_DICT}"'

    sed -E "s/^([[:space:]]*hex[[:space:]]*\([^)]*\)[[:space:]]*)\([^)]*\)/\1(${REGRESSION_CELLS})/" \
        "${BLOCKMESH_DICT}${BACKUP_SUFFIX}" > "${BLOCKMESH_DICT}"
    sed -E "s/^deltaT[[:space:]]+[^;]+;/deltaT          ${REGRESSION_DELTA_T};/" \
        "${CONTROL_DICT}${BACKUP_SUFFIX}" > "${CONTROL_DICT}"

    if ! grep -q "(${REGRESSION_CELLS})" "${BLOCKMESH_DICT}" \
        || ! grep -qE "^deltaT[[:space:]]+${REGRESSION_DELTA_T};" "${CONTROL_DICT}"; then
        regression_fail "could not set the regression mesh size or time step"
        regression_finish
    fi
fi

regression_run_or_fail

REGRESSION_SUMMARY_FILE="postProcessing/manufacturedElectromechanicsSummary.dat"
regression_require_output "manufactured electromechanics error summary" \
    "${REGRESSION_SUMMARY_FILE}" 'Field     L1-error'

regression_compare || true
regression_finish
