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

# The spherical-shell solution on the hex cubed sphere (N = 12, configuration
# ATI, exact current on the outer box), run in parallel on four ranks with the
# manual decomposition so that interface faces are both internal to a rank and
# on processor patches, and two ranks hold no heart cells. The verifier writes
# postProcessing/3D_hex_<N>_cells.dat; its summary values and error norms are
# compared against the reference.
regression_init "Bath-bidomain spherical-shell transmission manufactured-solution regression test" \
    regression/sphericalShell.reference "$@"

# regression/skips-on-openfoam: the regression is suppressed on OpenFOAM v2312,
# where the cellZone myocardium sub-mesh fails its fvSchemes lookups.
if [[ "${WM_PROJECT_VERSION:-unknown}" == *2312* ]]; then
    echo "SKIP: spherical-shell regression is suppressed on OpenFOAM v2312 (regression/skips-on-openfoam)."
    exit "${REGRESSION_SKIP_CODE}"
fi

regression_run_or_fail parallel -N 12 -endTime 0.02

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Bath-bidomain manufactured solution error summary' 'postProcessing/*_cells.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "bath-bidomain transmission manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

regression_compare || true
regression_finish
