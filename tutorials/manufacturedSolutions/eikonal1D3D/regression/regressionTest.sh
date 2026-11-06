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

# The case defaults for each seed profile (one terminal through the graph and the coupler, serial): the
# errors against the exact radial solution, postProcessing/radial_errors.dat. The
# mesh comes from gmsh, so the case skips where it is missing (regression/needs-tools).
# The errors depend on the gmsh mesh; the tolerances in the reference cover the spread
# over mesh seeds.
if ! command -v gmsh > /dev/null; then
    echo "SKIP: this case needs gmsh"
    exit "${REGRESSION_SKIP_CODE}"
fi

regression_init "Eikonal 1D-3D exact radial regression test" \
    regression/eikonal1D3D.reference "$@"

for seedProfile in front flat curvedFront; do
    refFile="regression/eikonal1D3D.reference"
    if [[ "${seedProfile}" != front ]]; then
        refFile="regression/eikonal1D3D.${seedProfile}.reference"
    fi

    echo "Seed profile: ${seedProfile}"
    # curvedFront needs its critical radius resolved: 2/Pe_L against h/L.
    runArgs=(-seedProfile "${seedProfile}")
    if [[ "${seedProfile}" == curvedFront ]]; then
        runArgs+=(-PeL 21 -N 20)
    fi

    # An existing run holds one profile, which is the only one to check.
    if (( REGRESSION_CHECK_ONLY )) \
        && ! grep -qx "# seedProfile ${seedProfile}" postProcessing/radial_errors.dat; then
        continue
    fi

    regression_run_or_fail "${runArgs[@]}"

    # A run that reaches the iteration cap still prints End
    regression_check_solver_logs -e 'PIMPLE: converged' || { regression_write_report; exit 1; }

    REGRESSION_SUMMARY_FILE="postProcessing/radial_errors.dat"
    regression_require_output "radial eikonal error summary" \
        "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

    regression_compare "${refFile}" || true
done
regression_finish
