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

# The 1D mesh (system/blockMeshDict.1D, 80 cells per block) matches the
# dimension "1D" entry in constant/electroProperties. The verifier writes
# postProcessing/<DIM>_<N>_cells.dat; its summary values and
# error norms are compared against the reference.
regression_init "Bath-bidomain manufactured-solution regression test" \
    regression/bathBidomainManufactured.reference "$@"

# OpenFOAM v2312's dictionary lookup of 'laplacian(conductivityIntracellular,Vm)'
# against system/fvSchemes fails for this case in both build modes, while
# v2412 and v2512 pass. The regression is suppressed on v2312 only.
if [[ "${WM_PROJECT_VERSION:-unknown}" == *2312* ]]; then
    echo "SKIP: bathBidomain regression is suppressed on OpenFOAM v2312 (regression/skips-on-openfoam)."
    exit "${REGRESSION_SKIP_CODE}"
fi

regression_run_or_fail

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Bath-bidomain manufactured solution error summary' 'postProcessing/*_cells.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "bath-bidomain manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

regression_compare || true
regression_finish
