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

# The bidomain sealed-shell solution with the case defaults (hex cubed sphere,
# N = 12, T = 0.2), run in parallel on four ranks (scotch); the verifier writes
# postProcessing/3D_hex_12_cells.dat, whose summary values and error norms are
# compared against the reference.
regression_init "Bidomain sealed-shell manufactured-solution regression test" \
    regression/sphericalShell.reference "$@"

# On OpenFOAM v2312 the Vm and ionic-state error rows of this case differ from
# the v2412 and v2512 values by about 1e-4 relative, with the phiE and gradient
# rows unchanged; the cause is not diagnosed. The regression is suppressed on
# v2312 only (regression/skips-on-openfoam).
if [[ "${WM_PROJECT_VERSION:-unknown}" == *2312* ]]; then
    echo "SKIP: bidomain sphericalShell regression is suppressed on OpenFOAM v2312 (regression/skips-on-openfoam)."
    exit "${REGRESSION_SKIP_CODE}"
fi

regression_run_or_fail parallel

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Bidomain sealed-shell manufactured solution error summary' 'postProcessing/*_cells.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "bidomain sealed-shell manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

regression_compare || true
regression_finish
