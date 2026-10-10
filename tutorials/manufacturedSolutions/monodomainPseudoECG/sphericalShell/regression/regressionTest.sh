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

# The monodomain sealed-shell solution on the hex cubed sphere (N = 12), run in
# parallel on four ranks (scotch); the verifier writes
# postProcessing/3D_hex_12_cells.dat, whose summary values and error norms are
# compared against the reference.
regression_init "Monodomain sealed-shell manufactured-solution regression test" \
    regression/sphericalShell.reference "$@"

regression_run_or_fail parallel -N 12 -endTime 0.02

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Monodomain sealed-shell manufactured solution error summary' 'postProcessing/*_cells.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "monodomain sealed-shell manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

regression_compare || true
regression_finish
