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

# Bidomain Niederer et al. (2011) slab: runs both conductivity sets on the
# 0.5 mm mesh (-dx 5e-4, 40 x 6 x 14 cells) to 15 ms at the default time step,
# equal in parallel and published in serial, and compares the activation
# times of the Niedererpoints probes and the Niedererlines diagonal against
# regression/<set>.reference.

regression_init "Niederer et al. (2011) bidomain activation regression test" \
    regression/equal.reference "$@"
regression_require_run_mode

for conductivity in equal published; do
    echo "------------------------------------------------------------"
    echo "Conductivity: ${conductivity}"
    echo "------------------------------------------------------------"

    runArgs=(-dx 5e-4 -endTime 0.015 -conductivity "${conductivity}")
    [[ "${conductivity}" == equal ]] && runArgs=(parallel "${runArgs[@]}")

    if ! regression_run "${runArgs[@]}"; then
        regression_fail "Allrun ${runArgs[*]} did not complete"
        continue
    fi
    regression_compare "regression/${conductivity}.reference" || true
    echo
done

regression_finish
