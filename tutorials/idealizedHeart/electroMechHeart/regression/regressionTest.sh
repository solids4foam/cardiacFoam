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

# Idealized-heart electromechanics. Needs the full solids4foam build: exits 77
# under CARDIAC_REGRESSION_BUILD_MODE=lightweight.
#
# The case runs to its own endTime of 0.02 s and probes two cells: one inside
# the apical stimulus region, which activates early, and one mid-wall, which
# does not activate within 20 ms. The reference pins the activation times,
# the displacement D at both, and the apical active tension.

regression_init "Idealized-heart electromechanics regression test" \
    regression/electroMechHeart.reference "$@"
regression_require_solids4foam

regression_run_or_fail parallel
regression_compare || true
regression_finish
