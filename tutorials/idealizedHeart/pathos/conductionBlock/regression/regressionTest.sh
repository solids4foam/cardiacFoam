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

# Confirms that severing a bundle branch blocks fast conduction. Two
# activation-time probes: the LV Purkinje-myocardial-junction site (node 211,
# the point electroHeart's regression requires to be activated by 28.9 ms)
# and an RV site. In lbbb the LV probe must still read -1 at t = 0.035 s
# while the RV probe has activated; in rbbb the roles swap. Myocardial spread
# reaches the blocked LV point at 43.2 ms, so the cutoff separates block from
# health with margin on both sides.
#
# The case's own controlDict runs to 0.7 s; the regression stops at 0.035 s
# and runs in parallel, as the other whole-heart regressions do.

VARIANTS=(lbbb rbbb)

regression_init "Idealized heart conduction-block regression test" \
    "regression/${VARIANTS[0]}.reference" "$@"
regression_require_run_mode

regression_set system/controlDict endTime 0.035

for variant in "${VARIANTS[@]}"; do
    echo "------------------------------------------------------------"
    echo "Variant: ${variant}"
    echo "------------------------------------------------------------"

    if ! regression_run "${variant}" parallel; then
        regression_fail "Allrun ${variant} parallel did not complete"
        continue
    fi
    regression_compare "regression/${variant}.reference" || true
    echo
done

regression_finish
