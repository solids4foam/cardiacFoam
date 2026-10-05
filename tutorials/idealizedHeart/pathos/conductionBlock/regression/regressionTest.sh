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
# The case's own controlDict runs to 0.7 s. endTime is set to 0.035 s for
# this script's invocation only and restored on exit.

VARIANTS=(lbbb rbbb)
CONTROL_DICT="system/controlDict"
CONTROL_DICT_BACKUP="system/controlDict.regressionTest.bak"

regression_init "Idealized heart conduction-block regression test" \
    "regression/${VARIANTS[0]}.reference" "$@"
regression_require_run_mode

cp "${CONTROL_DICT}" "${CONTROL_DICT_BACKUP}"
regression_add_exit_hook 'mv -f "${CONTROL_DICT_BACKUP}" "${CONTROL_DICT}"'
sed -E 's/^endTime[[:space:]]+[^;]+;/endTime    0.035;/' \
    "${CONTROL_DICT_BACKUP}" > "${CONTROL_DICT}"

for variant in "${VARIANTS[@]}"; do
    echo "------------------------------------------------------------"
    echo "Variant: ${variant}"
    echo "------------------------------------------------------------"

    if ! regression_run "${variant}"; then
        regression_fail "Allrun ${variant} did not complete"
        continue
    fi
    regression_compare "regression/${variant}.reference" || true
    echo
done

regression_finish
