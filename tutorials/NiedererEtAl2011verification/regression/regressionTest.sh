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

# Niederer et al. (2011) slab benchmark: compares the activation-time probes
# written by the Niedererpoints function object against the reference.

regression_init "Niederer et al. (2011) activation regression test" \
    regression/NiedererEtAl2011.reference "$@"

regression_run_or_fail parallel
regression_compare || true
regression_finish
