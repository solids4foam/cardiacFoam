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

# Single-cell traces: samples Vm, cai and the active tension at two times of
# the paced TWorld endocardial cell.

regression_init "Single-cell regression test" regression/singleCell.reference "$@"
REGRESSION_TIME_WINDOW=1e-4

regression_run_or_fail
regression_compare || true
regression_finish
