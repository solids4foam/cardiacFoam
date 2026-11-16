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

# Two runs of the tutorial (serial, seconds each):
# - as configured, the paced TWorld endocardial cell: Vm, cai and the active
#   tension at two times (regression/singleCell.reference);
# - the same protocol with ionicModel Stewart, tissue myocyte, the only
#   regression of the Stewart model: Vm on both beats
#   (regression/Stewart.reference). The ionic model and tissue are set for
#   this run only, as setup/sweep_ionic_model_tissue.json sets them.

regression_init "Single-cell regression test" regression/singleCell.reference "$@"
regression_require_run_mode
REGRESSION_TIME_WINDOW=1e-4

echo "Model: TWorld endocardialCells"
regression_run_or_fail
regression_compare regression/singleCell.reference || true

echo "Model: Stewart myocyte"
regression_set constant/electroProperties ionicModel Stewart
regression_set constant/electroProperties tissue myocyte
regression_run_or_fail
regression_compare regression/Stewart.reference || true

regression_finish
