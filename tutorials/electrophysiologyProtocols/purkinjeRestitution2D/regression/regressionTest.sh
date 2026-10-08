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

# Runs each variant in parallel and compares postProcessing/purkinjeNetwork.dat
# against regression/<variant>.reference. The monodomain variant also runs
# the graph-only runPurkinjeGraph utility and checks it against
# regression/monodomain.graphUtility.reference.
#
# Regression configuration (tutorials/README.md): retrograde stops at 0.6 s,
# once its second beat has crossed the network, instead of 0.9 s. antegrade,
# 2.45 s of which the tissue is idle until the 1.2 s escape beat, runs last,
# on a 75 x 75 slab instead of 150 x 150: what it checks are network
# activation times, which do not change with the slab resolution. The
# retrograde and monodomain variants need the fine slab, since the junctions
# exchange current with the tissue.

VARIANTS=(retrograde monodomain antegrade)

regression_init "purkinjeRestitution2D regression test" \
    "regression/${VARIANTS[0]}.reference" "$@"
regression_require_run_mode
REGRESSION_TIME_WINDOW=1e-6
regression_set system/controlDict.retrograde endTime 0.6

# macOS strips DYLD_LIBRARY_PATH from child processes; runPurkinjeGraph is
# called directly rather than through RunFunctions.
if [[ "$(uname -s)" == "Darwin" && -n "${WM_PROJECT_DIR:-}" && -n "${WM_OPTIONS:-}" ]]; then
    openfoamLibDir="${WM_PROJECT_DIR}/platforms/${WM_OPTIONS}/lib"
    if [[ -d "${openfoamLibDir}" ]]; then
        export DYLD_LIBRARY_PATH="${openfoamLibDir}:${DYLD_LIBRARY_PATH:-}"
    fi
fi

for variant in "${VARIANTS[@]}"; do
    echo "------------------------------------------------------------"
    echo "Variant: ${variant}"
    echo "------------------------------------------------------------"

    if [[ "${variant}" == antegrade ]]; then
        regression_edit system/blockMeshDict \
            's/^([[:space:]]*hex[[:space:]]*\([^)]*\)[[:space:]]*)\([^)]*\)/\1(75 75 1)/'
    fi

    if ! regression_run "${variant}" parallel; then
        regression_fail "Allrun ${variant} parallel did not complete"
        continue
    fi
    regression_compare "regression/${variant}.reference" || true

    if [[ "${variant}" == monodomain ]]; then
        rm -rf postProcessing
        if runPurkinjeGraph -case . -conductionDomain purkinjeNetwork > log.runPurkinjeGraph 2>&1; then
            regression_compare regression/monodomain.graphUtility.reference || true
        else
            regression_fail "runPurkinjeGraph exited non-zero"
            regression_log_tail runPurkinjeGraph log.runPurkinjeGraph
        fi
    fi
    echo
done

regression_finish
