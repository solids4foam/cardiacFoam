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

VARIANTS=(antegrade retrograde monodomain)

regression_init "purkinjeRestitution2D regression test" \
    "regression/${VARIANTS[0]}.reference" "$@"
regression_require_run_mode
REGRESSION_TIME_WINDOW=1e-6

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
