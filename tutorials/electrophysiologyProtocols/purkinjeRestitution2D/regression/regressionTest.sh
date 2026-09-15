#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# purkinjeRestitution2D regression test
# ============================================================
#
# Runs each variant in parallel and compares postProcessing/purkinjeNetwork.dat
# against regression/<variant>.reference (rows: file time column expected
# tolerance). The monodomain variant also runs the graph-only runPurkinjeGraph
# utility and checks it against regression/monodomain.graphUtility.reference.

VARIANTS=(antegrade retrograde monodomain)
ALLRUN_LOGFILE="log.Allrun"
GRAPH_LOGFILE="log.runPurkinjeGraph"

# macOS strips DYLD_LIBRARY_PATH from child processes; runPurkinjeGraph is
# called directly rather than through RunFunctions.
if [[ "$(uname -s)" == "Darwin" && -n "${WM_PROJECT_DIR:-}" && -n "${WM_OPTIONS:-}" ]]; then
    openfoamLibDir="${WM_PROJECT_DIR}/platforms/${WM_OPTIONS}/lib"
    if [[ -d "${openfoamLibDir}" ]]; then
        export DYLD_LIBRARY_PATH="${openfoamLibDir}:${DYLD_LIBRARY_PATH:-}"
    fi
fi

echo "============================================================"
echo "purkinjeRestitution2D regression test"
echo "============================================================"
echo

dumpLogTail()
{
    local label="$1"
    local logFile="$2"
    local maxLines="${3:-80}"

    if [[ -s "${logFile}" ]]; then
        echo "----- last ${maxLines} lines of ${label} (${logFile}) -----"
        tail -n "${maxLines}" "${logFile}"
        echo "----- end of ${label} -----"
    else
        echo "(no log file at ${logFile})"
    fi
}

# Compare one reference file's rows against postProcessing/.
# Echoes PASS/FAIL per row; returns the number of failures.
checkReference()
{
    local refFile="$1"
    local failures=0

    while IFS=' ' read -r fileName time column expected tolerance; do
        if [[ -z "${fileName}" || "${fileName}" == \#* ]]; then
            continue
        fi

        local dataFile="postProcessing/${fileName}"
        if [[ ! -f "${dataFile}" ]]; then
            echo "FAIL: missing output file ${dataFile}"
            failures=$((failures + 1))
            continue
        fi

        local actual
        actual="$(
            awk -v target="${time}" -v col="${column}" '
                BEGIN { bestDiff = 1e99; found = 0; actual = 0.0; }
                $1 !~ /^#/ && NF >= col {
                    d = $1 - target;
                    if (d < 0) d = -d;
                    if (d < bestDiff) {
                        bestDiff = d;
                        actual = $col;
                        found = 1;
                    }
                }
                END {
                    if (found && bestDiff <= 1e-6) {
                        print actual;
                        exit 0;
                    }
                    exit 1;
                }
            ' "${dataFile}"
        )" || true

        if [[ -z "${actual}" ]]; then
            echo "FAIL: ${dataFile} col=${column} at t=${time} not found"
            failures=$((failures + 1))
            continue
        fi

        local diffAbs
        diffAbs="$(
            awk -v a="${actual}" -v e="${expected}" \
                'BEGIN { d = a - e; if (d < 0) d = -d; print d; }'
        )"

        if awk -v d="${diffAbs}" -v t="${tolerance}" 'BEGIN {exit !(d <= t)}'; then
            printf "PASS: %s col=%s t=%s value=%.9g (difference = %.3g)\n" \
                "${dataFile}" "${column}" "${time}" "${actual}" "${diffAbs}"
        else
            printf "FAIL: %s col=%s t=%s value=%.9g expected=%.9g (difference = %.3g)\n" \
                "${dataFile}" "${column}" "${time}" "${actual}" "${expected}" "${diffAbs}"
            failures=$((failures + 1))
        fi
    done < "${refFile}"

    return "${failures}"
}

failures=0
failedVariants=()

for variant in "${VARIANTS[@]}"; do
    refFile="regression/${variant}.reference"

    echo "------------------------------------------------------------"
    echo "Variant: ${variant}"
    echo "------------------------------------------------------------"

    ./Allclean > /dev/null 2>&1 || true

    if ! ./Allrun "${variant}" parallel > "${ALLRUN_LOGFILE}" 2>&1 \
        || ! grep -q "^End" log.cardiacFoam; then
        echo "FAIL: Allrun ${variant} parallel did not complete. Surfacing logs:"
        dumpLogTail "Allrun" "${ALLRUN_LOGFILE}"
        dumpLogTail "cardiacFoam" "log.cardiacFoam"
        failures=$((failures + 1))
        failedVariants+=("${variant}")
        echo
        continue
    fi

    variantFailures=0
    checkReference "${refFile}" || variantFailures=$?

    if [[ "${variant}" == "monodomain" ]]; then
        rm -rf postProcessing
        if runPurkinjeGraph -case . -conductionDomain purkinjeNetwork \
            > "${GRAPH_LOGFILE}" 2>&1; then
            graphFailures=0
            checkReference "regression/monodomain.graphUtility.reference" \
                || graphFailures=$?
            variantFailures=$((variantFailures + graphFailures))
        else
            echo "FAIL: runPurkinjeGraph exited non-zero."
            dumpLogTail "runPurkinjeGraph" "${GRAPH_LOGFILE}"
            variantFailures=$((variantFailures + 1))
        fi
    fi

    if (( variantFailures > 0 )); then
        failures=$((failures + variantFailures))
        failedVariants+=("${variant}")
    fi
    echo
done

if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED (${#VARIANTS[@]} variants)"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures} check(s) in: ${failedVariants[*]})"
    echo "============================================================"
    exit 1
fi
