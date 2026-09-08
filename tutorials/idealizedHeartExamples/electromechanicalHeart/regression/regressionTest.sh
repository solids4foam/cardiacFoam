#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Electromechanical Niederer slab regression test
# ============================================================
#
# No reference file is checked in yet for this case - this script is a
# scaffold, not a wired-up regression test. It will report a controlled
# "FAIL: reference file not found" if run as-is.

REF_FILE="regression/electroMechHeterogeneity.reference"
ALLRUN_LOGFILE="log.Allrun"
PROBE_FILE="postProcessing/Taprobes/solid/0/Ta"
DISTINCT_EPS=1.0
SKIP_CODE=77
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

fullSolids4FoamAvailable()
{
    local solidsDir="${SOLIDS4FOAM_INST_DIR:-}"
    local candidate
    local builtMarker="src/solids4FoamModels/lnInclude/solidModel.H"

    if [[ -n "${solidsDir}" ]] \
        && [[ -f "${solidsDir}/src/solids4FoamModels/solidModels/solidModel/solidModel.H" ]]; then
        return 0
    fi

    for candidate in \
        "${HOME}/solids4foam" \
        "${WM_PROJECT_USER_DIR:-}/solids4foam" \
        "${REPO_ROOT}/modules/solids4foam"
    do
        if [[ -n "${candidate}" ]] && [[ -f "${candidate}/${builtMarker}" ]]; then
            export SOLIDS4FOAM_INST_DIR="${candidate}"
            return 0
        fi
    done

    return 1
}

electroMechanicalLibCompiled()
{
    local libDir="${FOAM_USER_LIBBIN:-}"
    local lib

    [[ -n "${libDir}" ]] || return 1

    for lib in "${libDir}"/libelectroMechanicalModels.*; do
        if [[ -e "${lib}" ]]; then
            return 0
        fi
    done

    return 1
}

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

absDiff()
{
    awk -v a="$1" -v e="$2" '
        BEGIN {
            d = a - e;
            if (d < 0) d = -d;
            print d;
        }
    '
}

checkWithinTolerance()
{
    local label="$1"
    local actual="$2"
    local expected="$3"
    local tolerance="$4"
    local diffAbs

    diffAbs="$(absDiff "${actual}" "${expected}")"

    if awk -v d="${diffAbs}" -v t="${tolerance}" 'BEGIN {exit !(d <= t)}'; then
        printf "PASS: %s actual=%.8g expected=%.8g difference=%.3g tolerance=%.3g\n" \
            "${label}" "${actual}" "${expected}" "${diffAbs}" "${tolerance}"
        return 0
    fi

    printf "FAIL: %s actual=%.8g expected=%.8g difference=%.3g tolerance=%.3g\n" \
        "${label}" "${actual}" "${expected}" "${diffAbs}" "${tolerance}"
    return 1
}

extractProbeValue()
{
    local dataFile="$1"
    local time="$2"
    local column="$3"

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
            if (found && bestDiff <= 2.5e-3) {
                print actual;
                exit 0;
            }
            exit 1;
        }
    ' "${dataFile}"
}

echo "============================================================"
echo "Electromechanical Niederer slab regression test"
echo "============================================================"
echo

if [[ "${CARDIAC_REGRESSION_BUILD_MODE:-}" == "lightweight" ]]; then
    echo "SKIP: electromechanical regression requires a full solids4foam build, but lightweight mode was specified."
    exit "${SKIP_CODE}"
fi

if ! fullSolids4FoamAvailable; then
    echo "SKIP: electromechanical regression requires a full solids4foam build."
    echo "      SOLIDS4FOAM_INST_DIR does not point to a compiled full solids4foam tree."
    exit "${SKIP_CODE}"
fi

if ! electroMechanicalLibCompiled; then
    echo "FAIL: full solids4foam is available, but libelectroMechanicalModels is not compiled."
    echo "      Rebuild cardiacFoam in full mode before running this regression."
    exit 1
fi

./Allclean > /dev/null 2>&1 || true

if ! ./Allrun parallel > "${ALLRUN_LOGFILE}" 2>&1; then
    echo "FAIL: Allrun exited non-zero. Surfacing logs:"
    dumpLogTail "Allrun" "${ALLRUN_LOGFILE}"
    for stage in blockMesh setExprFields decomposePar cardiacFoam; do
        dumpLogTail "${stage}" "log.${stage}"
    done
    exit 1
fi

if [[ ! -f "${REF_FILE}" ]]; then
    echo "FAIL: reference file not found: ${REF_FILE}"
    exit 1
fi

if [[ ! -f "${PROBE_FILE}" ]]; then
    echo "FAIL: probe output not found: ${PROBE_FILE}"
    exit 1
fi

failures=0
checks=0

while IFS=' ' read -r fileName time column expected tolerance; do
    if [[ -z "${fileName}" || "${fileName}" == \#* ]]; then
        continue
    fi

    dataFile="postProcessing/${fileName}"
    if [[ ! -f "${dataFile}" ]]; then
        echo "FAIL: missing output file ${dataFile}"
        failures=$((failures + 1))
        checks=$((checks + 1))
        continue
    fi

    actual="$(extractProbeValue "${dataFile}" "${time}" "${column}")" || true
    checks=$((checks + 1))

    if [[ -z "${actual}" ]]; then
        echo "FAIL: ${dataFile} col=${column} at t=${time} not found"
        failures=$((failures + 1))
        continue
    fi

    checkWithinTolerance "${dataFile} col=${column} t=${time}" \
        "${actual}" "${expected}" "${tolerance}" \
        || failures=$((failures + 1))
done < "${REF_FILE}"

endo="$(extractProbeValue "${PROBE_FILE}" "0.02" 2)" || {
    echo "FAIL: could not extract endocardial Ta value for heterogeneity check"
    exit 1
}
mid="$(extractProbeValue "${PROBE_FILE}" "0.02" 3)" || {
    echo "FAIL: could not extract mid-wall Ta value for heterogeneity check"
    exit 1
}
epi="$(extractProbeValue "${PROBE_FILE}" "0.02" 4)" || {
    echo "FAIL: could not extract epicardial Ta value for heterogeneity check"
    exit 1
}

if awk -v a="${endo}" -v b="${mid}" -v c="${epi}" -v eps="${DISTINCT_EPS}" '
    BEGIN {
        ab = a - b; if (ab < 0) ab = -ab;
        ac = a - c; if (ac < 0) ac = -ac;
        bc = b - c; if (bc < 0) bc = -bc;
        exit !((ab > eps) && (ac > eps) && (bc > eps));
    }
' ; then
    printf "PASS: heterogeneity check endo=%.8g mid=%.8g epi=%.8g\n" \
        "${endo}" "${mid}" "${epi}"
    checks=$((checks + 1))
else
    printf "FAIL: heterogeneity collapsed endo=%.8g mid=%.8g epi=%.8g\n" \
        "${endo}" "${mid}" "${epi}"
    failures=$((failures + 1))
    checks=$((checks + 1))
fi

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
fi

echo "============================================================"
echo "Regression test FAILED (${failures}/${checks} checks)"
echo "============================================================"
exit 1
