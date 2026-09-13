#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Idealized-heart electromechanics regression test
# ============================================================
#
# FULL SOLIDS4FOAM BUILD ONLY. Under
# CARDIAC_REGRESSION_BUILD_MODE=lightweight this exits 77 before touching the
# case and Alltest-regression reports an expected skip; in with-solids4foam
# mode a skip is a failure.
#
# The case runs to its own endTime of 0.02 s and probes two cells: one inside
# the apical stimulus region, which activates early, and one mid-wall, which
# does not activate within 20 ms. The reference pins the activation times, the
# displacement D at both, and the apical active tension.

REF_FILE="regression/electroMechHeart.reference"
ALLRUN_LOGFILE="log.Allrun"
SKIP_CODE=77
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

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
        "${REPO_ROOT}/modules/solids4foam" \
        "${HOME}/solids4foam" \
        "${WM_PROJECT_USER_DIR:-}/solids4foam"
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
    awk -v a="$1" -v e="$2" 'BEGIN { d = a - e; if (d < 0) d = -d; print d; }'
}

checkWithinTolerance()
{
    local label="$1" actual="$2" expected="$3" tolerance="$4" diffAbs

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

# Vector probes write "(x y z)" per location, so strip the parentheses before
# indexing: time is column 1, and a vector probe then occupies three columns.
extractProbeValue()
{
    local dataFile="$1" time="$2" column="$3"

    awk -v target="${time}" -v col="${column}" '
        BEGIN { bestDiff = 1e99; found = 0; actual = 0.0; }
        /^#/ { next }
        {
            gsub(/[()]/, "");
            if (NF < col) next;
            d = $1 - target;
            if (d < 0) d = -d;
            if (d < bestDiff) { bestDiff = d; actual = $col; found = 1; }
        }
        END {
            if (found && bestDiff <= 2.5e-3) { print actual; exit 0; }
            exit 1;
        }
    ' "${dataFile}"
}

echo "============================================================"
echo "Idealized-heart electromechanics regression test"
echo "============================================================"
echo

if [[ "${CARDIAC_REGRESSION_BUILD_MODE:-}" == "lightweight" ]]; then
    echo "SKIP: electromechanics regression requires a full solids4foam build, but lightweight mode was specified."
    exit "${SKIP_CODE}"
fi

if ! fullSolids4FoamAvailable; then
    echo "SKIP: electromechanics regression requires a full solids4foam build."
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

failures=0
checks=0

while IFS=' ' read -r fileName time column expected tolerance; do
    if [[ -z "${fileName}" || "${fileName}" == \#* ]]; then
        continue
    fi

    dataFile="postProcessing/${fileName}"
    checks=$((checks + 1))

    if [[ ! -f "${dataFile}" ]]; then
        echo "FAIL: missing output file ${dataFile}"
        failures=$((failures + 1))
        continue
    fi

    actual="$(extractProbeValue "${dataFile}" "${time}" "${column}")" || true

    if [[ -z "${actual}" ]]; then
        echo "FAIL: ${dataFile} col=${column} at t=${time} not found"
        failures=$((failures + 1))
        continue
    fi

    checkWithinTolerance "${dataFile} col=${column} t=${time}" \
        "${actual}" "${expected}" "${tolerance}" \
        || failures=$((failures + 1))
done < "${REF_FILE}"

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
