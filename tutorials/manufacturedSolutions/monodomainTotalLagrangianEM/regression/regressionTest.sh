#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Electromechanics manufactured-solution regression test
# ============================================================
#
# FULL SOLIDS4FOAM BUILD ONLY. This case needs cardiacFoam built with
# solids4foam (libelectroMechanicalModels). Under
# CARDIAC_REGRESSION_BUILD_MODE=lightweight it exits 77 before touching the
# case, and Alltest-regression reports an expected skip. In with-solids4foam
# mode a skip is a failure: no solids4foam tree found exits 77, which
# Alltest-regression counts as an unexpected skip, and a missing
# libelectroMechanicalModels fails outright.
#
# Runs the coupled monodomain + total-Lagrangian solid manufactured solution
# and compares the final-time error norms that
# manufacturedElectromechanicsVerifier writes for Vm, D, lambda and Ta to
# the reference.
#
# The tracked case runs 40^3 cells at the N=80 time step, far too slow for a
# regression. This script runs 20^3 cells at the matching N=20 time step from
# setup/driver_config.json instead: still a full 3D coupled solve, in about
# 45 steps. blockMeshDict and controlDict are rewritten for this invocation
# only and restored on exit - the tracked files are never left changed.
# Allrun itself is reused unmodified.

REF_FILE="regression/monodomainTotalLagrangianEM.reference"
ALLRUN_LOGFILE="log.Allrun"
SUMMARY_FILE="postProcessing/manufacturedElectromechanicsSummary.dat"
REGRESSION_CELLS="20 20 20"
REGRESSION_DELTA_T="0.00224215"
BLOCKMESH_DICT="system/blockMeshDict"
CONTROL_DICT="system/controlDict"
BACKUP_SUFFIX=".regressionTest.bak"
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

extractFinalTime()
{
    sed -nE 's/.*error summary \(t = ([^)]+)\):.*/\1/p' "$1" | head -n 1
}

extractErrorMetric()
{
    local summaryFile="$1"
    local field="$2"
    local metric="$3"
    local column

    case "${metric}" in
        L1) column=2 ;;
        L2) column=3 ;;
        Linf) column=4 ;;
        *) return 1 ;;
    esac

    awk -v field="${field}" -v column="${column}" '$1 == field {print $column; exit}' "${summaryFile}"
}

checkReferenceValues()
{
    local summaryFile="$1"
    local failures=0
    local checks=0
    local kind key metric expected tolerance actual

    if [[ ! -f "${REF_FILE}" ]]; then
        echo "FAIL: reference file not found: ${REF_FILE}"
        return 1
    fi

    while IFS=' ' read -r kind key metric expected tolerance; do
        if [[ -z "${kind}" || "${kind}" == \#* ]]; then
            continue
        fi

        case "${kind}" in
            summary)
                if [[ "${key}" == "finalTime" ]]; then
                    actual="$(extractFinalTime "${summaryFile}")"
                else
                    actual=""
                fi
                ;;
            error)
                actual="$(extractErrorMetric "${summaryFile}" "${key}" "${metric}")"
                ;;
            *)
                echo "FAIL: unknown reference kind '${kind}' in ${REF_FILE}"
                failures=$((failures + 1))
                continue
                ;;
        esac

        checks=$((checks + 1))

        if [[ -z "${actual}" ]]; then
            echo "FAIL: could not extract ${kind} ${key} ${metric}"
            failures=$((failures + 1))
            continue
        fi

        checkWithinTolerance "${kind} ${key} ${metric}" "${actual}" "${expected}" "${tolerance}" \
            || failures=$((failures + 1))
    done < "${REF_FILE}"

    echo "Electromechanics manufactured reference comparison: ${checks} checks, ${failures} failures"
    (( failures == 0 ))
}

echo "============================================================"
echo "Electromechanics manufactured-solution regression test"
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

cp "${BLOCKMESH_DICT}" "${BLOCKMESH_DICT}${BACKUP_SUFFIX}"
cp "${CONTROL_DICT}" "${CONTROL_DICT}${BACKUP_SUFFIX}"
restoreDicts()
{
    mv -f "${BLOCKMESH_DICT}${BACKUP_SUFFIX}" "${BLOCKMESH_DICT}"
    mv -f "${CONTROL_DICT}${BACKUP_SUFFIX}" "${CONTROL_DICT}"
}
trap restoreDicts EXIT

sed -E "s/^([[:space:]]*hex[[:space:]]*\([^)]*\)[[:space:]]*)\([^)]*\)/\1(${REGRESSION_CELLS})/" \
    "${BLOCKMESH_DICT}${BACKUP_SUFFIX}" > "${BLOCKMESH_DICT}"
sed -E "s/^deltaT[[:space:]]+[^;]+;/deltaT          ${REGRESSION_DELTA_T};/" \
    "${CONTROL_DICT}${BACKUP_SUFFIX}" > "${CONTROL_DICT}"

if ! grep -q "(${REGRESSION_CELLS})" "${BLOCKMESH_DICT}" \
    || ! grep -qE "^deltaT[[:space:]]+${REGRESSION_DELTA_T};" "${CONTROL_DICT}"; then
    echo "FAIL: could not set the regression mesh size or time step."
    exit 1
fi

if ! ./Allrun > "${ALLRUN_LOGFILE}" 2>&1; then
    echo "FAIL: Allrun exited non-zero. Surfacing logs:"
    dumpLogTail "Allrun" "${ALLRUN_LOGFILE}"
    for stage in blockMesh cardiacFoam; do
        dumpLogTail "${stage}" "log.${stage}"
    done
    exit 1
fi

if [[ ! -s "${SUMMARY_FILE}" ]] || ! grep -q 'Field     L1-error' "${SUMMARY_FILE}"; then
    echo "FAIL: manufactured electromechanics error summary not found in ${SUMMARY_FILE}"
    dumpLogTail "cardiacFoam" "log.cardiacFoam"
    exit 1
fi

echo "PASS: manufactured electromechanics error summary detected in ${SUMMARY_FILE}"

if ! checkReferenceValues "${SUMMARY_FILE}"; then
    echo
    echo "============================================================"
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
fi

echo
echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
exit 0
