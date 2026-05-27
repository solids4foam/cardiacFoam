#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Bath-bidomain manufactured-solution regression test
# ============================================================
#
# Uses the 1D mesh (system/blockMeshDict.1D, 80 cells per block) that
# matches the dimension "1D" entry in constant/electroProperties. The
# manufacturedFDABathBidomainVerifier writes a per-run summary file
# named postProcessing/bathBidomain_<DIM>_<N>_cells_<implicit|explicit>.dat
# whose values are compared against bathBidomainManufactured.reference.

ALLRUN_LOGFILE="log.Allrun"
REF_FILE="bathBidomainManufactured.reference"

echo "============================================================"
echo "Bath-bidomain manufactured-solution regression test"
echo "Mesh: system/blockMeshDict.1D, run mode: serial"
echo "============================================================"
echo

# Known-broken on OpenFOAM v2312 (both lightweight and with-solids4foam).
# v2312's dictionary lookup throws FATAL for
# 'laplacian(conductivityIntracellular,Vm)' against system/fvSchemes even
# when the literal entry is present and the heart sub-mesh is registered
# with the base mesh name (region0) and the field names exactly match.
# The same case passes on v2412 and v2512 in both modes. A diagnostic
# Info<< probe at extracellularPotentialDomain.C:647 confirmed the lookup
# key, sub-mesh name, and field names are correct; the bug is internal to
# v2312's schemesLookup machinery. Suppress the regression for v2312 only.
ofVersion="${WM_PROJECT_VERSION:-unknown}"
if [[ "${ofVersion}" == *2312* ]]; then
    echo "SKIP: bathBidomain regression is suppressed on OpenFOAM v2312."
    echo "      v2312 dictionary lookup of 'laplacian(conductivityIntracellular,Vm)'"
    echo "      fails for both lightweight and with-solids4foam modes; v2412 and"
    echo "      v2512 (both modes) exercise this test successfully."
    exit 0
fi

findFirstMatch()
{
    local pattern="$1"
    local candidate

    for candidate in ${pattern}
    do
        if [[ -s "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

findManufacturedErrorFile()
{
    local candidate

    for candidate in \
        postProcessing/bathBidomain_*_cells_*.dat \
        processor*/postProcessing/bathBidomain_*_cells_*.dat
    do
        if [[ -s "${candidate}" ]] \
            && grep -q 'Bath-bidomain manufactured solution error summary' "${candidate}"
        then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
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

extractSummaryValue()
{
    local summaryFile="$1"
    local key="$2"

    case "${key}" in
        cellsPerDirection)
            awk '/^# cellsPerDirection/ {print $3; exit}' "${summaryFile}"
            ;;
        finalTime)
            awk '/^# time/ {print $3; exit}' "${summaryFile}"
            ;;
        *)
            return 1
            ;;
    esac
}

extractErrorMetric()
{
    local summaryFile="$1"
    local field="$2"
    local metric="$3"
    local column

    case "${metric}" in
        L1)   column=2 ;;
        L2)   column=3 ;;
        Linf) column=4 ;;
        *) return 1 ;;
    esac

    awk -v field="${field}" -v column="${column}" \
        '$1 == field {print $column; exit}' "${summaryFile}"
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
                actual="$(extractSummaryValue "${summaryFile}" "${key}")"
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

    echo "Bath-bidomain manufactured reference comparison: ${checks} checks, ${failures} failures"
    return "${failures}"
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

./Allclean > /dev/null 2>&1 || true

if ! ./Allrun > "${ALLRUN_LOGFILE}" 2>&1; then
    echo "FAIL: Allrun exited non-zero. Surfacing logs:"
    dumpLogTail "Allrun" "${ALLRUN_LOGFILE}"
    for stage in blockMesh topoSet setTorsoOrganConductivityField decomposePar cardiacFoam reconstructPar; do
        dumpLogTail "${stage}" "log.${stage}"
    done
    exit 1
fi

errorFile="$(findManufacturedErrorFile)" || {
    echo "FAIL: bath-bidomain manufactured error summary file not found."
    echo "Surfacing recent logs to aid debugging:"
    dumpLogTail "Allrun" "${ALLRUN_LOGFILE}" 40
    dumpLogTail "cardiacFoam" "log.cardiacFoam"
    exit 1
}

if ! grep -q '^# field L1 L2 Linf' "${errorFile}"; then
    echo "FAIL: manufactured error table missing in ${errorFile}"
    exit 1
fi

echo "PASS: bath-bidomain manufactured error summary detected in ${errorFile}"
checkReferenceValues "${errorFile}"
echo
echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
exit 0
