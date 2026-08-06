#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Bidomain manufactured-solution regression test
# ============================================================

ALLRUN_LOGFILE="log.Allrun"
BLOCKMESH_LOGFILE="log.blockMesh"
REF_FILE="regression/bidomainManufactured.reference"

echo "============================================================"
echo "Bidomain manufactured-solution regression test"
echo "Manufactured error summary must be present"
echo "============================================================"
echo

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

    candidate="$(findFirstMatch 'postProcessing/*.dat')" || true
    if [[ -n "${candidate}" ]] && grep -Eq '(^| )([Bb]idomain )?[Mm]anufactured-solution error summary' "${candidate}"; then
        echo "${candidate}"
        return 0
    fi

    for candidate in postProcessing/*.dat processor*/postProcessing/*.dat
    do
        if [[ -s "${candidate}" ]] && grep -Eq '(^| )([Bb]idomain )?[Mm]anufactured-solution error summary' "${candidate}"; then
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
        cells)
            awk -F= '/Number of cells/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "${summaryFile}"
            ;;
        finalTime)
            awk -F= '/Final simulation time/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "${summaryFile}"
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

    echo "Bidomain manufactured reference comparison: ${checks} checks, ${failures} failures"
    return "${failures}"
}

./Allclean > /dev/null 2>&1 || true
blockMesh -dict system/blockMeshDict.3D > "${BLOCKMESH_LOGFILE}" 2>&1
./Allrun > "${ALLRUN_LOGFILE}" 2>&1

errorFile="$(findManufacturedErrorFile)" || {
    echo "FAIL: manufactured error summary file not found."
    exit 1
}

if ! grep -q 'Field     L1-error' "${errorFile}"; then
    echo "FAIL: manufactured error table missing in ${errorFile}"
    exit 1
fi

echo "PASS: manufactured error summary detected in ${errorFile}"
checkReferenceValues "${errorFile}"
echo
echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
exit 0
