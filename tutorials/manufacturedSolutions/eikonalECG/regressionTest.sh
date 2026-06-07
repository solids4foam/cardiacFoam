#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Eikonal ECG manufactured-solution regression test
# ============================================================

ALLRUN_LOGFILE="log.Allrun"
BLOCKMESH_LOGFILE="log.blockMesh"
EXPECTED_ELECTRODES=(E1 E2 E3 E4 E5)
REF_FILE="eikonalECG.reference"

echo "============================================================"
echo "Eikonal ECG manufactured-solution regression test"
echo "Manufactured field and ECG outputs must be present"
echo "Mesh: system/blockMeshDict.3D, run mode: parallel"
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
    if [[ -n "${candidate}" ]] && grep -q 'Eikonal manufactured activation-time summary' "${candidate}"; then
        echo "${candidate}"
        return 0
    fi

    for candidate in postProcessing/*.dat processor*/postProcessing/*.dat
    do
        if [[ -s "${candidate}" ]] && grep -q 'Eikonal manufactured activation-time summary' "${candidate}"; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

findPseudoECGFile()
{
    findFirstMatch 'postProcessing/eikonalECG.dat' \
        || findFirstMatch 'processor*/postProcessing/eikonalECG.dat'
}

findManufacturedPseudoECGSummary()
{
    findFirstMatch 'postProcessing/manufacturedEikonalECGSummary.dat' \
        || findFirstMatch 'processor*/postProcessing/manufacturedEikonalECGSummary.dat'
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
            awk '/^[[:space:]]*cells:/ {print $2; exit}' "log.blockMesh"
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

extractFinalPseudoECGValue()
{
    local dataFile="$1"
    local electrode="$2"

    awk -v electrode="${electrode}" '
        NR == 1 && $1 == "#" {
            for (i = 1; i <= NF; i++) {
                if ($i == "numeric_"electrode || $i == electrode) {
                    column = i - 1;
                }
            }
            next;
        }
        $1 !~ /^#/ && column > 0 && NF >= column {
            value = $column;
        }
        END {
            if (value != "") {
                print value;
                exit 0;
            }
            exit 1;
        }
    ' "${dataFile}"
}

checkReferenceValues()
{
    local summaryFile="$1"
    local pseudoECGFile="$2"
    local errorFile="$3"
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
                actual="$(extractErrorMetric "${errorFile}" "${key}" "${metric}")"
                ;;
            pseudoECG)
                actual="$(extractFinalPseudoECGValue "${pseudoECGFile}" "${key}")"
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

    echo "Eikonal ECG reference comparison: ${checks} checks, ${failures} failures"
    return "${failures}"
}

checkElectrodeConfiguration()
{
    local electroProperties="constant/electroProperties"
    local electrode

    if [[ ! -f "${electroProperties}" ]]; then
        echo "FAIL: electroProperties not found."
        return 1
    fi

    if ! grep -q 'ecgDomains' "${electroProperties}" || ! grep -q 'electrodePositions' "${electroProperties}"; then
        echo "FAIL: pseudo-ECG electrode configuration not found in ${electroProperties}"
        return 1
    fi

    for electrode in "${EXPECTED_ELECTRODES[@]}"
    do
        if ! grep -Eq "^[[:space:]]*${electrode}[[:space:]]*\\(" "${electroProperties}"; then
            echo "FAIL: electrode ${electrode} not found in ${electroProperties}"
            return 1
        fi
    done

    echo "PASS: pseudo-ECG electrode configuration contains E1 E2 E3 E4 E5"
    return 0
}

checkPseudoECGHeader()
{
    local dataFile="$1"
    local electrode

    if ! grep -q '^#.*time' "${dataFile}"; then
        echo "FAIL: eikonalECG header missing in ${dataFile}"
        return 1
    fi

    for electrode in "${EXPECTED_ELECTRODES[@]}"
    do
        if ! awk -v electrode="${electrode}" '
            NR == 1 && $1 == "#" {
                for (i = 1; i <= NF; i++) {
                    if ($i == "numeric_"electrode || $i == electrode) {
                        found = 1;
                    }
                }
            }
            END { exit !found }
        ' "${dataFile}"; then
            echo "FAIL: electrode ${electrode} missing from eikonalECG header in ${dataFile}"
            return 1
        fi
    done

    echo "PASS: eikonalECG output header contains E1 E2 E3 E4 E5"
    return 0
}

./Allclean > /dev/null 2>&1 || true
checkElectrodeConfiguration
blockMesh -dict system/blockMeshDict.3D > "${BLOCKMESH_LOGFILE}" 2>&1
./Allrun parallel > "${ALLRUN_LOGFILE}" 2>&1

errorFile="$(findManufacturedErrorFile)" || {
    echo "FAIL: manufactured error summary file not found."
    exit 1
}

if ! grep -q 'activationTime' "${errorFile}"; then
    echo "FAIL: manufactured error table missing activationTime in ${errorFile}"
    exit 1
fi

echo "PASS: manufactured error summary detected in ${errorFile}"

pseudoECGFile="$(findPseudoECGFile)" || {
    echo "FAIL: eikonalECG.dat not found in postProcessing/ or processor*/postProcessing/."
    exit 1
}

checkPseudoECGHeader "${pseudoECGFile}"

if ! awk '
    $1 !~ /^#/ && NF >= 2 { found = 1; exit 0 }
    END { exit !found }
' "${pseudoECGFile}"; then
    echo "FAIL: no numeric eikonalECG samples found in ${pseudoECGFile}"
    exit 1
fi

echo "PASS: eikonalECG output detected in ${pseudoECGFile}"

summaryFile="$(findManufacturedPseudoECGSummary)" || {
    echo "FAIL: manufacturedEikonalECGSummary.dat not found."
    exit 1
}

if ! grep -q '^Manufactured eikonal ECG summary' "${summaryFile}"; then
    echo "FAIL: manufactured eikonal ECG summary header missing in ${summaryFile}"
    exit 1
fi

echo "PASS: manufactured eikonal ECG summary detected in ${summaryFile}"

checkReferenceValues "${summaryFile}" "${pseudoECGFile}" "${errorFile}"
echo
echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
exit 0
