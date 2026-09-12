#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Idealized heart injection regression test
# ============================================================
#
# Confirms the stimulus/Purkinje-to-myocardium injection is correct across
# all three solver variants: probes activationTime at a
# Purkinje-myocardial-junction site on the LV free wall (node 211 in the
# shared purkinjeGraph). This is the same point conductionBlock's lbbb
# regression checks stays un-activated (severed LV subtree) - here it must
# have activated, in every variant.
#
# Each variant has its own reference because each reaches the probe
# differently: monodomain activates it at 19.4ms, hybrid's eikonal-1D
# Purkinje to 3D monodomain coupling at 31.3ms, and eikonal solves a single
# steady problem that writes only time 1. The sample time in each reference
# reflects that; the expected values were measured, not chosen.

VARIANTS=(monodomain eikonal hybrid)
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Idealized heart injection regression test"
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

# Compare one reference file's rows against this variant's output.
# Echoes PASS/FAIL per row; returns the number of failures.
checkReference()
{
    local refFile="$1"
    local variantFailures=0

    while IFS=' ' read -r fileName time column expected tolerance; do
        if [[ -z "${fileName}" || "${fileName}" == \#* ]]; then
            continue
        fi

        local dataFile="postProcessing/${fileName}"
        if [[ ! -f "${dataFile}" ]]; then
            echo "FAIL: missing output file ${dataFile}"
            variantFailures=$((variantFailures + 1))
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
                    if (found && bestDiff <= 2.5e-3) {
                        print actual;
                        exit 0;
                    }
                    exit 1;
                }
            ' "${dataFile}"
        )" || true

        if [[ -z "${actual}" ]]; then
            echo "FAIL: ${dataFile} col=${column} at t=${time} not found"
            variantFailures=$((variantFailures + 1))
            continue
        fi

        local diffAbs
        diffAbs="$(
            awk -v a="${actual}" -v e="${expected}" \
                'BEGIN { d = a - e; if (d < 0) d = -d; print d; }'
        )"

        if awk -v d="${diffAbs}" -v t="${tolerance}" 'BEGIN {exit !(d < t)}'; then
            printf "PASS: %s col=%s t=%s activationTime=%.7g (difference = %.3g)\n" \
                "${dataFile}" "${column}" "${time}" "${actual}" "${diffAbs}"
        else
            printf "FAIL: %s col=%s t=%s activationTime=%.7g (difference = %.3g)\n" \
                "${dataFile}" "${column}" "${time}" "${actual}" "${diffAbs}"
            variantFailures=$((variantFailures + 1))
        fi
    done < "${refFile}"

    return "${variantFailures}"
}

failures=0
failedVariants=()

for variant in "${VARIANTS[@]}"; do
    refFile="regression/injection.${variant}.reference"

    echo "------------------------------------------------------------"
    echo "Variant: ${variant}"
    echo "------------------------------------------------------------"

    if [[ ! -f "${refFile}" ]]; then
        echo "FAIL: reference file not found: ${refFile}"
        failures=$((failures + 1))
        failedVariants+=("${variant}")
        echo
        continue
    fi

    ./Allclean > /dev/null 2>&1 || true

    if ! ./Allrun "${variant}" > "${ALLRUN_LOGFILE}" 2>&1; then
        echo "FAIL: Allrun ${variant} exited non-zero. Surfacing logs:"
        dumpLogTail "Allrun" "${ALLRUN_LOGFILE}"
        dumpLogTail "cardiacFoam" "log.cardiacFoam"
        failures=$((failures + 1))
        failedVariants+=("${variant}")
        echo
        continue
    fi

    variantFailures=0
    checkReference "${refFile}" || variantFailures=$?

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
