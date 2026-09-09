#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Idealized heart injection regression test
# ============================================================
#
# Confirms the stimulus/Purkinje-to-myocardium injection is correct: probes
# activationTime at a Purkinje-myocardial-junction site on the LV free wall
# (node 211 in the shared purkinjeGraph) at t=0.02, the case's own endTime.
# This is the same point conductionBlock's lbbb regression checks stays
# un-activated (severed LV subtree) - here it must have activated.

REF_FILE="regression/injection.reference"
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

./Allclean > /dev/null 2>&1 || true
if ! ./Allrun > "${ALLRUN_LOGFILE}" 2>&1; then
    echo "FAIL: Allrun exited non-zero. Surfacing logs:"
    dumpLogTail "Allrun" "${ALLRUN_LOGFILE}"
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
    if [[ ! -f "${dataFile}" ]]; then
        echo "FAIL: missing output file ${dataFile}"
        failures=$((failures + 1))
        checks=$((checks + 1))
        continue
    fi

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

    checks=$((checks + 1))

    if [[ -z "${actual}" ]]; then
        echo "FAIL: ${dataFile} col=${column} at t=${time} not found"
        failures=$((failures + 1))
        continue
    fi

    diffAbs="$(
        awk -v a="${actual}" -v e="${expected}" '
            BEGIN {
                d = a - e;
                if (d < 0) d = -d;
                print d;
            }
        '
    )"

    if awk -v d="${diffAbs}" -v t="${tolerance}" 'BEGIN {exit !(d < t)}'; then
        printf "PASS: %s col=%s t=%s activationTime=%.7g (difference = %.3g)\n" \
            "${dataFile}" "${column}" "${time}" "${actual}" "${diffAbs}"
    else
        printf "FAIL: %s col=%s t=%s activationTime=%.7g (difference = %.3g)\n" \
            "${dataFile}" "${column}" "${time}" "${actual}" "${diffAbs}"
        failures=$((failures + 1))
    fi
done < "${REF_FILE}"

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures}/${checks} checks)"
    echo "============================================================"
    exit 1
fi
