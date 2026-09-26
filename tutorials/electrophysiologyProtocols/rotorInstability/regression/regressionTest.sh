#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Shared regression helpers (tutorials/regressionFunctions)
helperDir="$(cd "${BASH_SOURCE[0]%/*}" && pwd)"
until [[ -f "${helperDir}/regressionFunctions" || "${helperDir}" == / ]]
do
    helperDir="$(dirname "${helperDir}")"
done
. "${helperDir}/regressionFunctions"

# ============================================================
# Rotor Instability activation regression test
# ============================================================

REF_FILE="regression/rotorInstability.reference"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Rotor Instability activation regression test"
echo "Activation times within the tolerances in ${REF_FILE}"
echo "============================================================"
echo

./Allclean > /dev/null 2>&1 || true
./Allrun parallel > "${ALLRUN_LOGFILE}" 2>&1
checkSolverLogs log.cardiacFoam log.reconstructPar log.postProcess || exit 1

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
