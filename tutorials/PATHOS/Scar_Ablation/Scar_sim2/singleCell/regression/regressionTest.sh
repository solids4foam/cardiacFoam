#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Single-cell regression test
# ============================================================

VM_TOL=5e-3

REF_FILE="singleCell.reference"
DATA_FILE="postProcessing/AlievPanfilov_myocyte_S1_1000.txt"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Single-cell regression test"
echo "Vm difference < ${VM_TOL}"
echo "============================================================"
echo

./Allclean > /dev/null 2>&1 || true
./Allrun > "${ALLRUN_LOGFILE}" 2>&1

if [[ ! -f "${REF_FILE}" ]]; then
    echo "FAIL: reference file not found: ${REF_FILE}"
    exit 1
fi

if [[ ! -f "${DATA_FILE}" ]]; then
    echo "FAIL: output file not found: ${DATA_FILE}"
    exit 1
fi

failures=0
checks=0

while IFS=' ' read -r fileName time expected tolerance; do
    if [[ -z "${fileName}" || "${fileName}" == \#* ]]; then
        continue
    fi

    dataFile="${fileName}"
    if [[ ! -f "${dataFile}" ]]; then
        echo "FAIL: missing output file ${dataFile}"
        failures=$((failures + 1))
        checks=$((checks + 1))
        continue
    fi

    actual="$(
        awk -v target="${time}" '
            BEGIN { bestDiff = 1e99; found = 0; actual = 0.0; }
            NR > 1 && NF >= 2 {
                d = $1 - target;
                if (d < 0) d = -d;
                if (d < bestDiff) {
                    bestDiff = d;
                    actual = $2;
                    found = 1;
                }
            }
            END {
                if (found && bestDiff <= 1e-9) {
                    print actual;
                    exit 0;
                }
                exit 1;
            }
        ' "${dataFile}"
    )" || true

    checks=$((checks + 1))

    if [[ -z "${actual}" ]]; then
        echo "FAIL: ${dataFile} at t=${time} not found"
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
        printf "PASS: %s t=%s Vm=%.7g (difference = %.3g)\n" \
            "${dataFile}" "${time}" "${actual}" "${diffAbs}"
    else
        printf "FAIL: %s t=%s Vm=%.7g (difference = %.3g)\n" \
            "${dataFile}" "${time}" "${actual}" "${diffAbs}"
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
