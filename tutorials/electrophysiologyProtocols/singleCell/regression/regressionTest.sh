#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Single-cell regression test
# ============================================================

VM_TOL=5e-3

REF_FILE="regression/singleCell.reference"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Single-cell regression test"
echo "Comparing variables from reference file"
echo "============================================================"
echo

./Allclean > /dev/null 2>&1 || true
./Allrun > "${ALLRUN_LOGFILE}" 2>&1

if [[ ! -f "${REF_FILE}" ]]; then
    echo "FAIL: reference file not found: ${REF_FILE}"
    exit 1
fi



failures=0
checks=0

while IFS=' ' read -r fileName time variable expected tolerance; do
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
        awk -v target="${time}" -v varName="${variable}" '
            BEGIN { bestDiff = 1e99; found = 0; actual = 0.0; col = 0; }
            NR == 1 {
                # Look for column index in header
                for(i=1; i<=NF; i++) {
                    if($i == varName) { col = i; }
                }
                next;
            }
            col > 0 && NF >= col {
                d = $1 - target;
                if (d < 0) d = -d;
                if (d < bestDiff) {
                    bestDiff = d;
                    actual = $col;
                    found = 1;
                }
            }
            END {
                if (found && bestDiff <= 1e-4) {
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
        printf "PASS: %s %s t=%s val=%.7g (diff = %.3g)\n" \
            "${dataFile}" "${variable}" "${time}" "${actual}" "${diffAbs}"
    else
        printf "FAIL: %s %s t=%s val=%.7g (diff = %.3g)\n" \
            "${dataFile}" "${variable}" "${time}" "${actual}" "${diffAbs}"
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
