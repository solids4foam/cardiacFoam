#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Idealized heart conduction-block (lbbb) regression test
# ============================================================
#
# Confirms the LBB bridge severing actually blocks fast conduction: probes
# activationTime at the same Purkinje-myocardial-junction site (node 211 in
# the shared purkinjeGraph) that electrophysiologyHeart's own regression
# checks IS activated by t=0.02. Here, with the LV subtree disconnected
# from the root (a tree - no alternate path), it must stay un-activated
# (activationTime == -1) at the same time cutoff. Only the lbbb variant is
# covered - rbbb is exercised manually, not part of automated regression.
#
# conductionBlock's own controlDict runs to 0.7s (full ECG-scale, for real
# use of the tutorial), but this check only needs t=0.02 - the graph is a
# tree, so a severed subtree never receives Purkinje current regardless of
# how much longer the run continues. Running the full 0.7s here would only
# add ~45 minutes with no extra information for this check, so the run's
# own endTime is temporarily shortened to 0.02s just for this script's
# invocation, then restored - the tracked controlDict is never left
# changed. Allrun itself is reused unmodified.

REF_FILE="regression/lbbb.reference"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Idealized heart conduction-block (lbbb) regression test"
echo "============================================================"
echo

CONTROL_DICT="system/controlDict"
CONTROL_DICT_BACKUP="system/controlDict.regressionTest.bak"

cp "${CONTROL_DICT}" "${CONTROL_DICT_BACKUP}"
restoreControlDict()
{
    mv -f "${CONTROL_DICT_BACKUP}" "${CONTROL_DICT}"
}
trap restoreControlDict EXIT

sed -E 's/^endTime[[:space:]]+[^;]+;/endTime    0.02;/' \
    "${CONTROL_DICT_BACKUP}" > "${CONTROL_DICT}"

./Allclean > /dev/null 2>&1 || true
./Allrun lbbb > "${ALLRUN_LOGFILE}" 2>&1

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
