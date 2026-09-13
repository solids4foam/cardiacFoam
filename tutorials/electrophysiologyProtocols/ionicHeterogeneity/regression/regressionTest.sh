#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# ionicHeterogeneity regression test
#
# Compares postProcessing/ionicHeterogeneityProbe/*.csv columns from this
# case against regression/ionicHeterogeneity.reference. Unlike the
# whitespace-delimited postProcessing files used elsewhere in this repo,
# ionicHeterogeneityProbe writes comma-separated CSVs, so columns are
# matched with awk -F,.
# ============================================================

REF_FILE="regression/ionicHeterogeneity.reference"
ALLRUN_LOGFILE="log.Allrun"
CHECK_ONLY=0

usage() {
    cat <<'EOF'
Usage: regressionTest.sh [--check-only]

Without options, clean the case, run Allrun, and compare its outputs to the
reference data. --check-only performs the same comparison on existing
outputs without cleaning or running the case.
EOF
}

while (( $# > 0 )); do
    case "$1" in
        --check-only)
            CHECK_ONLY=1
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "FAIL: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

echo "============================================================"
echo "ionicHeterogeneity regression test"
echo "Comparing variables from reference file"
echo "============================================================"
echo

if [[ "${CHECK_ONLY}" -eq 1 ]]; then
    echo "Checking existing outputs (no cleanup or run)"
else
    ./Allclean > /dev/null 2>&1 || true
    ./Allrun > "${ALLRUN_LOGFILE}" 2>&1
fi

if [[ ! -f "${REF_FILE}" ]]; then
    echo "FAIL: reference file not found: ${REF_FILE}"
    exit 1
fi

failures=0
checks=0

while IFS=' ' read -r dataFile key variable expected tolerance; do
    if [[ -z "${dataFile}" || "${dataFile}" == \#* ]]; then
        continue
    fi

    if [[ ! -f "${dataFile}" ]]; then
        echo "FAIL: missing output file ${dataFile}"
        failures=$((failures + 1))
        checks=$((checks + 1))
        continue
    fi

    actual="$(
        awk -F, -v target="${key}" -v varName="${variable}" '
            BEGIN { bestDiff = 1e99; found = 0; actual = 0.0; col = 0; }
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == varName) { col = i; }
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
        echo "FAIL: ${dataFile} ${variable} at key=${key} not found"
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
        printf "PASS: %s %s key=%s val=%.7g (diff = %.3g)\n" \
            "${dataFile}" "${variable}" "${key}" "${actual}" "${diffAbs}"
    else
        printf "FAIL: %s %s key=%s val=%.7g (diff = %.3g)\n" \
            "${dataFile}" "${variable}" "${key}" "${actual}" "${diffAbs}"
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
