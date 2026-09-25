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
# Niederer activation regression test
# ============================================================

ACTIVATION_TOL=1e-4

REF_FILE="regression/NiedererEtAl2011.reference"
ALLRUN_LOGFILE="log.Allrun"
CHECK_ONLY=0
REPORT_PATH=""
REPORT_ROWS=""

# --check-only and --report exist for an out-of-tree consumer, not for this
# repository: omniDriver's cardiacfoam plugin runs the case itself and then
# calls this script as the case's own checker,
#
#     bash regression/regressionTest.sh --check-only --report regression/comparison-report.json
#
# and reads the JSON it writes. The split is the point: the solver side owns
# the pass/fail criteria -- this script and its reference file -- while the
# driver only runs them and repeats the verdict. Nothing in this repository
# passes either option; Alltest-regression calls the script with no arguments.
usage() {
    cat <<'EOF'
Usage: regressionTest.sh [--check-only] [--report PATH]

Without options, clean the case, run Allrun, and compare its outputs to the
solver-owned reference data.  --check-only performs the same comparison on
existing outputs without cleaning or running the case.  --report writes the
comparison evidence as JSON.
EOF
}

while (( $# > 0 )); do
    case "$1" in
        --check-only)
            CHECK_ONLY=1
            ;;
        --report)
            if (( $# < 2 )); then
                echo "FAIL: --report requires a path" >&2
                exit 2
            fi
            REPORT_PATH="$2"
            shift
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

json_string() {
    # Pure parameter expansion: BSD sed's N quits without printing on the last
    # line, so a sed-based escape returns an empty string on macOS.
    local s="$1"
    s=${s//\\/\\\\}
    s=${s//\"/\\\"}
    s=${s//$'\n'/\\n}
    printf '%s' "$s"
}

append_report_row() {
    if [[ -z "${REPORT_ROWS}" ]]; then
        return
    fi

    local result="$1"
    local reason="$2"
    local actual="$3"
    local difference="$4"
    printf '{"file":"%s","time":%s,"column":%s,"expected":%s,"tolerance":%s,"actual":%s,"difference":%s,"status":"%s","reason":"%s"}\n' \
        "$(json_string "${dataFile}")" "${time}" "${column}" \
        "${expected}" "${tolerance}" "${actual}" "${difference}" \
        "${result}" "${reason}" >> "${REPORT_ROWS}"
}

write_report() {
    if [[ -z "${REPORT_PATH}" ]]; then
        return
    fi

    mkdir -p "$(dirname "${REPORT_PATH}")"
    {
        printf '{\n'
        printf '  "schema_version": 1,\n'
        printf '  "mode": "%s",\n' "$([[ "${CHECK_ONLY}" -eq 1 ]] && printf 'check-only' || printf 'run-and-check')"
        printf '  "reference_file": "%s",\n' "$(json_string "${REF_FILE}")"
        printf '  "status": "%s",\n' "$([[ "${failures}" -eq 0 ]] && printf 'passed' || printf 'failed')"
        printf '  "checks": %s,\n' "${checks}"
        printf '  "failures": %s,\n' "${failures}"
        printf '  "results": ['
        if [[ -n "${REPORT_ROWS}" && -s "${REPORT_ROWS}" ]]; then
            paste -sd, "${REPORT_ROWS}"
        fi
        printf ']\n}\n'
    } > "${REPORT_PATH}"
}

if [[ -n "${REPORT_PATH}" ]]; then
    REPORT_ROWS="$(mktemp)"
    trap 'rm -f "${REPORT_ROWS}"' EXIT
fi

echo "============================================================"
echo "Niederer activation regression test"
echo "Activation-time difference < ${ACTIVATION_TOL}"
echo "============================================================"
echo

if [[ "${CHECK_ONLY}" -eq 1 ]]; then
    echo "Checking existing outputs (no cleanup or run)"
else
    ./Allclean > /dev/null 2>&1 || true
    ./Allrun parallel > "${ALLRUN_LOGFILE}" 2>&1
    checkSolverLogs log.cardiacFoam log.reconstructPar log.postProcess || exit 1
fi

if [[ ! -f "${REF_FILE}" ]]; then
    echo "FAIL: reference file not found: ${REF_FILE}"
    failures=1
    checks=0
    write_report
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
        append_report_row "failed" "missing-output" "null" "null"
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
        append_report_row "failed" "value-not-found" "null" "null"
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
        append_report_row "passed" "within-tolerance" "${actual}" "${diffAbs}"
    else
        printf "FAIL: %s col=%s t=%s activationTime=%.7g (difference = %.3g)\n" \
            "${dataFile}" "${column}" "${time}" "${actual}" "${diffAbs}"
        failures=$((failures + 1))
        append_report_row "failed" "outside-tolerance" "${actual}" "${diffAbs}"
    fi
done < "${REF_FILE}"

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    write_report
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures}/${checks} checks)"
    echo "============================================================"
    write_report
    exit 1
fi
