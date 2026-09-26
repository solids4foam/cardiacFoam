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
# 1D-3D monodomain manufactured-solution regression test
#
# Runs the coupled case at N_CELLS=20 with the default 41-node graph and
# compares, with regression/monodomain1D3D.reference:
#   error3D  <field> <L1|L2|Linf>  3-D error summary (3D_20_cells.dat)
#   error1D  <field> <L1|L2|Linf>  graph error summary (graph_1D_41_nodes.dat)
#   coupling <column>              final row of the 1D-3D coupling
#                                  diagnostics (coupled1D3DMonodomain_...csv)
# Each reference line gives the expected value and a relative (rel) or
# absolute (abs) tolerance.
# ============================================================

REF_FILE="regression/monodomain1D3D.reference"
SUMMARY_3D="postProcessing/3D_20_cells.dat"
SUMMARY_1D="postProcessing/graph_1D_41_nodes.dat"
DIAGNOSTICS="verification/coupled1D3DMonodomain_diagnostics.csv"

echo "============================================================"
echo "1D-3D monodomain manufactured-solution regression test"
echo "============================================================"
echo

./Allclean > /dev/null 2>&1 || true
N_CELLS=20 ./Allrun > log.Allrun 2>&1
checkSolverLogs log.cardiacFoam log.blockMesh || exit 1

for f in "${SUMMARY_3D}" "${SUMMARY_1D}" "${DIAGNOSTICS}"; do
    if [[ ! -s "${f}" ]]; then
        echo "FAIL: missing output file ${f}"
        exit 1
    fi
done

# Error table value: <file> <field> <L1|L2|Linf>
errorValue()
{
    local column
    case "$3" in
        L1) column=2 ;;
        L2) column=3 ;;
        Linf) column=4 ;;
        *) return 1 ;;
    esac
    awk -v field="$2" -v column="${column}" \
        '$1 == field {print $column; exit}' "$1"
}

# Final-row value of a named CSV column
csvValue()
{
    awk -F, -v name="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) if ($i == name) col = i; next }
        { last = $col }
        END { if (col) print last }
    ' "$1"
}

failures=0
checks=0

while IFS=' ' read -r kind key metric expected tolerance mode; do
    [[ -z "${kind}" || "${kind}" == \#* ]] && continue

    case "${kind}" in
        error3D) actual="$(errorValue "${SUMMARY_3D}" "${key}" "${metric}")" ;;
        error1D) actual="$(errorValue "${SUMMARY_1D}" "${key}" "${metric}")" ;;
        coupling) actual="$(csvValue "${DIAGNOSTICS}" "${key}")" ;;
        *)
            echo "FAIL: unknown reference kind '${kind}'"
            failures=$((failures + 1))
            continue
            ;;
    esac

    checks=$((checks + 1))
    label="${kind} ${key} ${metric}"

    if [[ -z "${actual}" ]]; then
        echo "FAIL: could not extract ${label}"
        failures=$((failures + 1))
        continue
    fi

    if awk -v a="${actual}" -v e="${expected}" -v t="${tolerance}" \
        -v mode="${mode}" 'BEGIN {
            d = a - e; if (d < 0) d = -d;
            s = e; if (s < 0) s = -s;
            exit !(mode == "abs" ? d <= t : d <= t*s)
        }'
    then
        echo "PASS: ${label} = ${actual} (expected ${expected}, ${mode} tol ${tolerance})"
    else
        echo "FAIL: ${label} = ${actual} (expected ${expected}, ${mode} tol ${tolerance})"
        failures=$((failures + 1))
    fi
done < "${REF_FILE}"

echo
echo "1D-3D reference comparison: ${checks} checks, ${failures} failures"
if (( failures > 0 )); then
    echo "Regression test FAILED"
    exit 1
fi
echo "Regression test PASSED"
exit 0
