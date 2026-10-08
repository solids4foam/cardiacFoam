#!/usr/bin/env bash
set -euo pipefail

# Shared regression library (tutorials/regression/lib.sh), found by walking up
# from this case; CARDIAC_REGRESSION_LIB overrides the lookup.
regressionLib="${CARDIAC_REGRESSION_LIB:-}"
if [[ -z "${regressionLib}" ]]; then
    dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    while [[ "${dir}" != / && ! -f "${dir}/regression/lib.sh" ]]; do
        dir="$(dirname "${dir}")"
    done
    regressionLib="${dir}/regression/lib.sh"
fi
. "${regressionLib}"

# 1D-3D monodomain manufactured solution at N_CELLS=20 with the default
# 41-node graph. Case-local reference kinds, each with a relative (rel) or
# absolute (abs) tolerance:
#     error3D  <field> <L1|L2|Linf> <expected> <tolerance> <rel|abs>
#         3-D error summary, postProcessing/3D_20_cells.dat
#     error1D  <field> <L1|L2|Linf> <expected> <tolerance> <rel|abs>
#         graph error summary, postProcessing/graph_1D_41_nodes.dat
#     coupling <column> - <expected> <tolerance> <rel|abs>
#         final row of the 1D-3D coupling diagnostics CSV

SUMMARY_3D="postProcessing/3D_20_cells.dat"
SUMMARY_1D="postProcessing/graph_1D_41_nodes.dat"
DIAGNOSTICS="verification/coupled1D3DMonodomain_diagnostics.csv"

# Final-row value of a named CSV column
csvFinalValue()
{
    awk -F, -v name="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) if ($i == name) col = i; next }
        { last = $col }
        END { if (col) print last }
    ' "$1"
}

regression_case_check()
{
    local kind="$1" key="${2:-}" metric="${3:-}" expected="${4:-}" tolerance="${5:-}" mode="${6:-abs}"
    local file actual tolAbs

    case "${kind}" in
        error3D)  file="${SUMMARY_3D}" ;;
        error1D)  file="${SUMMARY_1D}" ;;
        coupling) file="${DIAGNOSTICS}" ;;
        *) regression_fail "unknown reference kind '${kind}'"; return 1 ;;
    esac

    if [[ "${kind}" == coupling ]]; then
        actual="$(csvFinalValue "${file}" "${key}")" || actual=""
    else
        REGRESSION_SUMMARY_FILE="${file}"
        actual="$(regression_error_metric "${key}" "${metric}")" || actual=""
    fi

    if [[ "${mode}" == rel ]]; then
        tolAbs="$(awk -v t="${tolerance}" -v e="${expected}" 'BEGIN { if (e < 0) e = -e; print t*e }')"
    else
        tolAbs="${tolerance}"
    fi

    regression_check "${kind} ${key} ${metric} (${mode} ${tolerance})" "${actual}" "${expected}" "${tolAbs}" \
        "${kind}" "${file}" "" "${key} ${metric}"
}

regression_init "1D-3D monodomain manufactured-solution regression test" \
    regression/monodomain1D3D.reference "$@"

export N_CELLS=20
regression_run_or_fail

regression_require_output "3-D error summary" "${SUMMARY_3D}"
regression_require_output "graph error summary" "${SUMMARY_1D}"
regression_require_output "coupling diagnostics" "${DIAGNOSTICS}"

regression_compare || true
regression_finish
