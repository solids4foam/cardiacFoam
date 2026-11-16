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

# Runs the case defaults (slab, hex, N 10, monodomain) and the layered torso
# with a current (-problem layered, the same mesh resolution), each in
# parallel on two ranks, and compares the rows of the lead-field ECG
# verifier summary postProcessing/manufacturedLeadFieldECGSummary_ECG.dat
# against regression/leadFieldECG.reference and
# regression/leadFieldECGLayered.reference: the reciprocity gap
# leadField-direct and the static errors against the exact potential. The
# layered run's current crosses the heart wall and the torso layer
# interfaces.

# leadFieldStencil <probes|entries> <expected> <tolerance>: the probed
# explicit operator of buildLeadFields (log.buildLeadFields, "leadFieldAdjoint:
# ... <n> probes, <n> entries"). With sealedWallLeastSquares only the wall
# cells take point-neighbour rows; a point row on every cell stores the zero
# entries of all the others.
regression_case_check()
{
    local kind="$1" column="$2" expected="$3" tolerance="$4" file actual
    case "${kind}:${column}" in
        leadFieldStencil:probes|leadFieldStencil:entries) ;;
        *)
            regression_fail "unknown reference kind '${kind} ${column}'"
            return 1
            ;;
    esac
    file=log.buildLeadFields
    actual="$(sed -n "s/^leadFieldAdjoint: .* \([0-9][0-9]*\) probes, \([0-9][0-9]*\) entries\$/\1 \2/p" \
        "${file}" 2>/dev/null | tail -1)" || actual=""
    if [[ "${column}" == probes ]]; then
        actual="${actual%% *}"
    else
        actual="${actual##* }"
    fi
    regression_check "${kind} ${column}" "${actual}" "${expected}" "${tolerance}" \
        "${kind}" "${file}" "" "${column}"
}

regression_init "Lead-field ECG manufactured-solution regression test" \
    regression/leadFieldECG.reference "$@"
regression_require_solids4foam
regression_require_run_mode

for problem in slab layered; do
    echo "------------------------------------------------------------"
    echo "Problem: ${problem}"
    echo "------------------------------------------------------------"

    if ! regression_run parallel -problem "${problem}"; then
        regression_fail "Allrun parallel -problem ${problem} did not complete"
        continue
    fi

    REGRESSION_SUMMARY_FILE="$(regression_find_output \
        'Manufactured lead-field ECG summary' \
        'postProcessing/manufacturedLeadFieldECGSummary_*.dat')" \
        || REGRESSION_SUMMARY_FILE=""
    regression_require_output "lead-field ECG verifier summary" \
        "${REGRESSION_SUMMARY_FILE}" '# row L1 L2 Linf'

    if [[ "${problem}" == slab ]]; then
        regression_compare regression/leadFieldECG.reference || true
    else
        regression_compare regression/leadFieldECGLayered.reference || true
    fi
    echo
done

regression_finish
