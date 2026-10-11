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

# Runs the case defaults (3D hex, N 20, groundElectrode) in parallel and
# compares the summary values and error norms of the manufactured-solution
# error summary postProcessing/3D_20_cells.dat, and the interface rows of
# bathBidomainInterfaceMetrics (solved field) and
# bathBidomainInterfaceMetricsExactField (exact phiE), against the reference.

# interface|interfaceExact <column> <expected> <tolerance>: one column of the
# single data row of postProcessing/bathBidomainInterfaceMetrics[ExactField].csv
regression_case_check()
{
    local kind="$1" column="$2" expected="$3" tolerance="$4" file actual
    case "${kind}" in
        interface) file=postProcessing/bathBidomainInterfaceMetrics.csv ;;
        interfaceExact) file=postProcessing/bathBidomainInterfaceMetricsExactField.csv ;;
        *)
            regression_fail "unknown reference kind '${kind}'"
            return 1
            ;;
    esac
    actual="$(awk -F, -v name="${column}" '
        NR == 1 { for (i = 1; i <= NF; i++) if ($i == name) col = i; next }
        NR == 2 && col { print $col }' "${file}" 2>/dev/null)" || actual=""
    regression_check "${kind} ${column}" "${actual}" "${expected}" "${tolerance}" \
        "${kind}" "${file}" "" "${column}"
}

regression_init "Bath-bidomain manufactured-solution regression test" \
    regression/bathBidomainManufactured.reference "$@"

# OpenFOAM v2312's dictionary lookup of 'laplacian(conductivityIntracellular,Vm)'
# against system/fvSchemes fails for this case in both build modes, while
# v2412 and v2512 pass. The regression is suppressed on v2312 only.
if [[ "${WM_PROJECT_VERSION:-unknown}" == *2312* ]]; then
    echo "SKIP: bathBidomain regression is suppressed on OpenFOAM v2312 (regression/skips-on-openfoam)."
    exit "${REGRESSION_SKIP_CODE}"
fi

regression_run_or_fail parallel

REGRESSION_SUMMARY_FILE="$(regression_find_output \
    'Bath-bidomain manufactured solution error summary' 'postProcessing/*_cells.dat')" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "bath-bidomain manufactured error summary" \
    "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'

regression_compare || true
regression_finish
