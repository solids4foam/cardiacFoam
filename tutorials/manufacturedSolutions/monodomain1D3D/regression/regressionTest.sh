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

# Runs the case defaults (N 20, 21-node graph, unidirectional explicit PVJ) in
# serial and compares the 3D error summary, the PVJ source diagnostics and the
# graph error summary against the references.

regression_init "1D-3D monodomain manufactured-solution regression test" \
    regression/monodomain1D3D.reference "$@"

regression_run_or_fail

REGRESSION_SUMMARY_FILE="$(regression_resolve_file postProcessing/3D_20_cells.dat)" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "3D error summary" "${REGRESSION_SUMMARY_FILE}" 'Field     L1-error'
graphSummary="$(regression_resolve_file postProcessing/graph_1D_21_nodes.dat)" || graphSummary=""
regression_require_output "graph error summary" "${graphSummary}" 'Field     L1-error'
diagnostics="$(regression_resolve_file verification/coupled1D3DMonodomain_diagnostics.csv)" \
    || diagnostics=""
regression_require_output "PVJ source diagnostics" "${diagnostics}" 'sourceErrorL1'

regression_compare || true

# The Purkinje graph errors
REGRESSION_SUMMARY_FILE="${graphSummary}"
regression_compare regression/monodomain1D3D.graph.reference || true
regression_finish
