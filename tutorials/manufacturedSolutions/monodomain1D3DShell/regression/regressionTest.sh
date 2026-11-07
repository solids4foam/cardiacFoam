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

# Runs the case defaults (gmsh tetrahedra N 10 seed 1, two terminals on the
# endocardium, bidirectional transfer, continuum reference, linear kernel) in
# serial and compares the shell error summary, the graph error summary, the
# junction rows of each terminal and the linear-kernel map against the
# references. The mesh comes from gmsh, so the case skips where it is missing
# (regression/needs-tools).
if ! command -v gmsh > /dev/null; then
    echo "SKIP: this case needs gmsh"
    exit "${REGRESSION_SKIP_CODE}"
fi

regression_init "1D-3D monodomain sealed-shell manufactured-solution regression test" \
    regression/monodomain1D3DShell.reference "$@"

regression_run_or_fail

REGRESSION_SUMMARY_FILE="$(regression_resolve_file postProcessing/3D_tet_10_cells.dat)" \
    || REGRESSION_SUMMARY_FILE=""
regression_require_output "shell error summary" "${REGRESSION_SUMMARY_FILE}" '# field L1 L2 Linf'
graphSummary="$(regression_resolve_file postProcessing/graph_1D_11_nodes.dat)" || graphSummary=""
regression_require_output "graph error summary" "${graphSummary}" 'Field     L1-error'
localTable="$(regression_resolve_file verification/coupled1D3DMonodomain_local.csv)" || localTable=""
regression_require_output "junction local errors" "${localTable}" 'currentError'
mapTable="$(regression_resolve_file verification/coupled1D3DMonodomain_map.csv)" || mapTable=""
regression_require_output "junction map errors" "${mapTable}" 'normalisationError'

regression_compare || true

# The Purkinje graph errors
REGRESSION_SUMMARY_FILE="${graphSummary}"
regression_compare regression/monodomain1D3DShell.graph.reference || true

regression_finish
