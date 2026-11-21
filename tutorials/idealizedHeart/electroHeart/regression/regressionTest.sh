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

# Probes activationTime at a Purkinje-myocardial-junction site on the LV free
# wall (node 211 of the shared purkinjeGraph) and samples the pseudo-ECG. The
# same point is the one conductionBlock's lbbb check requires to stay
# un-activated; here it must have activated.
#
# Full tier (regression/full-only). The three solver variants, each with
# its own solver path, on the human Purkinje tree; the pig tree is only a
# different input and is not run. Each variant has its own reference,
# regression/injection.<variant>.reference: monodomain activates the probe at
# 28.6 ms, hybrid at 28.9 ms, and eikonal solves one steady problem that
# writes only time 1.

regression_init "Idealized heart injection regression test" \
    regression/injection.monodomain.reference "$@"
regression_require_run_mode

VARIANTS=(monodomain eikonal hybrid)
tree=human

for variant in "${VARIANTS[@]}"; do
    refFile="regression/injection.${variant}.reference"

    echo "------------------------------------------------------------"
    echo "Variant: ${variant}  Purkinje tree: ${tree}"
    echo "------------------------------------------------------------"

    if ! regression_run "${variant}" "${tree}" parallel; then
        regression_fail "Allrun ${variant} ${tree} parallel did not complete"
        continue
    fi
    # A run that reaches the outer-iteration cap still prints End
    if [[ "${variant}" == eikonal ]] \
        && ! regression_check_solver_logs -e 'PIMPLE: converged'; then
        regression_fail "eikonal ${tree}: the outer iteration did not converge"
    fi
    regression_compare "${refFile}" || true
    echo
done

regression_finish
