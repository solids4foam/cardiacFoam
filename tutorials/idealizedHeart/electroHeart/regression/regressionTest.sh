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
# CARDIAC_REGRESSION_SCOPE=standard (default): the monodomain variant on the
# human Purkinje tree, regression/injection.monodomain.reference.
# CARDIAC_REGRESSION_SCOPE=full: every solver variant (monodomain, eikonal,
# hybrid) on both trees (human, pig). Each combination has its own reference,
# injection.<variant>.reference for human and injection.<variant>.pig.reference
# for pig: monodomain activates the probe at 28.9 ms, hybrid at 30.1 ms, and
# eikonal solves one steady problem that writes only time 1.

regression_init "Idealized heart injection regression test" \
    regression/injection.monodomain.reference "$@"
regression_scope

VARIANTS=(monodomain)
TREES=(human)
if [[ "${REGRESSION_SCOPE}" == full ]]; then
    VARIANTS=(monodomain eikonal hybrid)
    TREES=(human pig)
fi

for tree in "${TREES[@]}"; do
for variant in "${VARIANTS[@]}"; do
    refFile="regression/injection.${variant}.reference"
    if [[ "${tree}" == pig ]]; then
        refFile="regression/injection.${variant}.pig.reference"
    fi

    echo "------------------------------------------------------------"
    echo "Variant: ${variant}  Purkinje tree: ${tree}"
    echo "------------------------------------------------------------"

    if ! regression_run "${variant}" "${tree}" parallel; then
        regression_fail "Allrun ${variant} ${tree} parallel did not complete"
        continue
    fi
    regression_compare "${refFile}" || true
    echo
done
done

regression_finish
