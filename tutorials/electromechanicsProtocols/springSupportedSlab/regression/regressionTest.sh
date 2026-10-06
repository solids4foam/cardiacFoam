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

# Spring-supported electromechanics slab. Needs the full solids4foam build:
# exits 77 under CARDIAC_REGRESSION_BUILD_MODE=lightweight.
#
# Besides the probe rows, the reference carries the case-local kind
#     springLaw <time> <tolerance [N]>
# which checks that the xMin force from the stress field equals
# -kEnds * A0 * <Dx> at that time, with kEnds read from system/caseParameters.

FORCE_FILE="postProcessing/0/solidForcesxMin.dat"
DISP_FILE="postProcessing/solid/D_xMin/0/surfaceFieldValue.dat"
END_AREA=2.1e-5    # end face, 3 mm x 7 mm (system/blockMeshDict)

regression_case_check()
{
    local kind="$1"
    local time="${2:-}" tolerance="${3:-}" force disp springForce

    if [[ "${kind}" != springLaw ]]; then
        regression_fail "unknown reference kind '${kind}'"
        return 1
    fi

    force="$(regression_probe_value "${FORCE_FILE}" "${time}" 2)" || force=""
    disp="$(regression_probe_value "${DISP_FILE}" "${time}" 2)" || disp=""
    if [[ -z "${force}" || -z "${disp}" ]]; then
        regression_fail "springLaw t=${time}: force or displacement not found"
        return 1
    fi

    springForce="$(awk -v k="${kEnds}" -v a="${END_AREA}" -v d="${disp}" 'BEGIN { print -k*a*d }')"
    regression_check "springLaw F_xMin vs -kEnds*A0*<Dx> t=${time}" \
        "${force}" "${springForce}" "${tolerance}" springLaw "${FORCE_FILE}" "${time}" 2
}

regression_init "Spring-supported electromechanics slab regression test" \
    regression/springSupportedSlab.reference "$@"
regression_require_solids4foam

# The tutorial runs the twitch to 0.25 s. The regression stops at 0.05 s, once
# activation has crossed the slab, Ta has risen in its middle and both
# spring-supported ends have moved well beyond their resting preload.
regression_set system/caseParameters endTime 0.05

regression_run_or_fail parallel

kEnds="$(awk '$1 == "kEnds" { sub(/;/, "", $2); print $2 }' system/caseParameters)"
if [[ -z "${kEnds}" ]]; then
    regression_fail "kEnds not found in system/caseParameters"
    regression_finish
fi

regression_compare || true
regression_finish
