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

# Sustained re-entry on the 2D slab: compares the activation-time probes
# written by the rotorpoints function object against the reference.

regression_init "Rotor instability activation regression test" \
    regression/rotorInstability.reference "$@"

# The tutorial runs 4 s. The S2 stimulus at 0.45 s starts the rotor, and by
# 1 s every probe has been activated again by the rotor alone (twice for
# some), so the regression stops there; the third stimulus, at 2 s, is not
# reached.
regression_set system/controlDict endTime 1

regression_run_or_fail parallel
regression_compare || true
regression_finish
