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

# The tutorial as it is (serial, about 10 s): ten S1 beats at 2000 ms, then
# two S2 stimuli at 250 ms intervals, Bueno-Orovio epicardial cell, Euler
# integrator. The first S2 falls in the refractory tail of the last S1 action
# potential and does not capture; the second does. Vm is checked on the tail
# of the last S1 beat (18.25 s), after the blocked S2 (18.3 s), and on the
# plateau and the repolarisation of the captured S2 beat (18.6 s, 18.7 s).

regression_init "S1-S2 restitution regression test" \
    regression/restitutionCurves_s1s2Protocol.reference "$@"
REGRESSION_TIME_WINDOW=1e-4

regression_run_or_fail
regression_compare || true
regression_finish
