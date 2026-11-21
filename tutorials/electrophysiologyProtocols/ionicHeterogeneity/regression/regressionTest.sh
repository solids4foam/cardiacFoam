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

# ionicHeterogeneityProbe writes comma-separated tables keyed by the
# transmural coordinate, hence the csv rows and the 1e-4 key window.
#
# Two runs, both seconds: the probe as configured
# (regression/ionicHeterogeneity.reference), then with
# constant/ionicConstantOverrides (./Allrun overrides,
# regression/overrides.reference). The second is the suite's check of
# ionicConstantOverrides: a global scale, a set and a scale on named regions.
# Its log must report each of the three.

regression_init "ionicHeterogeneity regression test" \
    regression/ionicHeterogeneity.reference "$@"
regression_require_run_mode
REGRESSION_TIME_WINDOW=1e-4

for variant in default overrides; do
    echo "------------------------------------------------------------"
    echo "Variant: ${variant}"
    echo "------------------------------------------------------------"
    if [[ "${variant}" == default ]]; then
        regression_run_or_fail
        regression_compare regression/ionicHeterogeneity.reference || true
        continue
    fi
    regression_run_or_fail overrides
    for applied in 'global: BuenoOrovio scaled tauSi' \
                   'mCells: BuenoOrovio set tauSo1' \
                   'epicardialCells: BuenoOrovio scaled tauW1Minus'; do
        grep -q "ionicConstantOverrides.${applied}" log.ionicHeterogeneityProbe \
            || regression_fail "overrides: log.ionicHeterogeneityProbe does not report ${applied}"
    done
    regression_compare regression/overrides.reference || true
done

regression_finish
