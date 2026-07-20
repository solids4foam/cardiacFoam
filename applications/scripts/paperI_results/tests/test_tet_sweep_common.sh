#!/usr/bin/env bash
# Unit tests for the OF-free helpers in tet_sweep_common.sh (mms_lc, dt_for_n).
# Run: bash tests/test_tet_sweep_common.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
. "$HERE/../tet_sweep_common.sh"

fail() { echo "FAIL: $1"; exit 1; }

# mms_lc must NOT reproduce the %.17g artifact (0.025000000000000001).
[ "$(mms_lc 40)" = "0.025" ]  || fail "mms_lc 40 -> $(mms_lc 40) (expected 0.025)"
[ "$(mms_lc 10)" = "0.1" ]    || fail "mms_lc 10 -> $(mms_lc 10) (expected 0.1)"
[ "$(mms_lc 20)" = "0.05" ]   || fail "mms_lc 20 -> $(mms_lc 20) (expected 0.05)"

# dt_for_n table is the single source of the transient time step.
[ "$(dt_for_n 10)" = "0.00892857" ]   || fail "dt_for_n 10 -> $(dt_for_n 10)"
[ "$(dt_for_n 40)" = "0.000560538" ]  || fail "dt_for_n 40 -> $(dt_for_n 40)"
[ "$(dt_for_n 80)" = "0.000140174" ]  || fail "dt_for_n 80 -> $(dt_for_n 80)"

# Unknown N must fail (return non-zero), not silently echo garbage.
if dt_for_n 999 >/dev/null 2>&1; then fail "dt_for_n 999 should return non-zero"; fi

echo "PASS (mms_lc, dt_for_n)"
