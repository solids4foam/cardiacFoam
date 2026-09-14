#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Static immersed-cylinder smoke regression test
#
# Runs the case and asserts that the immersed-boundary machinery
# is healthy. This is a smoke test, not a validation: see
# README.md for what it deliberately does not check.
# ============================================================

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Immersed cylinder (static) smoke regression test"
echo "============================================================"
echo

./Allclean > /dev/null 2>&1 || true

if ! ./Allrun > "${ALLRUN_LOGFILE}" 2>&1
then
    echo "FAIL: Allrun failed; last 40 lines of ${ALLRUN_LOGFILE}:"
    tail -40 "${ALLRUN_LOGFILE}"
    exit 1
fi

echo "Checking smoke criteria:"
if ! python3 regression/checkSmoke.py
then
    echo
    echo "FAIL: smoke criteria not met; last 40 lines of log.pimpleHFDIBFoam:"
    tail -40 log.pimpleHFDIBFoam 2>/dev/null || true
    exit 1
fi

echo
echo "Immersed cylinder smoke regression test passed"
