#!/usr/bin/env bash
set -euo pipefail
regressionLib="${CARDIAC_REGRESSION_LIB:-}"
if [[ -z "$regressionLib" ]]; then
    dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    while [[ "$dir" != / && ! -f "$dir/regression/lib.sh" ]]; do dir="$(dirname "$dir")"; done
    regressionLib="$dir/regression/lib.sh"
fi
. "$regressionLib"
regression_init "Idealized-heart ECG recording contracts" "" "$@"
regression_run_or_fail monodomain parallel
regression_require_output "matched samples" postProcessing/matchedLeadFieldECG.dat '^#.*time'
regression_require_output "raw torso samples" postProcessing/matchedTorsoPotentials.dat '^#.*time'
regression_require_output "online reciprocal samples" postProcessing/leadFieldECG.dat '^#.*time'
regression_require_output "online pseudo samples" postProcessing/pseudoECG.dat '^#.*time'
solverStatus=0
regression_check_solver_logs log.cardiacFoam log.replayLeadFieldECG || solverStatus=$?
regression_check "source/replay log completion (exit status)" "$solverStatus" 0 0 \
    contract log.replayLeadFieldECG || true
assessmentStatus=0
python3 setup/checkMatchedECG.py --online || assessmentStatus=$?
regression_check "source, Wilson and online/replay contracts (exit status)" \
    "$assessmentStatus" 0 0 contract postProcessing/matchedECGAssessment.json || true
# This reports the forward/reciprocal difference; no exact-conductor accuracy
# bound or physiological validation is inferred from an execution-contract pass.
regression_finish
