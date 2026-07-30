#!/bin/bash
# Thin wrapper: sets this study's own values, execs the shared
# run_study_sweep_common.sh. All logic lives there -- do not add anything
# tutorial-specific here beyond these two variables.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

export SWEEP_SPEC="$SCRIPT_DIR/studies/spatialConvergence/sweep_spatial_convergence.json"
export AGG_KEY="eikonal_hex"
exec "$REPO_ROOT/applications/scripts/paperI_results/run_study_sweep_common.sh"
