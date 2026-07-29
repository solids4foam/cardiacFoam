#!/bin/bash
# Shared hex-sweep runner for the four manufactured-solution tutorials
# (bidomain, monodomainPseudoECG, eikonalECG, bathBidomain). Each tutorial's
# own run_<x>_hex.sh is a thin wrapper that sets the three variables below
# and execs this file -- the logic itself (wipe the stale archive, sweep-run,
# aggregate) must be identical across all four. This is not a convenience:
# it's the whole point of the sweepCases pattern -- the action is agnostic to
# which tutorial it's running, so there is exactly one implementation to
# drift out of sync, not four hand-maintained copies that quietly diverge.
#
# Required environment variables (set by the calling wrapper):
#   TUTORIAL_DIR - the tutorial's case_root, e.g. .../manufacturedSolutions/bidomain
#   SWEEP_SPEC   - path to that tutorial's sweep_spatial_convergence.json
#   AGG_KEY      - aggregate.py case key, e.g. bidomain_hex
set -euo pipefail
: "${TUTORIAL_DIR:?TUTORIAL_DIR must be set by the calling run_<x>_hex.sh}"
: "${SWEEP_SPEC:?SWEEP_SPEC must be set by the calling run_<x>_hex.sh}"
: "${AGG_KEY:?AGG_KEY must be set by the calling run_<x>_hex.sh}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"

# One sweep-run call here is not one solve: SWEEP_SPEC expands into a LIST of
# resolved cases (1D/2D/3D x N=10/20/40/80, 12 in all for these four
# tutorials), each materialized and solved in turn, sequentially, in the
# shared case_root. Each case's raw output is archived into its own
# sweepCases/<case_id>/ subfolder for aggregate.py to read afterwards.
rm -rf "$TUTORIAL_DIR/sweepCases" "$TUTORIAL_DIR/sweepRun"
"$DRIVER" sweep-run --spec "$SWEEP_SPEC" --output-dir "$TUTORIAL_DIR/sweepRun"

python3 "$SCRIPT_DIR/aggregate.py" "$AGG_KEY" --repo-root "$REPO_ROOT" \
    || echo "WARN: paperI aggregate ($AGG_KEY) failed; native sweepCases output untouched" >&2

echo "$AGG_KEY sweep complete."
