#!/bin/bash
# Shared sweep runner for every manufactured-solution study (hex or tet,
# any tutorial). Each tutorial's own run_<x>_<hex|tet>.sh is a thin wrapper
# that sets the two variables below and execs this file -- the logic itself
# (wipe the stale archive, sweep-run, aggregate) must be identical across
# every one of them. This is not a convenience: it's the whole point of the
# sweepCases pattern -- the action is agnostic to which tutorial or mesh
# family it's running, so there is exactly one implementation to drift out
# of sync, not one hand-maintained copy per study that quietly diverges.
#
# The study's own directory (setup/studies/<name>/, containing the sweep.json
# itself) is derived from SWEEP_SPEC's own location, not passed separately --
# sweepCases/ and sweepRun/ are siblings of the sweep.json that drives them,
# one self-contained directory per study.
#
# Required environment variables (set by the calling wrapper):
#   SWEEP_SPEC - path to the study's own sweep_*.json, e.g.
#                .../bidomain/setup/studies/spatialConvergence/sweep_spatial_convergence.json
#   AGG_KEY    - aggregate.py case key, e.g. bidomain_hex / bidomain_tet
set -euo pipefail
: "${SWEEP_SPEC:?SWEEP_SPEC must be set by the calling run_<x>_<hex|tet>.sh}"
: "${AGG_KEY:?AGG_KEY must be set by the calling run_<x>_<hex|tet>.sh}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"
STUDY_DIR="$(cd "$(dirname "$SWEEP_SPEC")" && pwd)"

# One sweep-run call here is not one solve: SWEEP_SPEC expands into a LIST of
# resolved cases, each materialized and solved in turn, sequentially, in the
# tutorial's shared case_root. Each case's raw output is archived into its
# own sweepCases/<case_id>/ subfolder (a sibling of the sweep.json, under
# STUDY_DIR) for aggregate.py to read afterwards.
rm -rf "$STUDY_DIR/sweepCases" "$STUDY_DIR/sweepRun"
"$DRIVER" sweep-run --spec "$SWEEP_SPEC" --output-dir "$STUDY_DIR/sweepRun"

python3 "$SCRIPT_DIR/aggregate.py" "$AGG_KEY" --repo-root "$REPO_ROOT" \
    || echo "WARN: paperI aggregate ($AGG_KEY) failed; native sweepCases output untouched" >&2

echo "$AGG_KEY sweep complete."
