#!/bin/bash
# Reproduce the N=10,20,40 predictor-corrector ladder reported in Paper I.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../.." && pwd)"

RESOLUTIONS="10 20 40" VARIANTS="baseline predictor" \
  bash "$SCRIPT_DIR/studies/coupling/run_coupling_study.sh"

python3 "$REPO_ROOT/applications/scripts/paperI_results/aggregate.py" bath_tet \
  --repo-root "$REPO_ROOT"
