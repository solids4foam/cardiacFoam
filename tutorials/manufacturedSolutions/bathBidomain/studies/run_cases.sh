#!/bin/bash

set -euo pipefail

# Usage: ./run_cases.sh <caseDir> <1D|2D|3D>
CASE_DIR="${1:-}"
DIM="${2:-}"

if [[ -z "$CASE_DIR" || -z "$DIM" ]]; then
    echo "Usage: ./run_cases.sh <caseDir> <1D|2D|3D>" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
COMMON_RUNNER="$REPO_ROOT/applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh"

exec "$COMMON_RUNNER" \
    --case-dir "$CASE_DIR" \
    --dimension "$DIM" \
    --parallel
