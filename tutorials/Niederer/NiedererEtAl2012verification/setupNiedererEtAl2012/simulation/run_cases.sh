#!/bin/bash

set -euo pipefail

CASE_DIR="${1:-}"
if [[ -z "$CASE_DIR" ]]; then
    echo "Usage: ./run_cases.sh <caseDir>" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
COMMON_RUNNER="$REPO_ROOT/applications/scripts/driverFoam/openfoam_driver/scripts/run_case.sh"

if [[ ! -x "$COMMON_RUNNER" ]]; then
    echo "Common runner not found: $COMMON_RUNNER" >&2
    exit 1
fi

exec "$COMMON_RUNNER" \
    --case-dir "$CASE_DIR" \
    --parallel \
    --touch-case-foam
