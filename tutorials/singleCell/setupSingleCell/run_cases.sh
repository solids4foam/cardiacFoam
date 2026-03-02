#!/bin/bash

set -euo pipefail

CASE_DIR="${1:-}"
if [[ -z "$CASE_DIR" ]]; then
    echo "Usage: ./run_cases.sh <caseDir>" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON_RUNNER="$SCRIPT_DIR/../../openfoam_driver/scripts/run_case.sh"

exec "$COMMON_RUNNER" \
    --case-dir "$CASE_DIR"


