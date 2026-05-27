#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
CASE_ROOT="$REPO_ROOT/tutorials/manufacturedSolutions/bathBidomain"
CONFIG="$SCRIPT_DIR/driver_config.json"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"
ACTION="${1:-all}"

if [[ "$ACTION" != "sim" && "$ACTION" != "all" ]]
then
    echo "Usage: $0 [sim|all]" >&2
    exit 1
fi

if [[ -z "${WM_PROJECT_DIR:-}" ]]
then
    echo "OpenFOAM is not sourced. Run: source /Volumes/OpenFOAM-v2412/etc/bashrc" >&2
    exit 1
fi

rm -rf \
    "$CASE_ROOT/archivedPostProcessing" \
    "$CASE_ROOT/postProcessing" \
    "$CASE_ROOT/logs"

exec "$DRIVER" "$ACTION" \
    --entry manufacturedFDABathBidomain \
    --config "$CONFIG"
