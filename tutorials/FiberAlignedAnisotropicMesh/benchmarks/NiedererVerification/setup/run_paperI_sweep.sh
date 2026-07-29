#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"

for config in \
    paperI_dx05_implicit.json \
    paperI_dx02_explicit.json \
    paperI_dx01_explicit.json
do
    "$DRIVER" all \
        --entry niederer2012 \
        --config "$SCRIPT_DIR/$config"
done

echo "Paper I Niederer nine-case sweep complete."
