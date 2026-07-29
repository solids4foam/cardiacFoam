#!/bin/bash
# Full manufactured-FDA bidomain sweep via the driverFoam engine.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"
"$DRIVER" all --entry manufacturedFDABidomain --config "$SCRIPT_DIR/../config/driver_config.json"
echo "bidomain manufacturedFDABidomain sweep complete."
