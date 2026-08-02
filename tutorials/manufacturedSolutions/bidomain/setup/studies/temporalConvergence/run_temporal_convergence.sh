#!/bin/bash
# ---------------------------------------------------------------------------
# Bidomain Temporal Convergence Sweep
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

_PAPERI_ROOT="$(cd "$SCRIPT_DIR/../../../../../.." && pwd)"
DRIVER_FOAM_BIN="$_PAPERI_ROOT/applications/scripts/driverFoam/bin/driverFoam"

if [[ ! -x "$DRIVER_FOAM_BIN" ]]; then
    echo "ERROR: driverFoam not found or not executable at $DRIVER_FOAM_BIN" >&2
    exit 1
fi

"$DRIVER_FOAM_BIN" sweep-run \
    --spec sweep_temporal_convergence.json \
    --output-dir results/sweepCases \
    --max-cases 200

echo "Bidomain Temporal sweep complete."

# Refresh the Canonical CSV for Paper I
python3 "$_PAPERI_ROOT/applications/scripts/paperI_results/aggregate.py" bidomain_temporal \
    --repo-root "$_PAPERI_ROOT"
