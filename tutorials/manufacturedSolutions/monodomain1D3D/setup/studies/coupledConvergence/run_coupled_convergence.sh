#!/bin/bash
# ---------------------------------------------------------------------------
# Coupled 1D-3D Convergence Sweeps
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

_PAPERI_ROOT="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
DRIVER_FOAM_BIN="$_PAPERI_ROOT/applications/scripts/driverFoam/bin/driverFoam"

if [[ ! -x "$DRIVER_FOAM_BIN" ]]; then
    echo "ERROR: driverFoam not found or not executable at $DRIVER_FOAM_BIN" >&2
    exit 1
fi

echo "Running decoupled sweep..."
"$DRIVER_FOAM_BIN" sweep-run --spec sweep_decoupled.json --max-cases 200

echo "Running active (unidirectional) sweep..."
"$DRIVER_FOAM_BIN" sweep-run --spec sweep_active.json --max-cases 200

echo "Running bidirectional sweep..."
"$DRIVER_FOAM_BIN" sweep-run --spec sweep_bidirectional.json --max-cases 200

echo "Coupled 1D-3D sweeps complete."

# Refresh the Canonical CSV for Paper I
python3 "$_PAPERI_ROOT/applications/scripts/paperI_results/aggregate.py" purkinje_monodomain_coupled \
    --repo-root "$_PAPERI_ROOT"
