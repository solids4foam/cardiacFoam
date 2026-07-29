#!/bin/bash
# Full manufactured-FDA monodomain + pseudo-ECG sweep via the driverFoam engine.
# Produces driverPostProcessingArchive_postProcessing/*.dat for mono_spatial + pseudo_ecg_spatial.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"
"$DRIVER" all --entry manufacturedFDA --config "$SCRIPT_DIR/studies/spatialConvergence/config.json"
echo "monodomainPseudoECG manufacturedFDA sweep complete."
