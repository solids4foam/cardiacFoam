#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

for dim in 1D 2D 3D
do
    echo
    echo "Running manufactured eikonal ECG ${dim}"
    "$SCRIPT_DIR/run_cases.sh" "$CASE_DIR" "$dim"
done

# --- Paper I: persist canonical convergence CSV (additive; does not alter the sweep above) ---
_PAPERI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
python3 "$_PAPERI_ROOT/applications/scripts/paperI_results/aggregate.py" eikonal \
    --repo-root "$_PAPERI_ROOT" \
    || echo "WARN: paperI aggregate (eikonal) failed; native output untouched" >&2
