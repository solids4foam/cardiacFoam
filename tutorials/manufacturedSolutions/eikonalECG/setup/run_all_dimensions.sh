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
