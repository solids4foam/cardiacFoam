#!/usr/bin/env bash
# Collect all canonical Paper I result CSVs (+ Niederer cached points) into one
# staging tree for Zenodo upload. Read-only over the repo; writes only the bundle.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO_ROOT/paperI_data_bundle"
rm -rf "$OUT"; mkdir -p "$OUT/manufactured" "$OUT/niederer"

find "$REPO_ROOT/tutorials/manufacturedSolutions" \
     -path "*/setup/results/*_convergence.csv" -print0 \
  | while IFS= read -r -d '' f; do
        cp "$f" "$OUT/manufactured/$(basename "$f")"
    done

NIED="$REPO_ROOT/tutorials/NiedererEtAl2011/NiedererEtAl2011verification/setup/cachedCasePostProcessing"
find "$NIED" -name "*_points_DT0001_DX*.csv" -exec cp {} "$OUT/niederer/" \; 2>/dev/null || true

echo "Bundle staged at: $OUT"
find "$OUT" -type f | sort
