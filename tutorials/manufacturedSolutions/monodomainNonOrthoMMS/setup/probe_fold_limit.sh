#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/.venv/bin/python3
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

# Real-mesh calibration target from tutorials/PATHOS/LBBB and RBBB checkMesh
# (identical geometry, both give: max non-orthogonality 69.6317, average
# 15.1477, max skewness 1.17752). Stop as soon as we meet/exceed that, or
# the mesh folds (negative volume / checkMesh failure), whichever comes first.
TARGET_NON_ORTHO=69.6317

for A in 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.60 0.80 1.00 1.20 1.50 1.80 2.20; do
    ./Allclean
    sed -i '' 's/([0-9]* [0-9]* [0-9]*) simpleGrading/(10 10 10) simpleGrading/' system/blockMeshDict.3D
    blockMesh -dict system/blockMeshDict.3D > log.blockMesh 2>&1
    "$PY" setup/distort_mesh.py . -N 10 -A "$A" --dimension 3
    if checkMesh > log.checkMesh 2>&1; then
        if grep -q "cells with negative volumes\|has negative volume" log.checkMesh; then
            echo "A=$A: checkMesh ran but flagged negative volumes -> STOP, use previous A"
            break
        fi
        max_ortho=$(grep -o "non-orthogonality Max: [0-9.]*" log.checkMesh | grep -o "[0-9.]*$")
        max_skew=$(grep -o "Max skewness = [0-9.]*" log.checkMesh | grep -o "[0-9.]*$")
        echo "A=$A: OK, max non-orthogonality = ${max_ortho:-unknown}, max skewness = ${max_skew:-unknown}"
        if [[ -n "${max_ortho:-}" ]] && (( $(echo "$max_ortho >= $TARGET_NON_ORTHO" | bc -l) )); then
            echo "A=$A: reached/exceeded real-mesh non-orthogonality benchmark ($TARGET_NON_ORTHO) -> STOP, use this A"
            break
        fi
    else
        echo "A=$A: checkMesh FAILED -> STOP, use previous A"
        break
    fi
    git checkout -- system/blockMeshDict.3D
done

git checkout -- system/blockMeshDict.3D 2>/dev/null || true
./Allclean
