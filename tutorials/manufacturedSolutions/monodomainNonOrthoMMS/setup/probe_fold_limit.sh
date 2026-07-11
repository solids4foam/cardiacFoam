#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/.venv/bin/python3
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

for A in 0.10 0.15 0.20 0.25 0.30 0.35 0.40; do
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
        echo "A=$A: OK, max non-orthogonality = ${max_ortho:-unknown}"
    else
        echo "A=$A: checkMesh FAILED -> STOP, use previous A"
        break
    fi
    git checkout -- system/blockMeshDict.3D
done

./Allclean
