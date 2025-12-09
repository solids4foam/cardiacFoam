#!/bin/bash

# Usage: ./run_cases.sh <caseDir> <1D|2D|3D>

CASE_DIR="$1"
DIM="$2"

if [[ -z "$CASE_DIR" || -z "$DIM" ]]; then
    echo "Usage: ./run_cases.sh <caseDir> <1D|2D|3D>"
    exit 1
fi

cd "$CASE_DIR"

echo "🧹 Cleaning case"
./Allclean

# Source OpenFOAM environment
source /Volumes/OpenFOAM-v2412/etc/bashrc

DICT="system/blockMeshDict.$DIM"
echo "📐 Running blockMesh with: $DICT (output → log.blockMesh)"

# Run blockMesh silently, write all output to log.blockMesh
blockMesh -dict "$DICT" > log.blockMesh 2>&1
if [[ $? -ne 0 ]]; then
    echo "❌ blockMesh failed — see log.blockMesh"
    exit 1
fi

echo "🚀 Running case in parallel"
./Allrun parallel

echo "✍️ Extracting results to output folder"


