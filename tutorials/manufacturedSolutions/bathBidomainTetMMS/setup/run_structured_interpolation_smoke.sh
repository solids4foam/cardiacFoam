#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../.." && pwd)"
STRUCTURED="$REPO_ROOT/tutorials/manufacturedSolutions/bathBidomain"
WORK="$(mktemp -d)"
cleanup()
{
    rm -rf "$WORK"
}
trap cleanup EXIT

mkdir -p "$WORK/constant" "$WORK/system"
cp "$CASE_DIR/constant/physicsProperties" "$WORK/constant/"
cp "$CASE_DIR/constant/electroProperties" "$WORK/constant/"
cp "$CASE_DIR/system/controlDict" "$WORK/system/"
cp "$STRUCTURED/system/fvSchemes" "$WORK/system/"
cp "$CASE_DIR/system/fvSolution" "$WORK/system/"
cp "$STRUCTURED/system/topoSetDict" "$WORK/system/"
cp "$STRUCTURED/system/setTorsoOrganConductivityFieldDict" "$WORK/system/"
sed -E 's/\(80 80 80\)/(10 10 10)/g' \
    "$STRUCTURED/system/blockMeshDict.3D" \
    > "$WORK/system/blockMeshDict"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

cd "$WORK"
blockMesh > log.blockMesh 2>&1
topoSet > log.topoSet 2>&1

for method in \
    unweightedHarmonic \
    distanceWeightedHarmonic \
    naiveLinearSigmaTotal
do
    echo "=== structured interpolation smoke: $method ==="
    rm -rf 0 postProcessing [0-9]*
    foamDictionary constant/electroProperties \
        -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
        -set "$method"
    setTorsoOrganConductivityField > "log.setConductivity.$method" 2>&1
    cardiacFoam > "log.cardiacFoam.$method" 2>&1
    bathBidomainInterfaceMetrics -latestTime \
        > "log.interfaceMetrics.$method" 2>&1

    OUT_DIR="$CASE_DIR/setup/structuredInterpolationSmoke/$method"
    rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR"
    cp "log.cardiacFoam.$method" "log.interfaceMetrics.$method" "$OUT_DIR/"
    cp postProcessing/bathBidomain_3D_*_cells_implicit.dat "$OUT_DIR/summary.dat"
    cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/"
done

echo "Structured N=10 interpolation smoke matrix complete."
