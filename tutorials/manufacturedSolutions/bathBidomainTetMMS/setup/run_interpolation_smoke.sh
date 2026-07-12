#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

bash setup/run_mesh_gate.sh 10

ELECTRO_BACKUP="$(mktemp)"
cp constant/electroProperties "$ELECTRO_BACKUP"
restore_electro_properties()
{
    cp "$ELECTRO_BACKUP" constant/electroProperties
    rm -f "$ELECTRO_BACKUP"
}
trap restore_electro_properties EXIT

for method in \
    unweightedHarmonic \
    distanceWeightedHarmonic \
    naiveLinearSigmaTotal
do
    echo "=== interpolation smoke: $method ==="
    rm -rf 0 postProcessing [0-9]*
    foamDictionary constant/electroProperties \
        -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
        -set "$method"

    setTorsoOrganConductivityField > "log.setConductivity.$method" 2>&1
    cardiacFoam > "log.cardiacFoam.$method" 2>&1

    if ! grep -q '^End$' "log.cardiacFoam.$method"
    then
        echo "cardiacFoam did not complete for $method" >&2
        exit 1
    fi

    bathBidomainInterfaceMetrics -latestTime \
        > "log.interfaceMetrics.$method" 2>&1

    OUT_DIR="setup/interpolationSmoke/$method"
    rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR"
    cp \
        "log.cardiacFoam.$method" \
        "log.setConductivity.$method" \
        "log.interfaceMetrics.$method" \
        "$OUT_DIR/"
    cp postProcessing/bathBidomain_3D_*_cells_implicit.dat "$OUT_DIR/summary.dat"
    cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/"
done

restore_electro_properties
trap - EXIT

echo "Interpolation selector smoke matrix complete."
