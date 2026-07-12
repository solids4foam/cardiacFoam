#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../.." && pwd)"
STRUCTURED="$REPO_ROOT/tutorials/manufacturedSolutions/bathBidomain"
WORK="$(mktemp -d)"
cleanup()
{
    if [[ "${KEEP_WORK:-0}" == "1" ]]
    then
        echo "Preserving structured work directory: $WORK"
    else
        rm -rf "$WORK"
    fi
}
trap cleanup EXIT

dt_for_n()
{
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        *) echo "Unsupported N=$1" >&2; exit 2 ;;
    esac
}

steps_for_n()
{
    case "$1" in
        10) echo 2 ;;
        20) echo 9 ;;
        40) echo 36 ;;
        *) echo "Unsupported N=$1" >&2; exit 2 ;;
    esac
}

mkdir -p "$WORK/constant" "$WORK/system"
cp "$CASE_DIR/constant/physicsProperties" "$WORK/constant/"
cp "$CASE_DIR/constant/electroProperties" "$WORK/constant/"
cp "$CASE_DIR/system/controlDict" "$WORK/system/"
cp "$STRUCTURED/system/fvSchemes" "$WORK/system/"
cp "$CASE_DIR/system/fvSolution" "$WORK/system/"
cp "$STRUCTURED/system/topoSetDict" "$WORK/system/"
cp "$STRUCTURED/system/setTorsoOrganConductivityFieldDict" "$WORK/system/"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40}"
read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"
ASSEMBLY="${ASSEMBLY:-currentSplit}"
METHODS_STR="${METHODS:-unweightedHarmonic distanceWeightedHarmonic naiveLinearSigmaTotal}"
read -r -a METHODS <<< "$METHODS_STR"

cd "$WORK"
foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPotentialDomain.intracellularAssembly \
    -set "$ASSEMBLY"
for N in "${RESOLUTIONS[@]}"
do
    echo "=== structured N=$N ==="
    rm -rf constant/polyMesh
    sed -E "s/\(80 80 80\)/($N $N $N)/g" \
        "$STRUCTURED/system/blockMeshDict.3D" \
        > system/blockMeshDict
    blockMesh > "log.blockMesh.N$N" 2>&1
    topoSet > "log.topoSet.N$N" 2>&1
    checkMesh > "log.checkMesh.N$N" 2>&1

    foamDictionary system/controlDict -entry deltaT -set "$(dt_for_n "$N")"
    foamDictionary system/controlDict -entry endTime -set 0.02
    foamDictionary system/controlDict -entry writeControl -set timeStep
    foamDictionary system/controlDict -entry writeInterval -set "$(steps_for_n "$N")"

    for method in "${METHODS[@]}"
    do
        echo "--- $method ---"
        rm -rf 0 postProcessing [0-9]*
        foamDictionary constant/electroProperties \
            -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
            -set "$method"
        setTorsoOrganConductivityField > "log.setConductivity.$method.N$N" 2>&1
        cardiacFoam > "log.cardiacFoam.$method.N$N" 2>&1
        bathBidomainInterfaceMetrics -latestTime \
            > "log.interfaceMetrics.$method.N$N" 2>&1

        if [[ "$ASSEMBLY" == "currentSplit" ]]
        then
            OUT_DIR="$CASE_DIR/setup/structuredInterpolationSweep/$method/N$N"
        else
            OUT_DIR="$CASE_DIR/setup/structuredInterpolationSweep/$ASSEMBLY/$method/N$N"
        fi
        rm -rf "$OUT_DIR"
        mkdir -p "$OUT_DIR"
        cp "log.cardiacFoam.$method.N$N" "$OUT_DIR/"
        cp "log.interfaceMetrics.$method.N$N" "$OUT_DIR/"
        cp "log.checkMesh.N$N" "$OUT_DIR/log.checkMesh"
        cp postProcessing/bathBidomain_3D_*_cells_implicit.dat \
            "$OUT_DIR/summary.dat"
        cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/"
    done
done

echo "Structured interpolation sweep complete."
