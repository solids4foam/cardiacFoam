#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
# Darwin/SIP strips DYLD_LIBRARY_PATH across a fresh bash exec; RunFunctions
# restores it from FOAM_LD_LIBRARY_PATH (see its own "Darwin workaround" block).
# Without this, cardiacFoam aborts with "Library not loaded: @rpath/libOpenFOAM.dylib".
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1
set -eu

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

ELECTRO_BACKUP="$(mktemp)"
CONTROL_BACKUP="$(mktemp)"
SCHEMES_BACKUP="$(mktemp)"
cp constant/electroProperties "$ELECTRO_BACKUP"
cp system/controlDict "$CONTROL_BACKUP"
cp system/fvSchemes "$SCHEMES_BACKUP"
restore_inputs()
{
    cp "$ELECTRO_BACKUP" constant/electroProperties
    cp "$CONTROL_BACKUP" system/controlDict
    cp "$SCHEMES_BACKUP" system/fvSchemes
    rm -f "$ELECTRO_BACKUP" "$CONTROL_BACKUP" "$SCHEMES_BACKUP"
}
trap restore_inputs EXIT

# Merged case (bathBidomain): constant/electroProperties is the hex default
# (dimension "1D", no bathPotentialDomain interface-assembly entries). Activate
# this case's own tet overlay -- dimension "3D" plus the unstructured-interface
# settings the foamDictionary calls below refine -- before running.
cp setup/mesh/tet/electroProperties constant/electroProperties

foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
    -set distanceWeightedHarmonic

# Merged case (bathBidomain): system/fvSchemes is the hex default (Gauss
# linear); activate this case's own tet overlay (leastSquares) as the
# "corrected" baseline -- previously the standalone case's own permanent
# default -- before the "limitedCorrection" branch refines it further.
cp setup/mesh/tet/fvSchemes system/fvSchemes

SCHEME="${SCHEME:-corrected}"
if [[ "$SCHEME" == "limitedCorrection" ]]
then
    foamDictionary system/fvSchemes -entry gradSchemes.default -set leastSquares
    foamDictionary system/fvSchemes -entry laplacianSchemes.default \
        -set "Gauss linear limited 0.5"
    foamDictionary system/fvSchemes -entry snGradSchemes.default -set "limited 0.5"
elif [[ "$SCHEME" != "corrected" ]]
then
    echo "Unknown SCHEME=$SCHEME" >&2
    exit 2
fi

RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40}"
read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"
for N in "${RESOLUTIONS[@]}"
do
    echo "=== matched serial N=$N ==="
    bash setup/mesh/tet/run_mesh_gate.sh "$N"
    foamDictionary system/controlDict -entry deltaT -set "$(dt_for_n "$N")"
    foamDictionary system/controlDict -entry endTime -set 0.02
    foamDictionary system/controlDict -entry writeControl -set timeStep
    foamDictionary system/controlDict -entry writeInterval -set "$(steps_for_n "$N")"

    rm -rf 0 postProcessing [0-9]* processor*
    setTorsoOrganConductivityField > "log.setConductivity.matched.N$N" 2>&1
    cardiacFoam > "log.cardiacFoam.matched.N$N" 2>&1
    bathBidomainInterfaceMetrics -latestTime \
        > "log.interfaceMetrics.matched.N$N" 2>&1

    if [[ "$SCHEME" == "corrected" ]]
    then
        OUT_DIR="setup/mesh/tet/matchedSubmeshStudy/N$N"
    else
        OUT_DIR="setup/mesh/tet/matchedSubmeshLimitedStudy/N$N"
    fi
    rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR"
    cp "log.cardiacFoam.matched.N$N" "$OUT_DIR/log.cardiacFoam"
    cp "log.interfaceMetrics.matched.N$N" "$OUT_DIR/log.interfaceMetrics"
    cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/"
done

restore_inputs
trap - EXIT

echo "Matched serial sweep complete."
