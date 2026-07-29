#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

dt_for_n()
{
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        80) echo 0.0001401345 ;;
        *) echo "Unsupported N=$1" >&2; exit 2 ;;
    esac
}

steps_for_n()
{
    case "$1" in
        10) echo 2 ;;
        20) echo 9 ;;
        40) echo 36 ;;
        80) echo 143 ;;
        *) echo "Unsupported N=$1" >&2; exit 2 ;;
    esac
}

ELECTRO_BACKUP="$(mktemp)"
CONTROL_BACKUP="$(mktemp)"
cp constant/electroProperties "$ELECTRO_BACKUP"
cp system/controlDict "$CONTROL_BACKUP"
SCHEMES_BACKUP="$(mktemp)"
cp system/fvSchemes "$SCHEMES_BACKUP"
restore_inputs()
{
    cp "$ELECTRO_BACKUP" constant/electroProperties
    cp "$CONTROL_BACKUP" system/controlDict
    cp "$SCHEMES_BACKUP" system/fvSchemes
    rm -f "$SCHEMES_BACKUP"
    rm -f "$ELECTRO_BACKUP" "$CONTROL_BACKUP"
}
trap restore_inputs EXIT

RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40}"
read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"
NPROCS="${NPROCS:-6}"
ASSEMBLY="${ASSEMBLY:-currentSplit}"
METHODS_STR="${METHODS:-unweightedHarmonic distanceWeightedHarmonic}"
read -r -a METHODS <<< "$METHODS_STR"

# Merged case (bathBidomain): constant/electroProperties is the hex default
# (dimension "1D", no bathPotentialDomain interface-assembly entries).
# Activate this case's own tet overlay (dimension "3D" plus the
# unstructured-interface settings) before refining it below.
cp setup/mesh/tet/electroProperties constant/electroProperties
# The tet overlay must also replace system/fvSchemes. The case-root fvSchemes
# is the hexahedral one and uses Gauss linear cell gradients, which do not
# converge on this tetrahedral family; without this the sweep silently solves
# a non-convergent problem and reports errors flat under refinement.
cp setup/mesh/tet/fvSchemes system/fvSchemes

foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPotentialDomain.intracellularAssembly \
    -set "$ASSEMBLY"

for N in "${RESOLUTIONS[@]}"
do
    echo "=== parallel interface sweep N=$N ==="
    bash setup/mesh/tet/run_mesh_gate.sh "$N"

    BANK="setup/mesh/tet/interfaceMeshBank/N$N"
    rm -rf "$BANK"
    mkdir -p "$BANK"
    tar -czf "$BANK/polyMesh.tar.gz" constant/polyMesh
    find constant/polyMesh -type f -print0 \
        | sort -z \
        | xargs -0 shasum -a 256 \
        > "$BANK/polyMesh.sha256"

    foamDictionary system/controlDict -entry deltaT -set "$(dt_for_n "$N")"
    foamDictionary system/controlDict -entry endTime -set 0.02
    foamDictionary system/controlDict -entry writeControl -set timeStep
    foamDictionary system/controlDict -entry writeInterval -set "$(steps_for_n "$N")"

    for method in "${METHODS[@]}"
    do
        echo "--- $method ---"
        rm -rf 0 postProcessing [0-9]* processor*
        foamDictionary constant/electroProperties \
            -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
            -set "$method"
        setTorsoOrganConductivityField > "log.setConductivity.$method.N$N" 2>&1
        decomposePar -force > "log.decomposePar.$method.N$N" 2>&1
        mpirun -np "$NPROCS" cardiacFoam -parallel \
            > "log.cardiacFoam.$method.N$N" 2>&1
        reconstructPar -latestTime > "log.reconstructPar.$method.N$N" 2>&1
        bathBidomainInterfaceMetrics -latestTime \
            > "log.interfaceMetrics.$method.N$N" 2>&1
        cp postProcessing/bathBidomainInterfaceMetrics.csv \
            "bathBidomainInterfaceMetrics.numerical.$method.N$N.csv"

        if [[ "$ASSEMBLY" == "currentSplit" ]]
        then
            OUT_DIR="setup/mesh/tet/interfaceStudy/$method/N$N"
        else
            OUT_DIR="setup/mesh/tet/interfaceStudy/$ASSEMBLY/$method/N$N"
        fi
        rm -rf "$OUT_DIR"
        mkdir -p "$OUT_DIR"
        cp "log.cardiacFoam.$method.N$N" "$OUT_DIR/log.cardiacFoam"
        cp "log.interfaceMetrics.$method.N$N" "$OUT_DIR/log.interfaceMetrics"
        cp "log.reconstructPar.$method.N$N" "$OUT_DIR/log.reconstructPar"
        cp setup/mesh/tet/results/N$N/log.checkMesh "$OUT_DIR/"
        cp postProcessing/bathBidomain_3D_*_cells_implicit.dat \
            "$OUT_DIR/summary.dat"
        cp "bathBidomainInterfaceMetrics.numerical.$method.N$N.csv" \
            "$OUT_DIR/bathBidomainInterfaceMetrics.csv"
        cp "$BANK/polyMesh.sha256" "$OUT_DIR/"
    done
done

restore_inputs
trap - EXIT

echo "Parallel N=10/20/40 interface sweep complete."
