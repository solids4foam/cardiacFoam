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

N="${N:-10}"
NPROCS="${NPROCS:-6}"
METHOD="${METHOD:-distanceWeightedHarmonic}"
ASSEMBLY="${ASSEMBLY:-currentSplit}"

case "$N" in
    10) DELTA_T=0.00892857; STEPS=2 ;;
    20) DELTA_T=0.00224215; STEPS=9 ;;
    40) DELTA_T=0.000560538; STEPS=36 ;;
    *) echo "Unsupported N=$N" >&2; exit 2 ;;
esac

ELECTRO_BACKUP="$(mktemp)"
CONTROL_BACKUP="$(mktemp)"
cp constant/electroProperties "$ELECTRO_BACKUP"
cp system/controlDict "$CONTROL_BACKUP"
restore_inputs()
{
    cp "$ELECTRO_BACKUP" constant/electroProperties
    cp "$CONTROL_BACKUP" system/controlDict
    rm -f "$ELECTRO_BACKUP" "$CONTROL_BACKUP"
}
trap restore_inputs EXIT

# Merged case (bathBidomain): constant/electroProperties is the hex default
# (dimension "1D", no bathPotentialDomain interface-assembly entries).
# Activate this case's own tet overlay (dimension "3D" plus the
# unstructured-interface settings) before the foamDictionary calls below.
cp setup/mesh/tet/electroProperties constant/electroProperties

bash setup/mesh/tet/run_mesh_gate.sh "$N"
foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
    -set "$METHOD"
foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPotentialDomain.intracellularAssembly \
    -set "$ASSEMBLY"
foamDictionary system/controlDict -entry deltaT -set "$DELTA_T"
foamDictionary system/controlDict -entry endTime -set 0.02
foamDictionary system/controlDict -entry writeControl -set timeStep
foamDictionary system/controlDict -entry writeInterval -set "$STEPS"

OUT_DIR="setup/mesh/tet/parallelEquivalence/$ASSEMBLY/$METHOD/N$N"
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

# Serial and parallel solves deliberately reuse the same constant/polyMesh.
rm -rf 0 postProcessing [0-9]* processor*
setTorsoOrganConductivityField > "$OUT_DIR/log.setConductivity.serial" 2>&1
cardiacFoam > "$OUT_DIR/log.cardiacFoam.serial" 2>&1
bathBidomainInterfaceMetrics -latestTime \
    > "$OUT_DIR/log.interfaceMetrics.serial" 2>&1
cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/serial.csv"

rm -rf 0 postProcessing [0-9]* processor*
setTorsoOrganConductivityField > "$OUT_DIR/log.setConductivity.parallel" 2>&1
decomposePar -force > "$OUT_DIR/log.decomposePar" 2>&1
mpirun -np "$NPROCS" cardiacFoam -parallel \
    > "$OUT_DIR/log.cardiacFoam.parallel" 2>&1
reconstructPar -latestTime > "$OUT_DIR/log.reconstructPar" 2>&1
bathBidomainInterfaceMetrics -latestTime \
    > "$OUT_DIR/log.interfaceMetrics.parallel" 2>&1
cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/parallel.csv"

python3 setup/mesh/tet/compare_metrics.py \
    "$OUT_DIR/serial.csv" \
    "$OUT_DIR/parallel.csv" \
    --out "$OUT_DIR/comparison.csv"

restore_inputs
trap - EXIT

echo "Same-mesh serial/parallel equivalence complete for N=$N, $METHOD."
