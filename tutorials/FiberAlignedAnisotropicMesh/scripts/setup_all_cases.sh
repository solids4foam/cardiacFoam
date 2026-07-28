#!/bin/bash
set -euo pipefail
source /Volumes/OpenFOAM-v2412/etc/bashrc
export PATH=/Users/simaocastro/cobivecco-OpenFOam/platforms/darwin64ClangDPInt32Opt/bin:$PATH

# Canonical physical scale for BOTH pipelines (native ~4x5x4 units -> ~8x10x8 cm).
SCALE="(0.02 0.02 0.02)"

# Resolve the tutorial root from this script's location, independent of CWD.
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
TUT_DIR="$( cd "$SCRIPT_DIR/.." &> /dev/null && pwd )"
cd "$TUT_DIR"

for res in coarse medium fine; do
    echo "Setting up monodomain_${res}..."
    rm -rf "simulations/monodomain_${res}"
    cp -r templates/monodomainHeartTissue "simulations/monodomain_${res}"
    cd "simulations/monodomain_${res}"

    # Ensure UVC divSchemes are present in fvSchemes.
    python3 -c '
with open("system/fvSchemes", "r") as f: data = f.read()
if "div(phiTrajectoryDistance)" not in data:
    data = data.replace("div(phiU,psi)   Gauss linear;", "div(phiU,psi)   Gauss linear;\n    div(phiTrajectoryDistance) Gauss upwind;\n    div(phi,tmDistEpi)  Gauss upwind;\n    div(phi,tmDistEndo) Gauss upwind;\n    div(phi,abDistApex) Gauss upwind;\n    div(phi,abDistBase) Gauss upwind;")
    with open("system/fvSchemes", "w") as f: f.write(data)
'

    # Ensure UVC solvers are present in fvSolution.
    python3 -c '
with open("system/fvSolution", "r") as f: data = f.read()
if "tvLaplace" not in data:
    data = data.replace("solvers\n{", "solvers\n{\n    \"(tv|tvLaplace.*|tm.*|ab.*|rt.*|ridgeLaplace.*|dSeptPost|dSeptAnt|dFreePost|dFreeAnt|psiAbLvGuide|psiOtLvGuide|wLvGuide|psiAbRvGuide|psiOtRvGuide|wRvGuide)\"\n    {\n        solver          GAMG;\n        tolerance       1e-10;\n        relTol          0;\n        smoother        GaussSeidel;\n    }\n")
    with open("system/fvSolution", "w") as f: f.write(data)
'

    rm -rf constant/polyMesh 0 log.* results*
    gmshToFoam "../../meshes/gmsh/heart_aniso_${res}.msh" > log.gmshToFoam 2>&1
    polyDualMesh 75 -concaveMultiCells > log.polyDualMesh 2>&1

    # Generate content on the NATIVE mesh (scale-invariant; see README).
    fiberFoamComputeHeartAxes    > log.fiberFoam 2>&1
    fiberFoamComputeCoordinates >> log.fiberFoam 2>&1
    fiberFoamBuildFrames        >> log.fiberFoam 2>&1
    fiberFoamGenerateFibers     >> log.fiberFoam 2>&1
    setCardiacConductivity      >> log.fiberFoam 2>&1

    # Scale-last: single transformPoints as the final geometry step.
    transformPoints -scale "$SCALE" > log.transformPoints 2>&1

    # Keep only fields needed by the solver.
    echo "Cleaning intermediate fields for monodomain_${res}..."
    find 0 -maxdepth 1 -type f | grep -vE "fiber|sheet|normalDirection|Conductivity|sigma|cellToRegion|tags|uvc_transmural" | xargs rm -f

    touch results.foam

    if [ "$res" == "medium" ]; then
        sed -i '' 's/nNonOrthogonalCorrectors.*/nNonOrthogonalCorrectors 2;/g' system/fvSolution
    fi

    cd "$TUT_DIR"
done
echo "All cases setup successfully!"
