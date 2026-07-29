#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

N="${N:-20}"
NPROCS="${NPROCS:-6}"
SCHEMES_BACKUP="$(mktemp)"
cp system/fvSchemes "$SCHEMES_BACKUP"
restore_schemes()
{
    cp "$SCHEMES_BACKUP" system/fvSchemes
    rm -f "$SCHEMES_BACKUP"
}
trap restore_schemes EXIT

run_variant()
{
    local name="$1"
    local grad="$2"
    local laplacian="$3"
    local sn_grad="$4"

    cp "$SCHEMES_BACKUP" system/fvSchemes
    foamDictionary system/fvSchemes -entry gradSchemes.default -set "$grad"
    foamDictionary system/fvSchemes -entry laplacianSchemes.default -set "$laplacian"
    foamDictionary system/fvSchemes -entry snGradSchemes.default -set "$sn_grad"

    echo "=== $name: grad=[$grad], laplacian=[$laplacian], snGrad=[$sn_grad] ==="
    RESOLUTIONS="$N" METHODS=distanceWeightedHarmonic NPROCS="$NPROCS" \
        bash studies/run_parallel_interface_sweep.sh

    local source="studies/mesh/tet/interfaceStudy/distanceWeightedHarmonic/N$N"
    local target="studies/mesh/tet/studies/gradientScheme/$name/N$N"
    rm -rf "$target"
    mkdir -p "$target"
    cp -R "$source"/. "$target"/
    cp system/fvSchemes "$target/fvSchemes"
}

VARIANTS_STR="${VARIANTS:-current gaussLinear limitedCorrection orthogonalControl}"
for variant in $VARIANTS_STR
do
    case "$variant" in
        current)
            run_variant current leastSquares "Gauss linear corrected" corrected
            ;;
        gaussLinear)
            run_variant gaussLinear "Gauss linear" "Gauss linear corrected" corrected
            ;;
        limitedCorrection)
            run_variant limitedCorrection leastSquares \
                "Gauss linear limited 0.5" "limited 0.5"
            ;;
        orthogonalControl)
            run_variant orthogonalControl leastSquares \
                "Gauss linear orthogonal" orthogonal
            ;;
        *) echo "Unknown gradient variant '$variant'" >&2; exit 2 ;;
    esac
done

restore_schemes
trap - EXIT

echo "Gradient-scheme screen complete at N=$N."
