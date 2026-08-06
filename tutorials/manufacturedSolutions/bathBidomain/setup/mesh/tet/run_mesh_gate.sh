#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
N="${1:-10}"

if ! [[ "$N" =~ ^[1-9][0-9]*$ ]]
then
    echo "N must be a positive integer" >&2
    exit 2
fi

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

cd "$CASE_DIR"
OUT_DIR="setup/mesh/tet/results/N${N}"
rm -rf constant/polyMesh "$OUT_DIR"
mkdir -p "$OUT_DIR"

# Merged case (bathBidomain): system/controlDict's writeFormat is the hex
# default (binary); gmshToFoam writes constant/polyMesh in that same format,
# but verify_mesh.py (and checkMesh's own text parsing) need ASCII. Force it
# here -- the caller's own controlDict backup/restore already covers this
# field wholesale, so no separate backup is needed.
foamDictionary system/controlDict -entry writeFormat -set ascii >/dev/null

LC="$(awk -v n="$N" 'BEGIN { printf "%.17g", 1.0/n }')"
sed "s/__LC__/${LC}/" setup/mesh/tet/three_domain_box.geo.template \
    > setup/mesh/tet/three_domain_box.geo

gmsh -3 setup/mesh/tet/three_domain_box.geo -o three_domain_box.msh -format msh2 \
    > log.gmsh 2>&1
gmshToFoam three_domain_box.msh > log.gmshToFoam 2>&1
rm -f three_domain_box.msh

checkMesh > log.checkMesh 2>&1
# Keep the stricter diagnostic as an archived audit. OpenFOAM's extended cell
# determinant test flags edge-adjacent tetrahedra with sparse cell-centre
# stencils even when the standard topology and geometry checks pass.
checkMesh -allTopology -allGeometry > log.checkMesh.strict 2>&1 || true
python3 setup/mesh/tet/verify_mesh.py "$CASE_DIR" > "$OUT_DIR/mesh_manifest.txt"

cp log.gmsh log.gmshToFoam log.checkMesh log.checkMesh.strict "$OUT_DIR/"
cp setup/mesh/tet/three_domain_box.geo "$OUT_DIR/"

{
    echo "nominal_N=$N"
    echo "gmsh_lc=$LC"
    echo "gmsh_version=$(gmsh --version 2>&1 | tail -1)"
    echo "openfoam_version=${WM_PROJECT_VERSION:-unknown}"
    echo "git_sha=$(git -C "$CASE_DIR" rev-parse HEAD)"
    echo "strict_small_determinant_cells=$(sed -n 's/.*number of cells: \([0-9][0-9]*\).*/\1/p' log.checkMesh.strict | tail -1)"
} >> "$OUT_DIR/mesh_manifest.txt"

cat "$OUT_DIR/mesh_manifest.txt"
echo "Mesh gate passed. Review $OUT_DIR/log.checkMesh before solving."
