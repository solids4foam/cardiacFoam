#!/bin/bash
set -euo pipefail

cd /Users/simaocastro/noFrontendCardiacFoam_minor_errors/tutorials/manufacturedSolutions/bathBidomain

# 1. Create the optimised template
cat > setup/mesh/tet/three_domain_box.geo.template.optimised << 'EOF'
// Conformal tetrahedral bath--myocardium--bath mesh.
// __LC__ is replaced by setup/run_mesh_gate_optimised.sh.

SetFactory("OpenCASCADE");

lc = __LC__;
eps = 1e-7;

Box(1) = {-1, 0, 0, 1, 1, 1};
Box(2) = { 0, 0, 0, 1, 1, 1};
Box(3) = { 1, 0, 0, 1, 1, 1};

BooleanFragments{ Volume{1, 2, 3}; Delete; }{}

leftBath[] = Volume In BoundingBox{-1-eps, -eps, -eps, 0+eps, 1+eps, 1+eps};
heart[] = Volume In BoundingBox{0-eps, -eps, -eps, 1+eps, 1+eps, 1+eps};
rightBath[] = Volume In BoundingBox{1-eps, -eps, -eps, 2+eps, 1+eps, 1+eps};

Physical Volume("myocardium") = {heart[]};
Physical Volume("bath") = {leftBath[], rightBath[]};

xMin[] = Surface In BoundingBox{-1-eps, -eps, -eps, -1+eps, 1+eps, 1+eps};
xMax[] = Surface In BoundingBox{2-eps, -eps, -eps, 2+eps, 1+eps, 1+eps};
yMin[] = Surface In BoundingBox{-1-eps, -eps, -eps, 2+eps, eps, 1+eps};
yMax[] = Surface In BoundingBox{-1-eps, 1-eps, -eps, 2+eps, 1+eps, 1+eps};
zMin[] = Surface In BoundingBox{-1-eps, -eps, -eps, 2+eps, 1+eps, eps};
zMax[] = Surface In BoundingBox{-1-eps, -eps, 1-eps, 2+eps, 1+eps, 1+eps};

Physical Surface("xMin") = {xMin[]};
Physical Surface("xMax") = {xMax[]};
Physical Surface("sides") = {yMin[], yMax[], zMin[], zMax[]};

Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;
Mesh.Algorithm3D = 4;
Mesh.Optimize = 1;
Mesh.MshFileVersion = 2.2;
Mesh.OptimizeNetgen = 1;
Mesh.Smoothing = 100;
EOF

# 2. Create the optimised mesh gate
cat > setup/mesh/tet/run_mesh_gate_optimised.sh << 'EOF'
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
OUT_DIR="setup/mesh/tet/results/N${N}_optimised"
rm -rf constant/polyMesh "$OUT_DIR"
mkdir -p "$OUT_DIR"

foamDictionary system/controlDict -entry writeFormat -set ascii >/dev/null

LC="$(awk -v n="$N" 'BEGIN { printf "%.17g", 1.0/n }')"
sed "s/__LC__/${LC}/" setup/mesh/tet/three_domain_box.geo.template.optimised \
    > setup/mesh/tet/three_domain_box_optimised.geo

gmsh -3 setup/mesh/tet/three_domain_box_optimised.geo -o three_domain_box_optimised.msh -format msh2 \
    > log.gmsh 2>&1
gmshToFoam three_domain_box_optimised.msh > log.gmshToFoam 2>&1
rm -f three_domain_box_optimised.msh

checkMesh > log.checkMesh 2>&1
checkMesh -allTopology -allGeometry > log.checkMesh.strict 2>&1 || true
python3 setup/mesh/tet/verify_mesh.py "$CASE_DIR" > "$OUT_DIR/mesh_manifest.txt"

cp log.gmsh log.gmshToFoam log.checkMesh log.checkMesh.strict "$OUT_DIR/"
cp setup/mesh/tet/three_domain_box_optimised.geo "$OUT_DIR/"

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
EOF
chmod +x setup/mesh/tet/run_mesh_gate_optimised.sh

# 3. Create the optimised predictor script
cat > setup/mesh/tet/run_bath_tet_predictor_optimised.sh << 'EOF'
#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$CASE_DIR"

OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
set +eu
source "$OPENFOAM_BASHRC" > /dev/null
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" > /dev/null 2>&1
set -eu

RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40 80}"
NPROCS="${NPROCS:-6}"
PYTHON="${PYTHON:-python3}"

read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"

OVERLAY="$SCRIPT_DIR/electroProperties"
PREDICTOR_DIR="$SCRIPT_DIR/interfaceStudy/matchedSubmesh/distanceWeightedHarmonic_predictor_optimised"
RESULTS_DIR="$SCRIPT_DIR/results"
RESULTS_CSV="$RESULTS_DIR/bath_tet_predictor_optimised_convergence.csv"

dt_for_n() {
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        80) echo 0.0001401345 ;;
        *) echo "Unsupported N=$1" >&2; exit 2 ;;
    esac
}
steps_for_n() {
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
SCHEMES_BACKUP="$(mktemp)"
cp constant/electroProperties "$ELECTRO_BACKUP"
cp system/controlDict         "$CONTROL_BACKUP"
cp system/fvSchemes           "$SCHEMES_BACKUP"

restore_inputs() {
    cp "$ELECTRO_BACKUP" constant/electroProperties
    cp "$CONTROL_BACKUP" system/controlDict
    cp "$SCHEMES_BACKUP" system/fvSchemes
    rm -f "$ELECTRO_BACKUP" "$CONTROL_BACKUP" "$SCHEMES_BACKUP"
}
trap restore_inputs EXIT

cp "$SCRIPT_DIR/electroProperties" constant/electroProperties
cp "$SCRIPT_DIR/fvSchemes"         system/fvSchemes

foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPredictorCorrector \
    -set true > /dev/null

for N in "${RESOLUTIONS[@]}"; do
    echo "=== N=$N ==="
    BANK="$SCRIPT_DIR/interfaceMeshBankOptimised/N$N"
    if [[ -f "$BANK/polyMesh.tar.gz" ]]; then
        rm -rf constant/polyMesh
        tar -xzf "$BANK/polyMesh.tar.gz"
    else
        bash "$SCRIPT_DIR/run_mesh_gate_optimised.sh" "$N"
        mkdir -p "$BANK"
        tar -czf "$BANK/polyMesh.tar.gz" constant/polyMesh
        find constant/polyMesh -type f -print0 \
            | sort -z \
            | xargs -0 shasum -a 256 > "$BANK/polyMesh.sha256"
    fi

    foamDictionary system/controlDict -entry deltaT        -set "$(dt_for_n "$N")"  > /dev/null
    foamDictionary system/controlDict -entry endTime       -set 0.02                > /dev/null
    foamDictionary system/controlDict -entry writeControl  -set timeStep            > /dev/null
    foamDictionary system/controlDict -entry writeInterval -set "$(steps_for_n "$N")" > /dev/null

    foamDictionary constant/electroProperties \
        -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
        -set distanceWeightedHarmonic > /dev/null

    OUT_DIR="$PREDICTOR_DIR/N$N"
    rm -rf "$OUT_DIR"
    mkdir -p "$OUT_DIR"
    rm -rf 0 postProcessing [0-9]* processor*

    setTorsoOrganConductivityField > "$OUT_DIR/log.setConductivity" 2>&1
    if [[ "$NPROCS" -gt 1 ]]; then
        foamDictionary system/decomposeParDict \
            -entry numberOfSubdomains -set "$NPROCS" > /dev/null 2>&1 || true
        decomposePar -force > "$OUT_DIR/log.decomposePar" 2>&1
        mpirun -np "$NPROCS" cardiacFoam -parallel \
            > "$OUT_DIR/log.cardiacFoam" 2>&1
        reconstructPar -latestTime > "$OUT_DIR/log.reconstructPar" 2>&1
    else
        cardiacFoam > "$OUT_DIR/log.cardiacFoam" 2>&1
    fi

    bathBidomainInterfaceMetrics -latestTime \
        > "$OUT_DIR/log.interfaceMetrics" 2>&1

    cp postProcessing/bathBidomainInterfaceMetrics.csv \
        "$OUT_DIR/bathBidomainInterfaceMetrics.csv"
    cp postProcessing/bathBidomain_3D_*_cells_implicit.dat \
        "$OUT_DIR/summary.dat"
    cp "$BANK/polyMesh.sha256"  "$OUT_DIR/" 2>/dev/null || true
    cp "$RESULTS_DIR/N${N}_optimised/log.checkMesh" "$OUT_DIR/" 2>/dev/null || true
    cp "$RESULTS_DIR/N${N}_optimised/mesh_manifest.txt" "$OUT_DIR/" 2>/dev/null || true
done

SUMMARY_STAGING="$(mktemp -d)"
cleanup_staging() { rm -rf "$SUMMARY_STAGING"; }

restore_inputs_and_staging() {
    restore_inputs
    cleanup_staging
}
trap restore_inputs_and_staging EXIT

for n in "${RESOLUTIONS[@]}"; do
    SRC="$PREDICTOR_DIR/N$n"
    DST="$SUMMARY_STAGING/N$n"
    if [[ ! -d "$SRC" ]]; then
        continue
    fi
    mkdir -p "$DST"
    ln -sf "$SRC/mesh_manifest.txt" "$DST/mesh_manifest.txt"
    ln -sf "$SRC/log.checkMesh"     "$DST/log.checkMesh"
    ln -sf "$SRC/log.cardiacFoam"   "$DST/log.cardiacFoam"
    ln -sf "$SRC/summary.dat"       "$DST/summary.dat"
done

"$PYTHON" "$SCRIPT_DIR/summarize_tet.py" \
    "$SUMMARY_STAGING" \
    --resolutions "${RESOLUTIONS[@]}" \
    --out "$RESULTS_CSV"

EOF
chmod +x setup/mesh/tet/run_bath_tet_predictor_optimised.sh

# Let's run it in the background to avoid stalling the UI for 5 minutes!
nohup bash setup/mesh/tet/run_bath_tet_predictor_optimised.sh > log.optimised_run 2>&1 &
echo $! > run.pid
