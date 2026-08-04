#!/bin/bash
# Run the full conformal tetrahedral bath-bidomain convergence ladder
# (N=10, 20, 40, 80) with bathPredictorCorrector true, so the headline
# bath orders match the production solver default.
#
# Strategy: patch bathPredictorCorrector true in the tet electroProperties
# overlay, then call run_parallel_interface_sweep.sh (same physics as
# baseline) but post-copy results to a _predictor sub-tree.
#
# Results land in:
#   interfaceStudy/matchedSubmesh/distanceWeightedHarmonic_predictor/N*/
#   results/bath_tet_predictor_convergence.csv
#
# The registered baseline (…/distanceWeightedHarmonic/N*/) is not touched
# because the sweep overwrites that path only for the method name it loops
# over; we rename after each solve.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$CASE_DIR"

# ── OpenFOAM environment ────────────────────────────────────────────────────
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
if [[ ! -f "$OPENFOAM_BASHRC" ]]; then
    echo "ERROR: OpenFOAM bashrc not found: $OPENFOAM_BASHRC" >&2
    exit 2
fi
set +eu
source "$OPENFOAM_BASHRC" > /dev/null
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" > /dev/null 2>&1
set -eu

# ── Parameters ──────────────────────────────────────────────────────────────
RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40 80}"
NPROCS="${NPROCS:-6}"
PYTHON="${PYTHON:-python3}"

read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"

OVERLAY="$SCRIPT_DIR/electroProperties"
PREDICTOR_DIR="$SCRIPT_DIR/interfaceStudy/matchedSubmesh/distanceWeightedHarmonic_predictor"
RESULTS_DIR="$SCRIPT_DIR/results"
RESULTS_CSV="$RESULTS_DIR/bath_tet_predictor_convergence.csv"

# ── Helper: timestep schedule (matches run_parallel_interface_sweep.sh) ─────
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

# ── Back up files we will modify ────────────────────────────────────────────
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
    echo "Restored case inputs."
}
trap restore_inputs EXIT

# ── Apply the tet overlays (same as run_parallel_interface_sweep.sh does) ───
cp "$SCRIPT_DIR/electroProperties" constant/electroProperties
cp "$SCRIPT_DIR/fvSchemes"         system/fvSchemes

# ── Set the assembly to the valid matchedSubmesh value ──────────────────────
# ── Patch bathPredictorCorrector true ───────────────────────────────────────
foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPredictorCorrector \
    -set true > /dev/null

echo "=== bath_tet predictor-corrector ladder ==="
echo "bathPredictorCorrector: $(grep bathPredictorCorrector constant/electroProperties)"
echo "Resolutions: $RESOLUTIONS_STR"
echo "Output: $PREDICTOR_DIR"

# ── Run each resolution ──────────────────────────────────────────────────────
for N in "${RESOLUTIONS[@]}"; do
    echo ""
    echo "=== N=$N ==="

    # Re-use the pre-built mesh from the baseline run (identical geometry)
    BANK="$SCRIPT_DIR/interfaceMeshBank/N$N"
    if [[ -f "$BANK/polyMesh.tar.gz" ]]; then
        echo "Using banked mesh for N=$N"
        rm -rf constant/polyMesh
        tar -xzf "$BANK/polyMesh.tar.gz"
    else
        echo "No banked mesh for N=$N — running mesh gate"
        bash "$SCRIPT_DIR/run_mesh_gate.sh" "$N"
        # Bank it for possible reuse
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

    # Opt-in signed per-face dump (FACE_ERRORS=1). Off by default: the files are
    # one row per interface face, which is fine at N<=40 but large at N=80.
    FACE_FLAG=""
    if [[ "${FACE_ERRORS:-0}" == "1" ]]; then FACE_FLAG="-writeFaceErrors"; fi

    bathBidomainInterfaceMetrics -latestTime $FACE_FLAG \
        > "$OUT_DIR/log.interfaceMetrics" 2>&1

    cp postProcessing/bathBidomainInterfaceMetrics.csv \
        "$OUT_DIR/bathBidomainInterfaceMetrics.csv"

    # Exact-field isolation on the same converged case: the manufactured phiE
    # replaces the solved one, so the identical interface-flux construction is
    # evaluated with no solve, coupling-loop or algebraic error in its input.
    # Comparing the two files attributes the assembled-current-density stall to
    # either the interface discretisation or the solve.
    bathBidomainInterfaceMetrics -latestTime -exactFields $FACE_FLAG \
        > "$OUT_DIR/log.interfaceMetricsExactField" 2>&1
    cp postProcessing/bathBidomainInterfaceMetricsExactField.csv \
        "$OUT_DIR/bathBidomainInterfaceMetricsExactField.csv"
    if [[ -n "$FACE_FLAG" ]]; then
        cp postProcessing/bathBidomainFaceErrors*.csv "$OUT_DIR/" 2>/dev/null || true
    fi
    cp postProcessing/bathBidomain_3D_*_cells_implicit.dat \
        "$OUT_DIR/summary.dat"
    cp "$BANK/polyMesh.sha256"  "$OUT_DIR/" 2>/dev/null || true
    cp "$RESULTS_DIR/N$N/log.checkMesh" "$OUT_DIR/" 2>/dev/null || true
    cp "$RESULTS_DIR/N$N/mesh_manifest.txt" "$OUT_DIR/" 2>/dev/null || true

    echo "N=$N done."
done

echo ""
echo "=== All resolutions complete. Summarising… ==="

# ── Build summary CSV ────────────────────────────────────────────────────────
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
        echo "WARNING: output not found for N=$n, skipping in CSV" >&2
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

echo ""
echo "Done. Predictor-corrector convergence CSV:"
echo "  $RESULTS_CSV"
