#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/.venv/bin/python3
cd "$CASE_DIR"

source /Volumes/OpenFOAM-v2412/etc/bashrc

# N -> dt, reused verbatim from manufactured_fda.py's DT_VALUES so results
# are directly comparable to the orthogonal-mesh baseline sweep. A plain
# function (not an associative array) so this runs under macOS's stock
# bash 3.2, which has no associative-array support.
dt_for_n() {
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        80) echo 0.000140174 ;;
        *) echo "no dt configured for N=$1" >&2; exit 1 ;;
    esac
}

# AMPLITUDES arrives as a plain space-separated string (from the env var
# set in Task 7), not a bash array -- deliberately unquoted below so it
# word-splits into separate tokens.
AMPLITUDES="${AMPLITUDES:-0.0 0.05 0.10 0.15 0.20}"
RESOLUTIONS=(10 20 40 80)

for A in $AMPLITUDES; do
    OUT_DIR="setup/results/$A"
    mkdir -p "$OUT_DIR"
    for N in "${RESOLUTIONS[@]}"; do
        echo "=== A=$A N=$N ==="
        ./Allclean

        sed -E "s/\([0-9]+ [0-9]+ [0-9]+\) simpleGrading/($N $N $N) simpleGrading/" \
            system/blockMeshDict.3D > /tmp/blockMeshDict.3D.tmp
        mv /tmp/blockMeshDict.3D.tmp system/blockMeshDict.3D

        DT="$(dt_for_n "$N")"
        sed -E "s/^deltaT.*/deltaT    $DT;/" system/controlDict > /tmp/controlDict.tmp
        mv /tmp/controlDict.tmp system/controlDict

        blockMesh -dict system/blockMeshDict.3D > log.blockMesh 2>&1

        if [[ "$A" != "0.0" ]]; then
            "$PY" setup/distort_mesh.py . -N "$N" -A "$A" --dimension 3
        fi

        checkMesh > log.checkMesh 2>&1 || echo "checkMesh reported issues at A=$A N=$N (see log.checkMesh)"
        cp log.checkMesh "$OUT_DIR/log.checkMesh.$N"

        ./Allrun parallel  # decomposePar + runParallel cardiacFoam + reconstructPar, via RunFunctions

        cp "postProcessing/3D_${N}_cells_implicit.dat" "$OUT_DIR/"

        git checkout -- system/blockMeshDict.3D system/controlDict
    done
done

./Allclean
echo "Sweep complete. Summarizing..."
"$PY" setup/summarize_results.py setup/results \
    --amplitudes $AMPLITUDES \
    --resolutions "${RESOLUTIONS[@]}" \
    --out setup/results/summary.csv
