#!/bin/bash
set -euo pipefail

# Tetrahedral-mesh variant of the monodomain manufactured MMS sweep.
#
# Reuses monodomainPseudoECG's 3D exact solution, conductivity, and ionic
# forcing (identical constant/ and system/ dicts); only the mesh generator
# changes: gmsh builds a near-uniform tetrahedral
# mesh of the unit cube at characteristic length lc = 1/N, converted with
# gmshToFoam. This exercises the Gauss linear corrected Laplacian on a
# genuinely unstructured, high-non-orthogonality mesh -- the case a
# perturbed hex mesh cannot reach without folding.

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/.venv/bin/python3
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

# dt per nominal resolution N, reused verbatim from the orthogonal /
# non-orthogonal hex sweeps so the temporal error stays subdominant and the
# tet results are comparable to those baselines. N here is the *nominal*
# cells-per-side; the tet count is ~5-6x N^3 and the verifier reports the
# effective spacing dx = 1/cbrt(nCells).
dt_for_n() {
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        80) echo 0.000140174 ;;
        *) echo "no dt configured for N=$1" >&2; exit 1 ;;
    esac
}

# RESOLUTIONS and ENDTIME are overridable from the environment. ENDTIME
# defaults to a short window (temporal error is O(dt^2) with the backward
# ddt scheme, so a few steps suffice to keep tet runs tractable); set
# ENDTIME=0.2 to match the hex baseline exactly for the paper table.
RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40}"
ENDTIME="${ENDTIME:-0.02}"
read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"

for N in "${RESOLUTIONS[@]}"; do
    echo "=== tet N=$N (lc=1/$N) ==="
    OUT_DIR="setup/results/$N"
    mkdir -p "$OUT_DIR"
    ./Allclean

    # 1) build the gmsh geometry at this characteristic length
    LC="$("$PY" -c "print(1.0/$N)")"
    sed "s/__LC__/$LC/" setup/box.geo.template > setup/box.geo

    # 2) mesh with gmsh (legacy msh2) and import into OpenFOAM
    gmsh -3 setup/box.geo -o box.msh -format msh2 > log.gmsh 2>&1
    gmshToFoam box.msh > log.gmshToFoam 2>&1
    rm -f box.msh

    # 3) time controls for this resolution
    DT="$(dt_for_n "$N")"
    sed -E "s/^deltaT.*/deltaT    $DT;/" system/controlDict > /tmp/tet_controlDict.tmp
    sed -E "s/^endTime.*/endTime    $ENDTIME;/" /tmp/tet_controlDict.tmp > system/controlDict
    rm -f /tmp/tet_controlDict.tmp

    # 4) record mesh quality
    checkMesh > log.checkMesh 2>&1 \
        || echo "checkMesh reported issues at N=$N (see log.checkMesh)"
    cp log.checkMesh "$OUT_DIR/log.checkMesh"

    # 5) solve (serial: the verifier .dat lands in the top-level
    #    postProcessing/ without the decomposed-run processor0 quirk)
    ./Allrun

    # 6) collect the verifier summary. The filename embeds cbrt(nCells),
    #    not N, so glob for it.
    DAT="$(ls postProcessing/3D_*_cells_implicit.dat 2>/dev/null | head -1 || true)"
    if [[ -n "$DAT" ]]; then
        cp "$DAT" "$OUT_DIR/summary.dat"
    else
        echo "WARNING: no verifier .dat found for N=$N" >&2
    fi

    # pseudo-ECG manufactured summary (the ECG operator is gradient-based, so
    # it is a second, more grad-scheme-sensitive convergence check per mesh)
    if [[ -f postProcessing/manufacturedPseudoECGSummary.dat ]]; then
        cp postProcessing/manufacturedPseudoECGSummary.dat "$OUT_DIR/pseudoECG_summary.dat"
    fi

    git checkout -- system/controlDict 2>/dev/null || true
done

./Allclean
echo "Sweep complete. Summarizing..."
"$PY" setup/summarize_tet.py setup/results \
    --resolutions "${RESOLUTIONS[@]}" \
    --out setup/results/summary.csv
