#!/bin/bash
set -euo pipefail

# Unstructured (tetrahedral) Niederer slab dx-convergence sweep.
#
# Companion to the structured-hex verification case. Geometry, conductivity,
# stimulus, ionic model, fvSchemes and activation probes are reused verbatim
# from ../NiedererEtAl2011verification. The mesh is built with gmsh (tets at
# characteristic length lc = dx) instead of blockMesh, and the diffusion solve
# is switched to implicit + nNonOrthogonalCorrectors 2 (see README) so the
# non-orthogonal correction converges and the tet result stays second order.
#
# The timestep is LOCKED at the coarsest dt of the published Niederer sweep
# (0.01 ms = 1e-5 s) for every dx, so this is a pure spatial-convergence
# study: the corner activation times should converge toward the community
# reference (NiedererEtAl2012.reference) as dx -> 0. This is NOT a pass/fail
# regression -- the coarse tet meshes are expected to miss the hex-tuned
# tolerances; the point is the convergence trend, not a single-mesh match.

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$CASE_DIR"

# Source OpenFOAM if it is not already in the environment.
if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    set +eu
    source /Volumes/OpenFOAM-v2412/etc/bashrc
    set -eu
fi

# Locked timestep (coarsest dt of the Niederer sweep, 0.01 ms).
DELTAT="${DELTAT:-1e-5}"

# endTime per dx, taken from the driverFOAM niederer2012 defaults
# (openfoam_driver/core/defaults/niederer_2012.py: END_TIME_BY_DX). The coarser
# the mesh the slower the numerical conduction velocity, so the wave needs
# longer to activate the far corner: dx=0.5 -> 0.2 s, dx=0.2 -> 0.08 s,
# dx=0.1 -> 0.055 s. This is what lets the opposite corner activate at all.
end_time_for_dx() {
    case "$1" in
        0.5) echo 0.2   ;;
        0.2) echo 0.08  ;;
        0.1) echo 0.055 ;;
        *) echo "no endTime configured for dx=$1 mm" >&2; exit 1 ;;
    esac
}

# dx ladder in millimetres, mirroring the published Niederer dx sweep.
# Overridable from the environment, e.g. DX_VALUES="0.5 0.2" for a quick pass.
DX_VALUES_STR="${DX_VALUES:-0.5 0.2 0.1}"
read -r -a DX_VALUES <<< "$DX_VALUES_STR"

RESULTS_DIR="setup/results"
mkdir -p "$RESULTS_DIR"

for DX in "${DX_VALUES[@]}"; do
    # lc in metres = dx[mm] * 1e-3
    LC="$(awk -v d="$DX" 'BEGIN{printf "%g", d*1e-3}')"
    ENDTIME="$(end_time_for_dx "$DX")"
    echo "=== tet Niederer  dx=${DX} mm  (lc=${LC} m, dt=${DELTAT} s, endTime=${ENDTIME} s) ==="

    OUT_DIR="$RESULTS_DIR/dx_${DX}mm"
    mkdir -p "$OUT_DIR"
    ./Allclean > /dev/null 2>&1 || true

    # 1) mesh at this characteristic length
    sed "s/__LC__/$LC/" setup/slab.geo.template > setup/slab.geo
    gmsh -3 setup/slab.geo -o slab.msh -format msh2 > log.gmsh 2>&1
    gmshToFoam slab.msh > log.gmshToFoam 2>&1
    rm -f slab.msh

    # 2) lock time controls (dt fixed across all dx)
    sed -E "s/^deltaT.*/deltaT    $DELTAT;/"   system/controlDict > /tmp/nied_cd.tmp
    sed -E "s/^endTime.*/endTime    $ENDTIME;/" /tmp/nied_cd.tmp   > system/controlDict
    rm -f /tmp/nied_cd.tmp

    # 3) record mesh quality
    checkMesh > log.checkMesh 2>&1 \
        || echo "checkMesh reported issues at dx=${DX} (see log.checkMesh)"
    cp log.checkMesh "$OUT_DIR/log.checkMesh"

    # 4) solve (parallel: decomposePar -> runParallel -> reconstructPar) and
    #    extract the activation probes/lines. Parallelism is set by
    #    system/decomposeParDict (numberOfSubdomains). Set RUN_SERIAL=1 to
    #    force a serial run instead.
    if [[ "${RUN_SERIAL:-0}" == "1" ]]; then
        ./Allrun
    else
        ./Allrun parallel
    fi

    # 5) collect the activation-time output for this dx
    if [[ -f postProcessing/Niedererpoints/0/activationTime ]]; then
        cp postProcessing/Niedererpoints/0/activationTime "$OUT_DIR/points_activationTime.dat"
    else
        echo "WARNING: no Niedererpoints activationTime found for dx=${DX}" >&2
    fi
    if [[ -f postProcessing/Niedererlines/0/activationTime ]]; then
        cp postProcessing/Niedererlines/0/activationTime "$OUT_DIR/lines_activationTime.dat"
    fi

    # restore the committed controlDict
    git checkout -- system/controlDict 2>/dev/null || true
done

./Allclean > /dev/null 2>&1 || true
echo
echo "Sweep complete. Per-dx activation times under $RESULTS_DIR/dx_*mm/"
echo "Compare against NiedererEtAl2012.reference (convergence, not pass/fail)."
