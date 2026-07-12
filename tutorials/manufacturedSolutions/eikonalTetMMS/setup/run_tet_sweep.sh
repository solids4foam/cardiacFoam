#!/bin/bash
set -euo pipefail

# Tetrahedral-mesh variant of the eikonal manufactured MMS sweep.
#
# Mirrors monodomainTetMMS/setup/run_tet_sweep.sh: gmsh builds a
# near-uniform tetrahedral mesh of the unit cube at characteristic length
# lc = 1/N (same box.geo.template as the monodomain tet sweep), imported
# with gmshToFoam. Unlike the monodomain case, the eikonal solve is a
# single steady advection-diffusion solve (controlDict is startTime=0,
# endTime=1, deltaT=1 -- one step), so there is no per-N dt to pick and no
# temporal-error/spatial-error coupling to worry about.

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/.venv/bin/python3
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

# RESOLUTIONS is overridable from the environment. N here is the *nominal*
# cells-per-side; the tet count is ~5-6x N^3 and the verifier reports
# "Number of cells"; the summarizer computes the effective spacing
# dx = 1/round(cbrt(nCells)).
RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40}"
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

    # 3) record mesh quality
    checkMesh > log.checkMesh 2>&1 \
        || echo "checkMesh reported issues at N=$N (see log.checkMesh)"
    cp log.checkMesh "$OUT_DIR/log.checkMesh"

    # 4) solve (serial: the verifier .dat lands in the top-level
    #    postProcessing/ without the decomposed-run processor0 quirk)
    ./Allrun

    # archive the solver log so the outer-loop (PIMPLE/Picard) convergence
    # can be checked after the fact -- RunFunctions' runApplication writes
    # log.cardiacFoam in the case dir.
    if [[ -f log.cardiacFoam ]]; then
        cp log.cardiacFoam "$OUT_DIR/log.cardiacFoam"
        LAST_ITER="$(grep -c '^PIMPLE: iteration' log.cardiacFoam || true)"
        FINAL_RESID="$(grep 'Solving for activationTime' log.cardiacFoam | tail -1)"
        echo "  outer iterations: $LAST_ITER, final: $FINAL_RESID"
        if [[ "$LAST_ITER" -ge 2500 ]]; then
            echo "WARNING: N=$N hit the nOuterCorrectors cap (2500) -- outer loop may not have converged" >&2
        fi
    fi

    # 5) collect the activation-time verifier summary
    if [[ -f postProcessing/manufacturedEikonalActivationTime.dat ]]; then
        cp postProcessing/manufacturedEikonalActivationTime.dat "$OUT_DIR/summary.dat"
    else
        echo "WARNING: no activation-time verifier .dat found for N=$N" >&2
    fi

    # 6) pseudo-ECG manufactured summary (gradient-based: this is the
    #    grad-scheme-sensitive check, since gradVm is built analytically
    #    from fvc::grad(activationTime) computed once per solve)
    if [[ -f postProcessing/manufacturedEikonalECGSummary.dat ]]; then
        cp postProcessing/manufacturedEikonalECGSummary.dat "$OUT_DIR/pseudoECG_summary.dat"
    fi
done

./Allclean
echo "Sweep complete. Summarizing..."
"$PY" setup/summarize_tet.py setup/results \
    --resolutions "${RESOLUTIONS[@]}" \
    --out setup/results/summary.csv
