#!/bin/bash
# Is the tetrahedral N=80 bath point an algebraic artefact?
#
# The archived N=80 ladder point is internally inconsistent in a way that
# discretisation error cannot produce: x1AssembledFlux_L2 = 7.48e-1 against a
# physical flux scale of alpha = 1e-2 (an "error" 75x the quantity itself) and
# heartPhiE_L2 40x worse than N=40, while x0PhiE_L2 on the same mesh and the
# same solve keeps halving cleanly across all four levels. The suspected cause
# is system/fvSolution's
#
#     "phiE|phiEFinal|phiI|phiIFinal" { tolerance 1e-15; relTol 0; maxIter 5000; }
#
# PCG exhausting maxIter and exiting silently with whatever residual it has.
#
# Proving that needs no convergence study and no long run. Two short passes:
#
#   Tier 1  Run a handful of timesteps and read the phiE iteration count and
#           final residual straight off the solver log. maxIter reached with a
#           residual far above 1e-15 on the first step is direct evidence.
#
#   Tier 2  Run bathBidomainInterfaceMetrics -exactFields on the same case.
#           The exact-field probe substitutes the manufactured phiE, so it does
#           not depend on the solve at all. If its assembled-flux error lands on
#           the clean ladder while the solved value does not, the mesh and the
#           interface operator are exonerated and only the solve is at fault.
#
# NOTE on Tier 2's scope: from a truncated run only the assembled-flux columns
# are meaningful. The intracellular-leak columns read phiI, which has not
# converged after a few steps, and the potential columns describe a transient
# state. Do not quote anything else from this run.
#
# Usage:
#   ./run_bath_tet_solve_diagnostic.sh            # N=80, 2 steps, maxIter as-is
#   N=40 ./run_bath_tet_solve_diagnostic.sh       # control: a level that works
#   MAXITER=50000 ./run_bath_tet_solve_diagnostic.sh   # confirming pass
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$CASE_DIR"

OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
[[ -f "$OPENFOAM_BASHRC" ]] || { echo "ERROR: no OpenFOAM at $OPENFOAM_BASHRC" >&2; exit 2; }
set +eu
source "$OPENFOAM_BASHRC" > /dev/null
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" > /dev/null 2>&1
set -eu

N="${N:-80}"
STEPS="${STEPS:-2}"
NPROCS="${NPROCS:-6}"
MAXITER="${MAXITER:-}"

dt_for_n() { case "$1" in
    10) echo 0.00892857 ;; 20) echo 0.00224215 ;;
    40) echo 0.000560538 ;; 80) echo 0.0001401345 ;;
    *) echo "Unsupported N=$1" >&2; exit 2 ;;
esac; }

DT="$(dt_for_n "$N")"
END="$(awk -v d="$DT" -v s="$STEPS" 'BEGIN { printf "%.17g", d*s }')"

TAG="N${N}_${STEPS}step${MAXITER:+_maxIter$MAXITER}${PREDICTOR:+_predictor$PREDICTOR}"
OUT_DIR="$SCRIPT_DIR/interfaceStudy/solveDiagnostic/$TAG"
rm -rf "$OUT_DIR"; mkdir -p "$OUT_DIR"

# Restore every dictionary this script rewrites, whatever happens.
ELECTRO_BACKUP="$(mktemp)"; cp constant/electroProperties "$ELECTRO_BACKUP"
CONTROL_BACKUP="$(mktemp)"; cp system/controlDict         "$CONTROL_BACKUP"
SCHEMES_BACKUP="$(mktemp)"; cp system/fvSchemes           "$SCHEMES_BACKUP"
SOLUTION_BACKUP="$(mktemp)"; cp system/fvSolution         "$SOLUTION_BACKUP"
trap 'cp "$ELECTRO_BACKUP" constant/electroProperties; \
      cp "$CONTROL_BACKUP" system/controlDict; \
      cp "$SCHEMES_BACKUP" system/fvSchemes; \
      cp "$SOLUTION_BACKUP" system/fvSolution; \
      rm -f "$ELECTRO_BACKUP" "$CONTROL_BACKUP" "$SCHEMES_BACKUP" "$SOLUTION_BACKUP"' EXIT

# Same overlays and production default as the headline ladder, so the solve
# being diagnosed is the one the paper reports.
cp "$SCRIPT_DIR/electroProperties" constant/electroProperties
cp "$SCRIPT_DIR/fvSchemes"         system/fvSchemes
# PREDICTOR=false runs the baseline one-pass phi_e -> V_m coupling instead of
# the production predictor-corrector, which is how the corrector is tested as a
# candidate mechanism for the N=80 x=0 instability.
foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPredictorCorrector \
    -set "${PREDICTOR:-true}" > /dev/null
foamDictionary constant/electroProperties \
    -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
    -set distanceWeightedHarmonic > /dev/null

if [[ -n "$MAXITER" ]]; then
    foamDictionary system/fvSolution \
        -entry 'solvers."phiE|phiEFinal|phiI|phiIFinal".maxIter' \
        -set "$MAXITER" > /dev/null
fi

BANK="$SCRIPT_DIR/interfaceMeshBank/N$N"
if [[ -f "$BANK/polyMesh.tar.gz" ]]; then
    rm -rf constant/polyMesh
    tar -xzf "$BANK/polyMesh.tar.gz"
else
    echo "No banked mesh for N=$N. Run run_mesh_gate.sh $N first." >&2
    exit 1
fi

foamDictionary system/controlDict -entry deltaT        -set "$DT"    > /dev/null
foamDictionary system/controlDict -entry endTime       -set "$END"   > /dev/null
foamDictionary system/controlDict -entry writeControl  -set timeStep > /dev/null
foamDictionary system/controlDict -entry writeInterval -set "$STEPS" > /dev/null

echo "=== solve diagnostic: N=$N, $STEPS steps, dt=$DT, endTime=$END ==="
echo "maxIter: ${MAXITER:-as configured (5000)}"
echo "bathPredictorCorrector: ${PREDICTOR:-true}"

rm -rf 0 postProcessing [0-9]* processor*
setTorsoOrganConductivityField > "$OUT_DIR/log.setConductivity" 2>&1

if [[ "$NPROCS" -gt 1 ]]; then
    foamDictionary system/decomposeParDict \
        -entry numberOfSubdomains -set "$NPROCS" > /dev/null 2>&1 || true
    decomposePar -force > "$OUT_DIR/log.decomposePar" 2>&1
    mpirun -np "$NPROCS" cardiacFoam -parallel > "$OUT_DIR/log.cardiacFoam" 2>&1
    reconstructPar -latestTime > "$OUT_DIR/log.reconstructPar" 2>&1
else
    cardiacFoam > "$OUT_DIR/log.cardiacFoam" 2>&1
fi

# Guard: -latestTime silently selects time 0 if nothing past 0 was written.
if [[ -z "$(find . -maxdepth 1 -type d -regex '\./0\.[0-9].*' | head -1)" ]]; then
    echo "FAILED: solver wrote no time directory past 0" >&2
    exit 1
fi

bathBidomainInterfaceMetrics -latestTime \
    > "$OUT_DIR/log.interfaceMetrics" 2>&1
cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/"

bathBidomainInterfaceMetrics -latestTime -exactFields \
    > "$OUT_DIR/log.interfaceMetricsExactField" 2>&1
cp postProcessing/bathBidomainInterfaceMetricsExactField.csv "$OUT_DIR/"

# Tier 1 evidence, extracted so it does not have to be re-grepped by hand.
grep -E "Solving for phiE|Solving for phiI" "$OUT_DIR/log.cardiacFoam" \
    > "$OUT_DIR/phiE_solver_performance.txt" || true

echo
echo "=== Tier 1: phiE solver performance ==="
if [[ -s "$OUT_DIR/phiE_solver_performance.txt" ]]; then
    head -12 "$OUT_DIR/phiE_solver_performance.txt"
else
    echo "(no 'Solving for phiE' lines -- check log.cardiacFoam)"
fi
echo
echo "Results in $OUT_DIR"
