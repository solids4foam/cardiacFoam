#!/bin/bash
set +e
if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
  if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
    source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
  else
    echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
    exit 2
  fi
fi
# Darwin/SIP strips DYLD_LIBRARY_PATH across a fresh bash exec; RunFunctions
# restores it from FOAM_LD_LIBRARY_PATH (see its own "Darwin workaround" block).
# Without this, cardiacFoam aborts with "Library not loaded: @rpath/libOpenFOAM.dylib".
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1
PY="${PYTHON:-python3}"
# leastSquares sweep resolutions, overridable for a faster smoke run. Default
# keeps N=80 so the committed reference (which has N=80 rows) stays reproducible.
RESOLUTIONS="${RESOLUTIONS:-10 20 40 80}"
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
# Merged case (eikonalECG): system/ holds the shared dicts + the hex fvSolution
# default; activate the tet fvSolution overlay for this run and restore on exit.
_FVSOL_BAK="$(mktemp)"; cp system/fvSolution "$_FVSOL_BAK"
cp setup/mesh/tet/fvSolution system/fvSolution
# set_grad below rewrites system/fvSchemes in place. Restore it from the trap
# rather than only on normal completion: an aborted sweep would otherwise leave
# the case root on Gauss linear, and the Cartesian eikonal row would then run
# the wrong gradient scheme without any indication.
_FVSCH_BAK="$(mktemp)"; cp system/fvSchemes "$_FVSCH_BAK"
trap 'cp "$_FVSOL_BAK" system/fvSolution; cp "$_FVSCH_BAK" system/fvSchemes; rm -f "$_FVSOL_BAK" "$_FVSCH_BAK"' EXIT

# gradScheme A/B on the tet mesh: GaussLinear vs leastSquares for
# gradSchemes.default, which is what fvc::grad(activationTime) in
# eikonalECG.C uses to build the single upstream gradTau that gradVm is
# then built from analytically (gradVm = -dU/ds * gradTau). This isolates
# whether that one grad(activationTime) call is itself scheme-sensitive on
# a high-non-orthogonality tet mesh, the same question already answered
# for the monodomain path in ../../../../monodomainPseudoECG/setup/mesh/tet/run_scheme_study.sh.

set_grad(){
  if [ "$1" = "leastSquares" ]; then
    sed -E 's/^([[:space:]]*)default([[:space:]]+)Gauss linear;[[:space:]]*$/\1default\2leastSquares;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  else
    sed -E 's/^([[:space:]]*)default([[:space:]]+)leastSquares.*$/\1default\2Gauss linear;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  fi
}
activation_metrics(){ awk '/Number of cells/{n=$NF} /^activationTime /{l2=$3;li=$4} END{dx=1.0/int((n)^(1/3)+0.5); printf "%s %s %s",dx,l2,li}' postProcessing/manufacturedEikonalActivationTime.dat; }
ecg_metrics(){ awk 'NR>=8{if($3+0>a)a=$3; if($4+0>b)b=$4} END{printf "%s %s",a,b}' postProcessing/manufacturedEikonalECGSummary.dat; }
run_one(){ # $1 tag $2 N
  ./Allclean >/dev/null 2>&1
  LC=$($PY -c "print(1.0/$2)"); sed "s|__LC__|$LC|" setup/mesh/tet/box.geo.template > setup/mesh/tet/box.geo
  gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1; gmshToFoam box.msh >/dev/null 2>&1; rm -f box.msh
  # Parallel (same solver, verified serial/parallel-equivalent elsewhere in this
  # repo): decomposePar + mpirun is markedly faster than serial on N=80.
  decomposePar > log.cf 2>&1
  mpirun --oversubscribe -np 6 cardiacFoam -parallel >> log.cf 2>&1
  RC=$?
  reconstructPar >> log.cf 2>&1
  mkdir -p setup/results/logs
  cp log.cf "setup/results/logs/${1}_N${2}.log"
  ITERS="$(grep -c '^PIMPLE: iteration' log.cf || true)"
  if [[ "$ITERS" -ge 2500 ]]; then
    echo "WARNING: $1 N=$2 hit the nOuterCorrectors cap (2500) -- outer loop may not have converged, see setup/results/logs/${1}_N${2}.log" >&2
  fi
  if [[ $RC -ne 0 || ! -f postProcessing/manufacturedEikonalActivationTime.dat ]]; then
    echo "SKIPPING $1 N=$2: cardiacFoam exit=$RC or no verifier output (see setup/results/logs/${1}_N${2}.log)" >&2
    return
  fi
  M=$(activation_metrics); E=$(ecg_metrics)
  echo "$1,$2,${M// /,},${E// /,},${ITERS}" >> setup/results/scheme_study.csv
  echo "done $1 N=$2 : activationTime[$M] ecg[$E] outerIters=$ITERS"
}
echo "scheme,N,dx,activationTime_L2,activationTime_Linf,ecg_L2,ecg_Linf,outerIterations" > setup/results/scheme_study.csv
set_grad "GaussLinear"; for N in $RESOLUTIONS; do run_one GaussLinear $N; done
set_grad "leastSquares"; for N in $RESOLUTIONS; do run_one leastSquares $N; done
set_grad "leastSquares"   # restore tutorial default
./Allclean >/dev/null 2>&1
echo "=== scheme_study.csv ==="; cat setup/results/scheme_study.csv
