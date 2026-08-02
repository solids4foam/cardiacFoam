#!/bin/bash
# run_eikonal_tet_optimised.sh
#
# Runs the eikonal MMS ladder on the anatomically-optimised tetrahedral family
# (Frontal algorithm, Netgen optimisation, Smoothing=100) with least-squares
# gradient reconstruction only.
#
# Gauss-linear is intentionally omitted: its failure on non-orthogonal tets is
# already established on the generic Delaunay family and is not expected to
# change with improved bulk quality. The scientific question here is whether
# leastSquares recovers toward second order when the mean non-orthogonality
# matches the anatomical target (~14-16 deg vs ~18.6 deg for the generic family).
#
# Results are written to setup/results/scheme_study_optimised.csv so the
# original generic-family CSV is preserved alongside for direct comparison.
#
# Usage:
#   bash run_eikonal_tet_optimised.sh
#   RESOLUTIONS="10 20" bash run_eikonal_tet_optimised.sh   # smoke test

set +e
if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
  if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
    source /Volumes/OpenFOAM-v2412/etc/bashrc > /dev/null 2>&1
  else
    echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
    exit 2
  fi
fi
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" > /dev/null 2>&1

PY="${PYTHON:-python3}"
RESOLUTIONS="${RESOLUTIONS:-10 20 40 80}"

# The optimised template lives alongside the original
TPL_OPT="setup/mesh/tet/box.geo.template.optimised"
TPL_ORIG="setup/mesh/tet/box.geo.template"
OUT_CSV="setup/results/scheme_study_optimised.csv"

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
OUTER_ITERATION_CAP="$(awk '/nOuterCorrectors/{gsub(/;/, "", $2); print $2; exit}' setup/mesh/tet/fvSolution)"

# Sanity-check that the optimised template exists
if [[ ! -f "$TPL_OPT" ]]; then
  echo "ERROR: $TPL_OPT not found." >&2
  echo "Run characterize_optimised.py first — it creates the template as a side-effect." >&2
  exit 1
fi

# Activate tet fvSolution overlay and least-squares gradient scheme; restore on exit
_FVSOL_BAK="$(mktemp)"; cp system/fvSolution "$_FVSOL_BAK"
_FVSCH_BAK="$(mktemp)"; cp system/fvSchemes "$_FVSCH_BAK"
trap 'cp "$_FVSOL_BAK" system/fvSolution; cp "$_FVSCH_BAK" system/fvSchemes; rm -f "$_FVSOL_BAK" "$_FVSCH_BAK"' EXIT

cp setup/mesh/tet/fvSolution system/fvSolution

# Force least-squares gradient scheme
sed -E 's/^([[:space:]]*)default([[:space:]]+)Gauss linear;[[:space:]]*$/\1default\2leastSquares;/' \
    system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes

# Phase 1 was run separately. Proceeding directly to Phase 2 (Solver).



# ── Metric extraction helpers (identical to run_eikonal_tet.sh) ──────────────
activation_metrics(){
  awk '/Number of cells/{n=$NF} /^activationTime /{l2=$3;li=$4} \
       END{dx=1.0/int((n)^(1/3)+0.5); printf "%s %s %s",dx,l2,li}' \
       postProcessing/manufacturedEikonalActivationTime.dat
}
ecg_metrics(){
  awk 'NR>=8{if($3+0>a)a=$3; if($4+0>b)b=$4} END{printf "%s %s",a,b}' \
       postProcessing/manufacturedEikonalECGSummary.dat
}

# ── Single-level runner ───────────────────────────────────────────────────────
run_one(){
  local TAG="$1" N="$2"
  ./Allclean > /dev/null 2>&1

  # Instantiate the OPTIMISED template (not the original)
  LC=$($PY -c "print(1.0/$N)")
  sed "s|__LC__|$LC|" "$TPL_OPT" > setup/mesh/tet/box.geo

  gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 > /dev/null 2>&1
  gmshToFoam box.msh > /dev/null 2>&1
  rm -f box.msh

  decomposePar > log.cf 2>&1
  mpirun --oversubscribe -np 6 cardiacFoam -parallel >> log.cf 2>&1
  RC=$?
  reconstructPar >> log.cf 2>&1

  mkdir -p setup/results/logs
  cp log.cf "setup/results/logs/${TAG}_optimised_N${N}.log"

  ITERS="$(grep -c '^PIMPLE: iteration' log.cf || true)"
  if [[ -n "$OUTER_ITERATION_CAP" && "$ITERS" -ge "$OUTER_ITERATION_CAP" ]]; then
    echo "WARNING: $TAG N=$N hit the nOuterCorrectors cap ($OUTER_ITERATION_CAP) -- check logs" >&2
  fi

  if [[ $RC -ne 0 || ! -f postProcessing/manufacturedEikonalActivationTime.dat ]]; then
    echo "SKIPPING $TAG N=$N: cardiacFoam exit=$RC or no verifier output" >&2
    return
  fi

  M=$(activation_metrics); E=$(ecg_metrics)
  echo "$TAG,$N,${M// /,},${E// /,},$ITERS" >> "$OUT_CSV"
  echo "done $TAG N=$N : activationTime[$M] ecg[$E] outerIters=$ITERS"
}

# ── Main sweep: leastSquares only ────────────────────────────────────────────
echo "scheme,N,dx,activationTime_L2,activationTime_Linf,ecg_L2,ecg_Linf,outerIterations" > "$OUT_CSV"

for N in $RESOLUTIONS; do
  run_one leastSquares_optimised "$N"
done

# Restore tutorial default gradient scheme (leastSquares is already set, but be explicit)
sed -E 's/^([[:space:]]*)default([[:space:]]+)leastSquares.*$/\1default\2leastSquares;/' \
    system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes

./Allclean > /dev/null 2>&1

echo ""
echo "=== scheme_study_optimised.csv ==="
cat "$OUT_CSV"
echo ""
echo "Compare with generic Delaunay results in: setup/results/scheme_study.csv"
