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
# Merged case (monodomainPseudoECG): system/ holds the shared dicts + the hex
# fvSolution/fvSchemes/controlDict defaults. This sweep activates the tet
# fvSolution overlay and repeatedly overwrites fvSchemes (gradScheme A/B) and
# controlDict (per-N deltaT/endTime) in place -- back up all three and restore
# byte-for-byte on exit so the shared hex sweeps (mono_spatial,
# pseudo_ecg_spatial) never see a leaked tet-run setting.
_FVSOL_BAK="$(mktemp)"; cp system/fvSolution "$_FVSOL_BAK"
cp setup/mesh/tet/fvSolution system/fvSolution
_FVSCH_BAK="$(mktemp)"; cp system/fvSchemes "$_FVSCH_BAK"
_CD_BAK="$(mktemp)"; cp system/controlDict "$_CD_BAK"
trap 'cp "$_FVSOL_BAK" system/fvSolution; rm -f "$_FVSOL_BAK";
      cp "$_FVSCH_BAK" system/fvSchemes; rm -f "$_FVSCH_BAK";
      cp "$_CD_BAK" system/controlDict; rm -f "$_CD_BAK"' EXIT
dt_for_n(){ case "$1" in 10) echo 0.00892857;;20) echo 0.00224215;;40) echo 0.000560538;;80) echo 0.000140174;;esac; }
set_grad(){
  # system/fvSchemes' gradSchemes.default line carries a trailing
  # "//leastSquares;" documentation comment -- match past it (as the
  # else-branch already did) so the substitution isn't a silent no-op.
  if [ "$1" = "leastSquares" ]; then
    sed -E 's/^([[:space:]]*)default([[:space:]]+)Gauss linear;.*$/\1default\2leastSquares;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  else
    sed -E 's/^([[:space:]]*)default([[:space:]]+)leastSquares.*$/\1default\2Gauss linear;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  fi
}
mono_metrics(){ D=$(ls postProcessing/3D_*_cells_implicit.dat 2>/dev/null|head -1); awk '/Grid spacing/{dx=$NF} /^Vm /{l2=$3;li=$4} END{printf "%s %s %s",dx,l2,li}' "$D"; }
# Fail loudly rather than emitting empty ECG columns. This script runs the
# committed case directory, whose constant/electroProperties must therefore
# configure ecgDomains; driverFOAM runs are unaffected because the driver
# generates its own configuration from the manufacturedFDA entry defaults.
ecg_metrics(){
  local f=postProcessing/manufacturedPseudoECGSummary.dat
  if [[ ! -f "$f" ]]; then
    echo "ERROR: $f not written. The case's constant/electroProperties has no" >&2
    echo "       ecgDomains block, so no pseudo-ECG was computed. Either add it" >&2
    echo "       or drive this study through driverFOAM (entry manufacturedFDA)," >&2
    echo "       whose defaults already carry the five electrodes." >&2
    exit 1
  fi
  awk 'NR>=7{if($3+0>a)a=$3; if($4+0>b)b=$4} END{printf "%s %s",a,b}' "$f"
}
run_one(){ # $1 tag $2 N
  ./Allclean >/dev/null 2>&1
  LC=$($PY -c "print(1.0/$2)"); sed "s|__LC__|$LC|" setup/mesh/tet/box.geo.template > setup/mesh/tet/box.geo
  gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1; gmshToFoam box.msh >/dev/null 2>&1; rm -f box.msh
  DT=$(dt_for_n $2); sed -E "s/^deltaT.*/deltaT    $DT;/; s/^endTime.*/endTime    0.2;/" system/controlDict>/tmp/cd && mv /tmp/cd system/controlDict
  # Parallel (same solver, verified serial/parallel-equivalent elsewhere in this
  # repo): decomposePar + runParallel is markedly faster than serial on N=80.
  runApplication decomposePar
  runParallel cardiacFoam
  runApplication reconstructPar
  M=$(mono_metrics); E=$(ecg_metrics)
  echo "$1,$2,${M// /,},${E// /,}" >> setup/results/scheme_study.csv
  echo "done $1 N=$2 : mono[$M] ecg[$E]"
}
echo "scheme,N,dx,mono_L2,mono_Linf,ecg_L2,ecg_Linf" > setup/results/scheme_study.csv
set_grad "GaussLinear"; for N in $RESOLUTIONS; do run_one GaussLinear $N; done
set_grad "leastSquares"; for N in $RESOLUTIONS; do run_one leastSquares $N; done
./Allclean >/dev/null 2>&1
echo "=== scheme_study.csv ==="; cat setup/results/scheme_study.csv

# --- Paper I: persist canonical convergence CSV (additive; does not alter the sweep above) ---
_PAPERI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../../.." && pwd)"
python3 "$_PAPERI_ROOT/applications/scripts/paperI_results/aggregate.py" mono_tet \
    --repo-root "$_PAPERI_ROOT" \
    || echo "WARN: paperI aggregate (tet) failed; native output untouched" >&2
