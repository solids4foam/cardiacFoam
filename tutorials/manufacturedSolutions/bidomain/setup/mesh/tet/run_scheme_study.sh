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
# Gradient schemes to sweep, and the phiE/phiI linear-solver tolerance to use.
# Both default to the committed study settings, so an unqualified invocation
# reproduces the archived scheme_study.csv exactly. PHI_TOL=1e-15 runs against
# the case's own committed tolerance instead of the loosened study value, which
# is how the reported maximum-norm behaviour is tested for tolerance sensitivity.
SCHEMES="${SCHEMES:-GaussLinear leastSquares}"
PHI_TOL="${PHI_TOL:-1e-6}"
# Output CSV, overridable so a partial sweep cannot clobber the full study.
OUT="${OUT:-setup/results/scheme_study.csv}"
export PHI_TOL
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
# system/ holds the shared dicts + the hex fvSchemes/controlDict/fvSolution
# defaults; this sweep repeatedly overwrites fvSchemes (gradScheme A/B),
# controlDict (per-N deltaT), and fvSolution (phiE/phiI tolerance, below) in
# place -- back all three up and restore byte-for-byte on exit so the case's
# own hex convergence sweep (committed at tolerance 1e-15) never sees a
# leaked tet-run setting.
# nOuterCorrectors is left at the case's own default (2): the corrector-loop
# sensitivity study (setup/mesh/tet/studies/corrector) already established
# that two outer sweeps are within 1% of the fully converged block, so the
# reported spatial order here is not iteration-error-limited.
_FVSCH_BAK="$(mktemp)"; cp system/fvSchemes "$_FVSCH_BAK"
_CD_BAK="$(mktemp)"; cp system/controlDict "$_CD_BAK"
_FVSOL_BAK="$(mktemp)"; cp system/fvSolution "$_FVSOL_BAK"
trap 'cp "$_FVSCH_BAK" system/fvSchemes; rm -f "$_FVSCH_BAK";
      cp "$_CD_BAK" system/controlDict; rm -f "$_CD_BAK";
      cp "$_FVSOL_BAK" system/fvSolution; rm -f "$_FVSOL_BAK"' EXIT
# The committed phiE/phiI tolerance (1e-15, relTol 0) is tighter than this
# study needs: an isolated N=20 A/B test found identical reported Vm/phiE
# errors to 6 significant figures at tolerance 1e-6 (only ~15-105 PCG
# iterations/step vs 243-250 at 1e-15), confirming the extra iterations
# beyond ~1e-6 just chase the residual toward machine precision without
# changing the manufactured-solution error this study reports. Loosening it
# here cuts the dominant per-timestep cost (the phiE elliptic solve) without
# affecting what gets reported.
"$PY" - <<'PY'
import os, re
from pathlib import Path
tol = os.environ.get("PHI_TOL", "1e-6")
path = Path("system/fvSolution")
text = path.read_text()
pattern = r'("phiE\|phiEFinal\|phiI\|phiIFinal"\s*\{[^}]*?tolerance\s+)1e-15(\s*;)'
n = len(re.findall(pattern, text))
if n != 1:
    raise SystemExit(f"expected 1 committed 1e-15 tolerance block, found {n}")
if tol == "1e-15":
    print("phiE/phiI tolerance: committed 1e-15 retained (no rewrite)")
else:
    path.write_text(re.sub(pattern, r"\g<1>" + tol + r"\g<2>", text))
    print(f"phiE/phiI tolerance: rewritten 1e-15 -> {tol} for this sweep")
PY
dt_for_n(){ case "$1" in 10) echo 0.00892857;;20) echo 0.00224215;;40) echo 0.000560538;;80) echo 0.000140174;;esac; }
set_grad(){
  # system/fvSchemes' gradSchemes.default line carries a trailing
  # "//leastSquares;" documentation comment -- match past it so the
  # substitution isn't a silent no-op.
  if [ "$1" = "leastSquares" ]; then
    sed -E 's/^([[:space:]]*)default([[:space:]]+)Gauss linear;.*$/\1default\2leastSquares;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  else
    sed -E 's/^([[:space:]]*)default([[:space:]]+)leastSquares.*$/\1default\2Gauss linear;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  fi
}
bidomain_metrics(){ D=$(ls postProcessing/3D_*_cells_implicit.dat 2>/dev/null|head -1); awk '/Grid spacing/{dx=$NF} /^Vm /{vl2=$3;vli=$4} /^phiE_gauge /{pl2=$3;pli=$4} END{printf "%s %s %s %s %s",dx,vl2,vli,pl2,pli}' "$D"; }
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
  M=$(bidomain_metrics)
  echo "$1,$2,${M// /,}" >> "$OUT"
  echo "done $1 N=$2 : bidomain[$M]"
}
mkdir -p "$(dirname "$OUT")"
echo "scheme,N,dx,vm_L2,vm_Linf,phiE_L2,phiE_Linf" > "$OUT"
for S in $SCHEMES; do
  set_grad "$S"
  for N in $RESOLUTIONS; do run_one "$S" "$N"; done
done
./Allclean >/dev/null 2>&1
echo "=== $OUT ==="; cat "$OUT"

# --- Paper I: persist canonical convergence CSV (additive; does not alter the sweep above) ---
# Only aggregate when this was the full committed study. A partial or
# alternative-tolerance sweep writes elsewhere and must not overwrite the
# canonical Paper I CSV with a subset of the ladder.
if [ "$OUT" = "setup/results/scheme_study.csv" ] \
   && [ "$PHI_TOL" = "1e-6" ] \
   && [ "$SCHEMES" = "GaussLinear leastSquares" ] \
   && [ "$RESOLUTIONS" = "10 20 40 80" ]; then
  _PAPERI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../../.." && pwd)"
  python3 "$_PAPERI_ROOT/applications/scripts/paperI_results/aggregate.py" bidomain_tet \
      --repo-root "$_PAPERI_ROOT" \
      || echo "WARN: paperI aggregate (bidomain_tet) failed; native output untouched" >&2
else
  echo "NOTE: partial/alternative sweep -- canonical Paper I CSV left untouched." >&2
fi
