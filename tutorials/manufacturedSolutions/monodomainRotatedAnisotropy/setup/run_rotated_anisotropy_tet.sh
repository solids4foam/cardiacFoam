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

source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1

PY="${PYTHON:-python3}"
RESOLUTIONS="${RESOLUTIONS:-10 20 40 80}"
cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p setup/results

_FVSCH_BAK="$(mktemp)"
_CD_BAK="$(mktemp)"
cp system/fvSchemes "$_FVSCH_BAK"
cp system/controlDict "$_CD_BAK"
trap 'cp "$_FVSCH_BAK" system/fvSchemes; rm -f "$_FVSCH_BAK";
      cp "$_CD_BAK" system/controlDict; rm -f "$_CD_BAK"' EXIT

dt_for_n()
{
  case "$1" in
    10) echo 0.00892857;;
    20) echo 0.00224215;;
    40) echo 0.000560538;;
    80) echo 0.000140174;;
    *)  echo "no dt configured for N=$1" >&2; return 1;;
  esac
}

set_grad()
{
  if [ "$1" = "leastSquares" ]; then
    sed -E 's/^([[:space:]]*)default([[:space:]]+)Gauss linear;.*$/\1default\2leastSquares;/' \
      system/fvSchemes > system/fvSchemes.active
  else
    sed -E 's/^([[:space:]]*)default([[:space:]]+)leastSquares;.*$/\1default\2Gauss linear;/' \
      system/fvSchemes > system/fvSchemes.active
  fi
  mv system/fvSchemes.active system/fvSchemes
}

metrics()
{
  dat="$(ls postProcessing/rotatedAnisotropy_3D_*_cells_implicit.dat 2>/dev/null | head -1)"
  awk '
    /^Vm / { l1=$2; l2=$3; linf=$4 }
    /^Grid spacing/ { h=$NF }
    /^Time step/ { dt=$NF }
    /^Final simulation time/ { t=$NF }
    END { printf "%s %s %s %s %s %s", h, dt, t, l1, l2, linf }
  ' "$dat"
}

ecg_metrics()
{
  awk '
    /^[[:space:]]*#/ { next }
    NF >= 6 { e1=$2; e2=$3; e3=$4; e4=$5; e5=$6 }
    END { printf "%s %s %s %s %s", e1, e2, e3, e4, e5 }
  ' postProcessing/pseudoECG.dat 2>/dev/null
}

run_one()
{
  tag="$1"
  N="$2"

  ./Allclean >/dev/null 2>&1
  LC="$($PY -c "print(1.0/${N})")"
  sed "s|__LC__|${LC}|" setup/mesh/tet/box.geo.template > setup/mesh/tet/box.geo
  gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1
  gmshToFoam box.msh >/dev/null 2>&1
  rm -f box.msh setup/mesh/tet/box.geo

  DT="$(dt_for_n "$N")"
  sed -E "s/^deltaT.*/deltaT          ${DT};/; s/^endTime.*/endTime         0.2;/" \
    system/controlDict > system/controlDict.active
  mv system/controlDict.active system/controlDict

  runApplication decomposePar
  runParallel cardiacFoam
  runApplication reconstructPar

  M="$(metrics)"
  E="$(ecg_metrics)"
  echo "delaunayTet,${tag},${N},${M// /,},${E// /,}" \
    >> setup/results/rotated_anisotropy_tet_convergence.csv
  echo "done ${tag} N=${N}: Vm[${M}] ECG[${E}]"
}

echo "mesh,gradient,N,h,dt,t,Vm_L1,Vm_L2,Vm_Linf,ECG_E1,ECG_E2,ECG_E3,ECG_E4,ECG_E5" \
  > setup/results/rotated_anisotropy_tet_convergence.csv

set_grad "GaussLinear"
for N in $RESOLUTIONS; do
  run_one GaussLinear "$N"
done

set_grad "leastSquares"
for N in $RESOLUTIONS; do
  run_one leastSquares "$N"
done

./Allclean >/dev/null 2>&1
echo "=== rotated_anisotropy_tet_convergence.csv ==="
cat setup/results/rotated_anisotropy_tet_convergence.csv
