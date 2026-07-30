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

RESOLUTIONS="${RESOLUTIONS:-10 20 40 80}"
cd "$(dirname "${BASH_SOURCE[0]}")/.."
mkdir -p setup/results

_CD_BAK="$(mktemp)"
cp system/controlDict "$_CD_BAK"
trap 'cp "$_CD_BAK" system/controlDict; rm -f "$_CD_BAK"' EXIT

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

metrics()
{
  dat="$(ls postProcessing/rotatedAnisotropy_3D_*_cells_implicit.dat 2>/dev/null | head -1)"
  awk '
    /^Vm / { l1=$2; l2=$3; linf=$4 }
    /^Time step/ { dt=$NF }
    /^Final simulation time/ { t=$NF }
    END { printf "%s %s %s %s %s", dt, t, l1, l2, linf }
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
  N="$1"
  ./Allclean >/dev/null 2>&1

  sed "s/NCELLS/${N}/g" system/blockMeshDict.template > system/blockMeshDict.active
  DT="$(dt_for_n "$N")"
  sed -E "s/^deltaT.*/deltaT          ${DT};/; s/^endTime.*/endTime         0.2;/" \
    system/controlDict > system/controlDict.active
  mv system/controlDict.active system/controlDict

  runApplication blockMesh -dict system/blockMeshDict.active
  rm -f system/blockMeshDict.active
  runApplication cardiacFoam

  M="$(metrics)"
  E="$(ecg_metrics)"
  H="$(awk "BEGIN { printf \"%.12g\", 1.0/${N} }")"

  echo "hex,default,${N},${H},${M// /,},${E// /,}" >> setup/results/rotated_anisotropy_hex_convergence.csv
  echo "done hex N=${N}: Vm[${M}] ECG[${E}]"
}

echo "mesh,gradient,N,h,dt,t,Vm_L1,Vm_L2,Vm_Linf,ECG_E1,ECG_E2,ECG_E3,ECG_E4,ECG_E5" \
  > setup/results/rotated_anisotropy_hex_convergence.csv

for N in $RESOLUTIONS; do
  run_one "$N"
done

./Allclean >/dev/null 2>&1
echo "=== rotated_anisotropy_hex_convergence.csv ==="
cat setup/results/rotated_anisotropy_hex_convergence.csv
