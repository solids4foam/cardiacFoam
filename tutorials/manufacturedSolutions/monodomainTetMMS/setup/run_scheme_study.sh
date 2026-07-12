#!/bin/bash
set +e
source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
PY=/Users/simaocastro/noFrontendCardiacFoam_minor_errors/.venv/bin/python3
cd "$(dirname "${BASH_SOURCE[0]}")/.."
dt_for_n(){ case "$1" in 10) echo 0.00892857;;20) echo 0.00224215;;40) echo 0.000560538;;80) echo 0.000140174;;esac; }
set_grad(){
  if [ "$1" = "leastSquares" ]; then
    sed -E 's/^([[:space:]]*)default([[:space:]]+)Gauss linear;[[:space:]]*$/\1default\2leastSquares;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  else
    sed -E 's/^([[:space:]]*)default([[:space:]]+)leastSquares.*$/\1default\2Gauss linear;/' system/fvSchemes > /tmp/fvs && mv /tmp/fvs system/fvSchemes
  fi
}
mono_metrics(){ D=$(ls postProcessing/3D_*_cells_implicit.dat 2>/dev/null|head -1); awk '/Grid spacing/{dx=$NF} /^Vm /{l2=$3;li=$4} END{printf "%s %s %s",dx,l2,li}' "$D"; }
ecg_metrics(){ awk 'NR>=7{if($3+0>a)a=$3; if($4+0>b)b=$4} END{printf "%s %s",a,b}' postProcessing/manufacturedPseudoECGSummary.dat; }
run_one(){ # $1 tag $2 N
  ./Allclean >/dev/null 2>&1
  LC=$($PY -c "print(1.0/$2)"); sed "s|__LC__|$LC|" setup/box.geo.template > setup/box.geo
  gmsh -3 setup/box.geo -o box.msh -format msh2 >/dev/null 2>&1; gmshToFoam box.msh >/dev/null 2>&1; rm -f box.msh
  DT=$(dt_for_n $2); sed -E "s/^deltaT.*/deltaT    $DT;/; s/^endTime.*/endTime    0.2;/" system/controlDict>/tmp/cd && mv /tmp/cd system/controlDict
  cardiacFoam > log.cf 2>&1
  M=$(mono_metrics); E=$(ecg_metrics)
  echo "$1,$2,${M// /,},${E// /,}" >> setup/results/scheme_study.csv
  echo "done $1 N=$2 : mono[$M] ecg[$E]"
}
echo "scheme,N,dx,mono_L2,mono_Linf,ecg_L2,ecg_Linf" > setup/results/scheme_study.csv
set_grad "GaussLinear"; for N in 10 20 40; do run_one GaussLinear $N; done
set_grad "leastSquares"; for N in 10 20 40 80; do run_one leastSquares $N; done
set_grad "leastSquares"   # restore tutorial default
./Allclean >/dev/null 2>&1
echo "=== scheme_study.csv ==="; cat setup/results/scheme_study.csv
