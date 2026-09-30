#!/bin/bash
#SBATCH --job-name=mms-euler-substeps
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=/tmp/cardiac_manufactured_euler_substeps_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}
export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_manufactured_euler_substeps_20260929
if [ -e "$OUT" ]; then
    echo "Refusing to overwrite existing output: $OUT" >&2
    exit 2
fi

cd src/ionicModels
wclean libso
wmake libso > /tmp/cardiac_manufactured_euler_substeps_cuda_build.log 2>&1
cd /home/simao/cardiacFoam
python3 tutorials/ioniGPUKernel/run_manufactured_euler_substep_sweep.py \
    "$OUT" --require-gpu
