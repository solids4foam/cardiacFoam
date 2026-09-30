#!/bin/bash
#SBATCH --job-name=cardiac-sbdf2-baseline
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=/tmp/cardiac_sbdf2_baseline_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
export FOAM_SIGFPE=false
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam/src/ionicModels
wmake libso > /tmp/cardiac_sbdf2_ionic_build.log 2>&1

cd /home/simao/cardiacFoam
python3 tutorials/ioniGPUKernel/run_sbdf2_baseline.py \
    /tmp/cardiac_gpu_sbdf2_baseline_20260928 \
    --models TNNP TWorld \
    --delta-t 2e-6 \
    --substeps 5 \
    --end-time 0.015 \
    --require-gpu
