#!/bin/bash -l
#SBATCH --job-name=manufactured-cuda-smoke
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=dev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00

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

cd /home/simao/cardiacFoam
echo "Running the scalar and required-CUDA manufactured cases on $(hostname)"
python3 tutorials/ioniGPUKernel/run_manufactured_cuda_smoke.py \
    /tmp/cardiac_manufactured_cuda_smoke_20260928 --require-gpu
