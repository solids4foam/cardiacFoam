#!/bin/bash
#SBATCH --job-name=cardiac-heavy-model-dt
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=/tmp/cardiac_heavy_model_dt_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export NVARCH=75
export FOAM_SIGFPE=false
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
python3 tutorials/gpuValidation/run_slab_integrator_study.py \
    /tmp/cardiac_gpu_tnnp_tworld_coarse_dt_20260928 \
    --models TNNP TWorld \
    --end-time 0.015 \
    --steps 5 10 25 \
    --dt-multipliers 10 20 \
    --integrators rushLarsen \
    --require-gpu
python3 tutorials/gpuValidation/summarize_integrator_study.py \
    /tmp/cardiac_gpu_tnnp_tworld_coarse_dt_20260928
