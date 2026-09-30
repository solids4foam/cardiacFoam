#!/bin/bash
#SBATCH --job-name=cardiac-niederer-tuned
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=/tmp/cardiac_niederer_tuned_heavy_%j.log

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
python3 tutorials/ioniGPUKernel/run_niederer_model_matrix.py \
    /tmp/cardiac_niederer_tuned_heavy_20260928/TNNP \
    --models TNNP --delta-t 1e-5 --substeps 1 --end-time 0.015
python3 tutorials/ioniGPUKernel/run_niederer_model_matrix.py \
    /tmp/cardiac_niederer_tuned_heavy_20260928/TWorld \
    --models TWorld --delta-t 1e-5 --substeps 5 --end-time 0.015
