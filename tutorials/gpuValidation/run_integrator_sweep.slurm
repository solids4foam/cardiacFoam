#!/bin/bash
#SBATCH --job-name=cardiac-gpu-dt-sweep
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=/tmp/cardiac_gpu_integrator_sweep_%j.log

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
    /tmp/cardiac_gpu_integrator_sweep_20260928_nonfabbri \
    --models AlievPanfilov BuenoOrovio Courtemanche Gaur Grandi PerisYague \
        Stewart TNNP TWorld ToRORd_dynCl Trovato \
    --end-time 0.015 \
    --steps 1 5 10 25 \
    --dt-multipliers 1 2 5 \
    --integrators rushLarsen \
    --require-gpu
python3 tutorials/gpuValidation/summarize_integrator_study.py \
    /tmp/cardiac_gpu_integrator_sweep_20260928_nonfabbri
