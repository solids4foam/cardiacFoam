#!/bin/bash
#SBATCH --job-name=cardiac-gpu-dt-summary
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=/tmp/cardiac_gpu_integrator_summary_%j.log

set -euo pipefail
cd /home/simao/cardiacFoam
python3 tutorials/ioniGPUKernel/slab2D/setup/post_processing_integrator_sweep_summary.py \
    /tmp/cardiac_gpu_integrator_sweep_20260929
