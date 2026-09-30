#!/bin/bash
#SBATCH --job-name=fabbri-cell-parity
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=/tmp/cardiac_fabbri_singlecell_parity_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_fabbri_singlecell_parity_${SLURM_JOB_ID}
mkdir -p "$OUT"

# Build CPU-only and retain scalar plus CPU-batched single-cell traces.
unset CARDIAC_ENABLE_CUDA CUDA_HOME NVARCH CUDA_HOST_CXX
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cpu.log" 2>&1
cd /home/simao/cardiacFoam
python3 tutorials/ioniGPUKernel/run_single_cell_smoke.py \
    --output-dir "$OUT/cpu" \
    --models Fabbri \
    --end-time 2.0 \
    --delta-t 2e-6 \
    --substeps 5 \
    --integrator rushLarsen \
    --all-variables \
    --write-frequency 1e-4

# Rebuild for CUDA and require a real device for the batched case.
export CARDIAC_ENABLE_CUDA=1 CUDA_HOME=/usr NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cuda.log" 2>&1
cd /home/simao/cardiacFoam
nvidia-smi
python3 tutorials/ioniGPUKernel/run_single_cell_smoke.py \
    --output-dir "$OUT/cuda" \
    --models Fabbri \
    --end-time 2.0 \
    --delta-t 2e-6 \
    --substeps 5 \
    --integrator rushLarsen \
    --all-variables \
    --write-frequency 1e-4 \
    --require-gpu

grep -q 'using CUDA device' "$OUT/cuda/batched/Fabbri/log.cardiacFoam"
python3 tutorials/ioniGPUKernel/compare_singlecell_rates.py \
    "$OUT/cpu/batched/Fabbri" "$OUT/cuda/batched/Fabbri" \
    --start-time 1.0 \
    | tee "$OUT/Fabbri_cpu_vs_cuda.txt"
echo "Fabbri single-cell CPU/CUDA backend parity over 1-2 s completed: $OUT"
nvidia-smi
