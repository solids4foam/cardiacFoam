#!/bin/bash
#SBATCH --job-name=manufactured-batch-parity
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --dependency=afterany:9359
#SBATCH --output=/tmp/cardiac_manufactured_batched_temporal_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_manufactured_batched_temporal_20260928
if [ -e "$OUT" ]; then
    echo "Refusing to overwrite existing output: $OUT" >&2
    exit 2
fi
mkdir -p "$OUT"

# First exercise the registered CPU-batched implementation.
unset CARDIAC_ENABLE_CUDA
cd src/ionicModels
wclean libso
wmake libso > /tmp/cardiac_manufactured_batched_cpu_build.log 2>&1
cd /home/simao/cardiacFoam
set +e
python3 tutorials/ioniGPUKernel/run_manufactured_batched_temporal.py \
    "$OUT/cpu" --backend cpu --substeps 25
CPU_STATUS=$?
set -e

# Rebuild cleanly and require actual CUDA device selection for every case.
export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
cd src/ionicModels
wclean libso
wmake libso > /tmp/cardiac_manufactured_batched_cuda_build.log 2>&1
cd /home/simao/cardiacFoam
set +e
python3 tutorials/ioniGPUKernel/run_manufactured_batched_temporal.py \
    "$OUT/gpu" --backend gpu --substeps 25 --require-gpu
GPU_STATUS=$?
set -e

set +e
python3 tutorials/ioniGPUKernel/compare_manufactured_batched_parity.py \
    "$OUT/cpu" "$OUT/gpu" > "$OUT/cpu_vs_gpu_parity.txt" 2>&1
PARITY_STATUS=$?
set -e
cat "$OUT/cpu_vs_gpu_parity.txt"
echo "CPU temporal-ladder status: $CPU_STATUS"
echo "CUDA temporal-ladder status: $GPU_STATUS"
echo "CPU/GPU parity status: $PARITY_STATUS"
test "$CPU_STATUS" -eq 0
test "$GPU_STATUS" -eq 0
test "$PARITY_STATUS" -eq 0
