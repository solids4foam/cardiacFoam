#!/bin/bash
#SBATCH --job-name=cardiac-sbdf2-parity
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --dependency=afterany:9356
#SBATCH --output=/tmp/cardiac_sbdf2_parity_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_gpu_sbdf2_batched_parity_20260928
if [ -e "$OUT" ]; then
    echo "Refusing to overwrite existing output: $OUT" >&2
    exit 2
fi

# Build without CUDA to exercise the registered batched host implementation.
unset CARDIAC_ENABLE_CUDA
cd src/ionicModels
wclean libso
wmake libso > /tmp/cardiac_sbdf2_parity_cpu_build.log 2>&1
cd /home/simao/cardiacFoam
python3 tutorials/gpuValidation/run_sbdf2_baseline.py "$OUT/cpu" \
    --models TNNP TWorld --delta-t 2e-6 --substeps 5 --end-time 0.015 \
    --write-precision 15

# Rebuild from clean objects with CUDA, then run the identical cases on-device.
export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
cd src/ionicModels
wclean libso
wmake libso > /tmp/cardiac_sbdf2_parity_cuda_build.log 2>&1
cd /home/simao/cardiacFoam
python3 tutorials/gpuValidation/run_sbdf2_baseline.py "$OUT/gpu" \
    --models TNNP TWorld --delta-t 2e-6 --substeps 5 --end-time 0.015 \
    --require-gpu --write-precision 15

PARITY_STATUS=0
for model in TNNP TWorld; do
    grep -Eq "^${model},sbdf2,batched,completed-CPU," "$OUT/cpu/summary.csv"
    grep -Eq "^${model},sbdf2,batched,completed-GPU," "$OUT/gpu/summary.csv"
    mkdir -p "$OUT/$model"
    if python3 tutorials/gpuValidation/compare_batched_parity.py "$model" \
        "$OUT/cpu/$model/sbdf2_batched" \
        "$OUT/gpu/$model/sbdf2_batched" \
        --times 0.005 0.01 0.015 \
        > "$OUT/$model/cpu_batched_vs_gpu_batched_sbdf2.txt"; then
        echo "$model CPU/GPU batched parity: PASS"
    else
        echo "$model CPU/GPU batched parity: FAIL; see $OUT/$model/cpu_batched_vs_gpu_batched_sbdf2.txt"
        PARITY_STATUS=1
    fi
done

echo "SBDF2 CPU-batched/GPU-batched parity results: $OUT"
exit "$PARITY_STATUS"
