#!/bin/bash
#SBATCH --job-name=sbdf2-all-model-parity
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=/tmp/cardiac_sbdf2_all_model_parity_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_sbdf2_all_model_parity_${SLURM_JOB_ID}
mkdir -p "$OUT/cpu" "$OUT/gpu"
MODELS=(AlievPanfilov BuenoOrovio Courtemanche Gaur Grandi PerisYague Stewart TNNP TWorld ToRORd_dynCl Trovato)

unset CARDIAC_ENABLE_CUDA
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cpu.log" 2>&1
cd /home/simao/cardiacFoam

for model in "${MODELS[@]}"; do
    case_dir="$OUT/cpu/$model"
    python3 tutorials/gpuValidation/prepare_slab_case.py \
        "$model" batched "$case_dir" \
        --delta-t 2e-6 --substeps 5 --integrator rushLarsen \
        --time-coupling-scheme sbdf2 --end-time 0.015 \
        --write-precision 15
    (cd "$case_dir" && ./Allrun > allrun.stdout 2>&1)
    grep -q '^End$' "$case_dir/log.cardiacFoam"
    echo "CPU-batched SBDF2 completed: $model"
done

export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cuda.log" 2>&1
cd /home/simao/cardiacFoam

for model in "${MODELS[@]}"; do
    case_dir="$OUT/gpu/$model"
    python3 tutorials/gpuValidation/prepare_slab_case.py \
        "$model" batched "$case_dir" \
        --delta-t 2e-6 --substeps 5 --integrator rushLarsen \
        --time-coupling-scheme sbdf2 --end-time 0.015 \
        --write-precision 15
    (cd "$case_dir" && ./Allrun > allrun.stdout 2>&1)
    grep -q '^End$' "$case_dir/log.cardiacFoam"
    grep -q 'using CUDA device' "$case_dir/log.cardiacFoam"
    echo "CUDA-batched SBDF2 completed: $model"
done

PARITY_STATUS=0
for model in "${MODELS[@]}"; do
    mkdir -p "$OUT/$model"
    if python3 tutorials/gpuValidation/compare_sbdf2_batched_backends.py \
        "$model" "$OUT/cpu/$model" "$OUT/gpu/$model" \
        --times 0.005 0.01 0.015 \
        > "$OUT/$model/cpu_batched_vs_cuda_batched.txt"; then
        echo "$model SBDF2 CPU/CUDA parity: PASS"
    else
        echo "$model SBDF2 CPU/CUDA parity: FAIL"
        PARITY_STATUS=1
    fi
done

echo "SBDF2 all-model CPU/CUDA parity results: $OUT"
exit "$PARITY_STATUS"
