#!/bin/bash
#SBATCH --job-name=court-torord-slab-parity
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=/tmp/cardiac_court_torord_slab_parity_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_court_torord_slab_parity_${SLURM_JOB_ID}
MODELS=(BuenoOrovio Courtemanche PerisYague ToRORd_dynCl)
SCHEMES=(godunov sbdf2)
mkdir -p "$OUT/cpu" "$OUT/gpu"

unset CARDIAC_ENABLE_CUDA
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cpu.log" 2>&1
cd /home/simao/cardiacFoam

for scheme in "${SCHEMES[@]}"; do
    for model in "${MODELS[@]}"; do
        case_dir="$OUT/cpu/$scheme/$model"
        python3 tutorials/gpuValidation/prepare_slab_case.py \
            "$model" batched "$case_dir" \
            --delta-t 2e-6 --substeps 5 --integrator rushLarsen \
            --time-coupling-scheme "$scheme" --end-time 0.015 \
            --write-precision 15
        (cd "$case_dir" && ./Allrun > allrun.stdout 2>&1)
        grep -q '^End$' "$case_dir/log.cardiacFoam"
        echo "CPU-batched $scheme completed: $model"
    done
done

export CARDIAC_ENABLE_CUDA=1 CUDA_HOME=/usr NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cuda.log" 2>&1
cd /home/simao/cardiacFoam

for scheme in "${SCHEMES[@]}"; do
    for model in "${MODELS[@]}"; do
        case_dir="$OUT/gpu/$scheme/$model"
        python3 tutorials/gpuValidation/prepare_slab_case.py \
            "$model" batched "$case_dir" \
            --delta-t 2e-6 --substeps 5 --integrator rushLarsen \
            --time-coupling-scheme "$scheme" --end-time 0.015 \
            --write-precision 15
        (cd "$case_dir" && ./Allrun > allrun.stdout 2>&1)
        grep -q '^End$' "$case_dir/log.cardiacFoam"
        grep -q 'using CUDA device' "$case_dir/log.cardiacFoam"
        echo "CUDA-batched $scheme completed: $model"
    done
done

STATUS=0
for scheme in "${SCHEMES[@]}"; do
    for model in "${MODELS[@]}"; do
        report="$OUT/${scheme}_${model}_cpu_vs_cuda.txt"
        if python3 tutorials/gpuValidation/compare_sbdf2_batched_backends.py \
            "$model" "$OUT/cpu/$scheme/$model" "$OUT/gpu/$scheme/$model" \
            --scheme "$scheme" --times 0.005 0.01 0.015 > "$report"; then
            echo "$model $scheme CPU/CUDA slab parity: PASS"
        else
            echo "$model $scheme CPU/CUDA slab parity: FAIL"
            STATUS=1
        fi
    done
done

echo "Courtemanche/ToRORd slab backend parity results: $OUT"
exit "$STATUS"
