#!/bin/bash
#SBATCH --job-name=full-batched-backend
#SBATCH --partition=verylong
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=25-00:00:00
#SBATCH --output=/tmp/cardiac_full_backend_validation_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_full_backend_validation_${SLURM_JOB_ID}
mkdir -p "$OUT"
MODELS=(AlievPanfilov BuenoOrovio Courtemanche Gaur Grandi PerisYague Stewart TNNP TWorld ToRORd_dynCl Trovato)
RUN_STATUS=0

echo "Host=$(hostname) CPUs=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
nvidia-smi

unset CARDIAC_ENABLE_CUDA
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cpu.log" 2>&1
cd /home/simao/cardiacFoam
for geometry in slab niederer; do
    python3 tutorials/ioniGPUKernel/run_long_batched_backend_matrix.py \
        "$geometry" cpu "$OUT/cpu/$geometry" --models "${MODELS[@]}" \
        --delta-t 1e-5 --substeps 25 --end-time 0.08 --scheme sbdf2 \
        || RUN_STATUS=1
done

export CARDIAC_ENABLE_CUDA=1 CUDA_HOME=/usr NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cuda.log" 2>&1
cd /home/simao/cardiacFoam
for geometry in slab niederer; do
    python3 tutorials/ioniGPUKernel/run_long_batched_backend_matrix.py \
        "$geometry" cuda "$OUT/cuda/$geometry" --models "${MODELS[@]}" \
        --delta-t 1e-5 --substeps 25 --end-time 0.08 --scheme sbdf2 \
        || RUN_STATUS=1
done

mkdir -p "$OUT/comparisons/slab" "$OUT/comparisons/niederer"
for model in "${MODELS[@]}"; do
    if [[ -d "$OUT/cpu/slab/$model" && -d "$OUT/cuda/slab/$model" ]]; then
        python3 tutorials/ioniGPUKernel/compare_batched_waveforms.py \
            "$OUT/cpu/slab/$model" "$OUT/cuda/slab/$model" \
            --write-interval 1e-3 \
            > "$OUT/comparisons/slab/$model.txt" || RUN_STATUS=1
    fi
    if [[ -d "$OUT/cpu/niederer/$model" && -d "$OUT/cuda/niederer/$model" ]]; then
        python3 tutorials/ioniGPUKernel/compare_batched_waveforms.py \
            "$OUT/cpu/niederer/$model" "$OUT/cuda/niederer/$model" \
            --write-interval 1e-3 \
            > "$OUT/comparisons/niederer/$model.txt" || RUN_STATUS=1
    fi
done

echo "Completed full backend validation at $(date -Is); outputs=$OUT; status=$RUN_STATUS"
nvidia-smi
exit "$RUN_STATUS"
