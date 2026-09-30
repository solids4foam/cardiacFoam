#!/bin/bash
#SBATCH --job-name=equation-parity-fix
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=/tmp/cardiac_equation_alignment_singlecell_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_equation_alignment_singlecell_${SLURM_JOB_ID}
mkdir -p "$OUT"
MODELS=(BuenoOrovio Courtemanche PerisYague ToRORd_dynCl)

# First compare scalar against CPU-batched from a clean non-CUDA build.
unset CARDIAC_ENABLE_CUDA CUDA_HOME NVARCH CUDA_HOST_CXX
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cpu.log" 2>&1
cd /home/simao/cardiacFoam
python3 tutorials/gpuValidation/run_single_cell_smoke.py \
    --output-dir "$OUT/cpu" \
    --models "${MODELS[@]}" \
    --end-time 0.05 \
    --delta-t 2e-6 \
    --substeps 5 \
    --integrator rushLarsen \
    --all-variables \
    --write-frequency 1e-4

# Rebuild with CUDA and repeat the same scalar/batched comparison on device.
export CARDIAC_ENABLE_CUDA=1 CUDA_HOME=/usr NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
cd src/ionicModels
wclean libso
wmake libso > "$OUT/build_cuda.log" 2>&1
cd /home/simao/cardiacFoam
nvidia-smi
python3 tutorials/gpuValidation/run_single_cell_smoke.py \
    --output-dir "$OUT/cuda" \
    --models "${MODELS[@]}" \
    --end-time 0.05 \
    --delta-t 2e-6 \
    --substeps 5 \
    --integrator rushLarsen \
    --all-variables \
    --write-frequency 1e-4 \
    --require-gpu

for model in "${MODELS[@]}"; do
    grep -q 'using CUDA device' "$OUT/cuda/batched/$model/log.cardiacFoam"
    python3 tutorials/gpuValidation/compare_singlecell_rates.py \
        "$OUT/cpu/batched/$model" "$OUT/cuda/batched/$model" \
        | tee "$OUT/${model}_rates_cpu_vs_cuda.txt"
done
echo "Equation-alignment scalar/CPU-batched and scalar/CUDA-batched comparisons completed: $OUT"
nvidia-smi
