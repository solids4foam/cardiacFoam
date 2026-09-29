#!/bin/bash
#SBATCH --job-name=cardiac-sbdf2-mms
#SBATCH --partition=main
#SBATCH --nodelist=xenosim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --dependency=afterany:9358
#SBATCH --output=/tmp/cardiac_sbdf2_mms_%j.log

set -euo pipefail
set +e +u
source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null
set -e -u
unset CARDIAC_ENABLE_CUDA
export FOAM_SIGFPE=false
export OMP_NUM_THREADS=1

cd /home/simao/cardiacFoam
OUT=/tmp/cardiac_sbdf2_mms_temporal_20260928
if [ -e "$OUT" ]; then
    echo "Refusing to overwrite existing output: $OUT" >&2
    exit 2
fi
python3 tutorials/gpuValidation/run_sbdf2_manufactured_temporal.py "$OUT"
echo "SBDF2 manufactured temporal results: $OUT"
